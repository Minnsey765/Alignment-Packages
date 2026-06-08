# _phylo_verify.py

import os
import json
import random
import subprocess
import requests
import time
import numpy as np

from pathlib import Path
from Bio import SeqIO, AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from ..verify_funcs._score_hit import (recommended_subsample_length,
                                            recommended_n_samples)
    from ..verify_funcs._species_parser import load_taxonomy
except ImportError:
    from _score_hit import recommended_subsample_length, recommended_n_samples
    from _species_parser import load_taxonomy

try:
    from ..alignment_funcs._clustalOmega import (run_clustalo,
                                                  build_short_ids,
                                                  restore_headers_in_file)
except ImportError:
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    '..', 'alignment_funcs'))
    from _clustalOmega import (run_clustalo, build_short_ids,
                                restore_headers_in_file)

BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"


# ─────────────────────────────────────────────────────────────────────────────
# TAXONOMY HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def get_expected_order(taxon: str, taxonomy: dict) -> str:
    """
    Return the taxonomic order for a taxon string of the form
    Genus_species, looked up from the taxonomy dict.
    """
    genus    = taxon.split("_")[0]
    tax_info = taxonomy.get(genus)
    if tax_info:
        return tax_info.get("order")
    return None


def get_ingroup_genera(order: str, taxonomy: dict) -> set:
    """
    Return the set of all genera belonging to the given order
    according to the taxonomy dict.
    """
    return {
        genus for genus, info in taxonomy.items()
        if info.get("order") == order
    }


# ─────────────────────────────────────────────────────────────────────────────
# ALIGNMENT TRIMMING
# ─────────────────────────────────────────────────────────────────────────────

def trim_to_shared_region(fasta_alignment: str,
                           output_path: str,
                           min_coverage: float = 0.5) -> dict:
    """
    Trim an alignment by removing leading and lagging columns where
    coverage falls below min_coverage. Internal columns are kept
    regardless of their gap content, preserving genuine indels that
    may be phylogenetically informative.

    Only the ends of the alignment are trimmed — this corrects for
    partial sequences that do not cover the full gene length without
    penalising real insertion/deletion events within the alignment.

    Sequences that are entirely gaps or Ns after trimming are dropped.

    Parameters
    ----------
    fasta_alignment : str    path to aligned FASTA file
    output_path     : str    path to write trimmed FASTA output
    min_coverage    : float  minimum fraction of non-gap characters
                             required to keep a terminal column.
                             Default 0.5 — appropriate for small
                             datasets (~12 taxa) with distant outgroups.
                             Use 0.8 for large datasets with closely
                             related sequences only.

    Returns
    -------
    dict with:
        original_length  : int   alignment length before trimming
        trimmed_length   : int   alignment length after trimming
        trim_start       : int   first kept column index
        trim_end         : int   last kept column index (exclusive)
        columns_removed  : int   total columns removed from ends
        sequences_kept   : int   number of sequences written
        sequences_dropped: list  IDs of sequences dropped (all gaps)
    """
    alignment = AlignIO.read(fasta_alignment, "fasta")
    n_seq     = len(alignment)
    aln_len   = alignment.get_alignment_length()

    # Calculate per-column coverage (fraction of non-gap characters)
    coverage = np.array([
        1.0 - (alignment[:, i].count("-") / n_seq)
        for i in range(aln_len)
    ])

    covered = coverage >= min_coverage

    if not covered.any():
        print(f"  WARNING: No columns meet min_coverage {min_coverage} "
              f"— returning untrimmed alignment")
        first = 0
        last  = aln_len
    else:
        # Find the first and last column that meets the threshold —
        # everything between these two positions is kept, preserving
        # internal gaps from genuine insertions/deletions
        first = int(np.argmax(covered))
        last  = int(aln_len - np.argmax(covered[::-1]))

    keep_cols = list(range(first, last))

    n_removed = aln_len - len(keep_cols)
    print(f"  End-trimming: {aln_len}bp → {len(keep_cols)}bp "
          f"(removed {first} leading + {aln_len - last} trailing columns, "
          f"{len(keep_cols)/aln_len*100:.1f}% retained)")

    if len(keep_cols) < 50:
        print(f"  WARNING: Only {len(keep_cols)}bp retained after "
              f"trimming. Consider lowering min_coverage or checking "
              f"that scaffold sequences overlap sufficiently.")

    trimmed_records = []
    dropped         = []

    for record in alignment:
        trimmed_seq = "".join(str(record.seq)[i] for i in keep_cols)
        if set(trimmed_seq) <= {"-", "N", "n"}:
            dropped.append(record.id)
            print(f"  Dropping {record.id} — "
                  f"no real sequence in shared region")
            continue
        trimmed_records.append(
            SeqRecord(Seq(trimmed_seq), id=record.id, description="")
        )

    SeqIO.write(trimmed_records, output_path, "fasta")

    return {
        "original_length":   aln_len,
        "trimmed_length":    len(keep_cols),
        "trim_start":        first,
        "trim_end":          last,
        "columns_removed":   n_removed,
        "sequences_kept":    len(trimmed_records),
        "sequences_dropped": dropped,
    }


# ─────────────────────────────────────────────────────────────────────────────
# CONSTRAINT TREE
# ─────────────────────────────────────────────────────────────────────────────

def write_constraint_tree(taxonomy: dict,
                           scaffold_fasta_aligned: str,
                           output_path: str) -> str:
    """
    Write a Newick constraint tree for IQ-TREE that enforces
    order-level monophyly based on the taxonomy dict. Each order
    present in the scaffold alignment forms a constrained clade.
    Internal relationships within orders are left free.

    This prevents long-branch attraction artefacts (e.g. Mus grouping
    with fast-evolving soricids) while still allowing branch lengths
    and within-order topology to be estimated from the data.

    Parameters
    ----------
    taxonomy               : dict  from load_taxonomy()
    scaffold_fasta_aligned : str   aligned scaffold FASTA — used to
                                   get the taxon IDs actually present
    output_path            : str   path to write constraint Newick

    Returns
    -------
    str   path to the written constraint tree file
    """
    from collections import defaultdict

    records = list(SeqIO.parse(scaffold_fasta_aligned, "fasta"))
    taxa    = [r.id for r in records]

    order_groups = defaultdict(list)
    unplaced     = []

    for taxon_id in taxa:
        # Handle both pipe-delimited and plain headers
        genus    = taxon_id.split("|")[0].split("_")[0]
        tax_info = taxonomy.get(genus)
        if tax_info:
            order_groups[tax_info["order"]].append(taxon_id)
        else:
            unplaced.append(taxon_id)
            print(f"  Warning: '{taxon_id}' not found in taxonomy — "
                  f"left unconstrained")

    print(f"\nConstraint tree order groups:")
    for order, members in order_groups.items():
        print(f"  {order}: {len(members)} taxa")
    if unplaced:
        print(f"  Unconstrained: {unplaced}")

    # Each order forms a monophyletic clade
    # Inter-order relationships are unresolved (polytomy at root)
    order_clades = []
    for order, members in order_groups.items():
        if len(members) == 1:
            order_clades.append(members[0])
        else:
            inner = ",".join(members)
            order_clades.append(f"({inner})")

    all_clades        = order_clades + unplaced
    constraint_newick = ("(" + ",".join(all_clades) + ");"
                         if len(all_clades) > 1
                         else all_clades[0] + ";")

    with open(output_path, "w") as f:
        f.write(constraint_newick)

    print(f"  Constraint tree written: {output_path}")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# IQ-TREE
# ─────────────────────────────────────────────────────────────────────────────

def run_iqtree(nexus_alignment: str,
               output_dir: str,
               iqtree_bin: str = "iqtree2",
               model: str = "GTR+G",
               bootstrap: int = 1000,
               n_threads: int = 4,
               seed: int = 12345,
               constraint_tree: str = None) -> dict:
    """
    Run IQ-TREE on a Nexus alignment with a fixed random seed for
    reproducibility.

    Parameters
    ----------
    nexus_alignment : str   path to Nexus alignment. Must use short
                            sequence IDs — IQ-TREE fails on long IDs
                            or IDs containing colons.
    output_dir      : str   directory for IQ-TREE output files
    iqtree_bin      : str   path to iqtree2/iqtree3 binary
    model           : str   substitution model. Use "GTR+G" for speed
                            or "TEST" for automatic model selection.
    bootstrap       : int   ultrafast bootstrap replicates (min 1000)
    n_threads       : int   CPU threads
    seed            : int   random seed — fix this for reproducibility.
                            Without a fixed seed the tree topology may
                            change between runs due to the heuristic
                            search algorithm.
    constraint_tree : str   optional path to a Newick constraint tree.
                            If provided, IQ-TREE only searches topologies
                            consistent with the constraint (-g flag).
                            Use write_constraint_tree() to generate one.

    Returns
    -------
    dict with:
        treefile : str   path to .treefile (Newick with bootstrap)
        log      : str   path to IQ-TREE log file
    """
    os.makedirs(output_dir, exist_ok=True)

    prefix = os.path.join(
        output_dir,
        os.path.splitext(os.path.basename(nexus_alignment))[0]
    )

    cmd = [
        iqtree_bin,
        "-s",        nexus_alignment,
        "-m",        model,
        "-B",        str(bootstrap),
        "-T",        str(n_threads),
        "--prefix",  prefix,
        "--seed",    str(seed),
        "--redo",
    ]

    if constraint_tree:
        cmd += ["-g", constraint_tree]
        print(f"  Using constraint tree: {constraint_tree}")

    print(f"\nRunning IQ-TREE...")
    print(f"  Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"IQ-TREE failed:\n{result.stderr}\n{result.stdout}"
        )

    treefile = f"{prefix}.treefile"
    logfile  = f"{prefix}.log"

    if not os.path.exists(treefile):
        raise FileNotFoundError(
            f"IQ-TREE finished but treefile not found: {treefile}"
        )

    print(f"  ✓ IQ-TREE complete")
    print(f"  Tree: {treefile}")

    return {
        "treefile": treefile,
        "log":      logfile,
    }


# ─────────────────────────────────────────────────────────────────────────────
# SCAFFOLD TREE  (run once per gene, reused for all queries)
# ─────────────────────────────────────────────────────────────────────────────

def build_scaffold_tree(scaffold_nexus: str,
                         scaffold_tree_dir: str,
                         taxonomy: dict = None,
                         scaffold_fasta_aligned: str = None,
                         iqtree_bin: str = "iqtree2",
                         model: str = "GTR+G",
                         bootstrap: int = 1000,
                         n_threads: int = 4,
                         seed: int = 12345) -> dict:
    """
    Build a phylogenetic tree from an already-aligned scaffold Nexus
    file. Run this once per gene — the output is reused for all
    subsequent query verifications via verify_by_phylogeny().

    If taxonomy and scaffold_fasta_aligned are both provided, a
    constraint tree enforcing order-level monophyly is automatically
    generated and applied. This prevents long-branch attraction
    artefacts (e.g. fast-evolving soricids grouping with Mus) and is
    strongly recommended when the scaffold contains distant outgroups.

    Parameters
    ----------
    scaffold_nexus         : str   path to aligned scaffold Nexus.
                                   Must be the short-ID version
                                   (_short_aln.nex) to avoid IQ-TREE
                                   parsing errors from colons in IDs.
    scaffold_tree_dir      : str   output directory for scaffold tree.
                                   Keep this separate from query output
                                   directories.
    taxonomy               : dict  from load_taxonomy(). Required if
                                   using constraint tree.
    scaffold_fasta_aligned : str   path to full-header aligned FASTA
                                   (_aln.fasta). Required if using
                                   constraint tree — used to read taxon
                                   IDs for the constraint.
    iqtree_bin             : str   IQ-TREE binary path
    model                  : str   substitution model
    bootstrap              : int   ultrafast bootstrap replicates
    n_threads              : int   CPU threads
    seed                   : int   random seed

    Returns
    -------
    dict with:
        treefile         : str   path to scaffold .treefile
        log              : str   path to IQ-TREE log
        scaffold_nexus   : str   input Nexus path
        scaffold_tree_dir: str   output directory
        constraint_tree  : str   path to constraint Newick (or None)
    """
    os.makedirs(scaffold_tree_dir, exist_ok=True)

    print(f"\nBuilding scaffold tree from: {scaffold_nexus}")
    print(f"Output directory:            {scaffold_tree_dir}")

    # Generate constraint tree if taxonomy provided
    constraint_path = None
    if taxonomy is not None and scaffold_fasta_aligned is not None:
        constraint_path = os.path.join(scaffold_tree_dir, "constraint.nwk")
        write_constraint_tree(
            taxonomy               = taxonomy,
            scaffold_fasta_aligned = scaffold_fasta_aligned,
            output_path            = constraint_path,
        )
    elif taxonomy is not None or scaffold_fasta_aligned is not None:
        print("  Warning: both taxonomy and scaffold_fasta_aligned must "
              "be provided to use a constraint tree. Skipping constraint.")

    iqtree_result = run_iqtree(
        nexus_alignment = scaffold_nexus,
        output_dir      = scaffold_tree_dir,
        iqtree_bin      = iqtree_bin,
        model           = model,
        bootstrap       = bootstrap,
        n_threads       = n_threads,
        seed            = seed,
        constraint_tree = constraint_path,
    )

    return {
        **iqtree_result,
        "scaffold_nexus":    scaffold_nexus,
        "scaffold_tree_dir": scaffold_tree_dir,
        "constraint_tree":   constraint_path,
    }


# ─────────────────────────────────────────────────────────────────────────────
# PROFILE ALIGNMENT  (add query subsamples to scaffold alignment)
# ─────────────────────────────────────────────────────────────────────────────

def add_queries_to_alignment(scaffold_fasta_aligned: str,
                              query_records: list,
                              output_dir: str,
                              query_taxon: str,
                              query_accession: str,
                              email: str = "om380@cam.ac.uk",
                              max_polls: int = 120,
                              poll_interval: int = 10,
                              min_coverage: float = 0.5) -> dict:
    """
    Add query subsamples to an existing scaffold alignment using
    ClustalOmega profile alignment. The scaffold alignment is not
    re-aligned — query sequences are simply placed into it, preserving
    the scaffold topology and alignment quality.

    The combined alignment is then end-trimmed to remove leading and
    lagging columns where coverage falls below min_coverage. Internal
    gap columns are preserved.

    Parameters
    ----------
    scaffold_fasta_aligned : str    path to aligned scaffold FASTA
                                    (full-header _aln.fasta, not short)
    query_records          : list   SeqRecord objects for query subsamples
    output_dir             : str    directory to save output files
    query_taxon            : str    e.g. "Solenodon_paradoxus"
    query_accession        : str    e.g. "AY530070"
    email                  : str    EBI ClustalOmega API email
    max_polls              : int    maximum status poll attempts
    poll_interval          : int    seconds between status polls
    min_coverage           : float  end-trimming coverage threshold.
                                    Default 0.5 for small diverse
                                    scaffolds.

    Returns
    -------
    dict with:
        combined_fasta       : str   full-header combined aligned FASTA
        combined_nexus       : str   full-header combined aligned Nexus
        short_combined_nexus : str   short-ID combined aligned Nexus
        trimmed_fasta        : str   end-trimmed full-header FASTA
        trimmed_nexus_short  : str   end-trimmed short-ID Nexus
                                     (pass this to run_iqtree)
        id_map               : dict  short_id -> full_id
        id_map_path          : str   path to saved id_map JSON
        trim_stats           : dict  trimming statistics
        n_scaffold           : int   number of scaffold sequences
        n_queries            : int   number of query subsamples
    """
    os.makedirs(output_dir, exist_ok=True)
    base = f"{query_taxon}_{query_accession}"

    # ── Write query subsamples to FASTA ───────────────────────────────────────
    query_fasta = os.path.join(output_dir, f"{base}_queries.fasta")
    SeqIO.write(query_records, query_fasta, "fasta")

    # ── Build short IDs for both scaffold and query sequences ─────────────────
    # ClustalOmega and IQ-TREE both struggle with long or
    # colon-containing IDs, so we remap everything before submitting
    scaffold_records = list(SeqIO.parse(scaffold_fasta_aligned, "fasta"))
    all_records      = scaffold_records + query_records
    n_scaffold       = len(scaffold_records)

    short_all, id_map, _ = build_short_ids(all_records)
    short_scaffold        = short_all[:n_scaffold]
    short_queries         = short_all[n_scaffold:]

    short_scaffold_path = os.path.join(output_dir,
                                        f"{base}_short_scaffold.fasta")
    short_query_path    = os.path.join(output_dir,
                                        f"{base}_short_queries.fasta")

    SeqIO.write(short_scaffold, short_scaffold_path, "fasta")
    SeqIO.write(short_queries,  short_query_path,    "fasta")

    map_path = os.path.join(output_dir, f"{base}_id_map.json")
    with open(map_path, "w") as f:
        json.dump(id_map, f, indent=2)

    # ── Submit profile alignment to EBI ClustalOmega ──────────────────────────
    print(f"\nSubmitting profile alignment to EBI ClustalOmega...")
    print(f"  Scaffold sequences : {len(short_scaffold)}")
    print(f"  Query subsamples   : {len(short_queries)}")

    with open(short_scaffold_path, "r") as f:
        scaffold_data = f.read()
    with open(short_query_path, "r") as f:
        query_data = f.read()

    results = {}

    for fmt, result_key, suffix in [
        ("fasta", "aln-fasta", "_combined_aln.fasta"),
        ("nexus", "aln-nexus", "_combined_aln.nex"),
    ]:
        params = {
            "sequence":  scaffold_data,
            "sequence2": query_data,
            "stype":     "dna",
            "email":     email,
            "outfmt":    fmt,
            "order":     "aligned",
        }

        response = requests.post(BASE_URL + "/run", data=params)
        response.raise_for_status()
        job_id = response.text.strip()
        print(f"  Job submitted ({fmt}). ID: {job_id}")

        # Poll until finished
        status = "RUNNING"
        polls  = 0
        while status in ("RUNNING", "PENDING") and polls < max_polls:
            time.sleep(poll_interval)
            status = requests.get(
                BASE_URL + f"/status/{job_id}"
            ).text.strip()
            polls += 1
            print(f"  [{polls}] {status}")

        if status != "FINISHED":
            raise RuntimeError(
                f"ClustalOmega job {job_id} ended with status: {status}"
            )

        print(f"  Fetching {fmt} result...")
        short_out = os.path.join(output_dir, f"{base}_short{suffix}")
        full_out  = os.path.join(output_dir, f"{base}{suffix}")

        try:
            res = requests.get(BASE_URL + f"/result/{job_id}/{result_key}")
            res.raise_for_status()
            with open(short_out, "w") as f:
                f.write(res.text)
        except requests.exceptions.HTTPError as e:
            if fmt == "nexus":
                raise
            # Fasta fetch failed — convert from nexus locally
            print(f"  FASTA fetch failed ({e}) — converting locally")
            aln = AlignIO.read(results["nexus"]["short"], "nexus")
            AlignIO.write(aln, short_out, "fasta")

        restore_headers_in_file(short_out, id_map, full_out)
        print(f"  Saved: {full_out}")
        results[fmt] = {"short": short_out, "full": full_out}

    # ── End-trim to shared region ─────────────────────────────────────────────
    trimmed_fasta = os.path.join(output_dir, f"{base}_trimmed.fasta")
    trim_stats    = trim_to_shared_region(
        fasta_alignment = results["fasta"]["full"],
        output_path     = trimmed_fasta,
        min_coverage    = min_coverage,
    )

    # ── Convert trimmed FASTA to short-ID Nexus for IQ-TREE ──────────────────
    trimmed_nexus_short = trimmed_fasta.replace(".fasta", "_short.nex")
    short_trimmed_fasta = os.path.join(output_dir,
                                        f"{base}_trimmed_short.fasta")

    reverse_map     = {v: k for k, v in id_map.items()}
    trimmed_records = list(SeqIO.parse(trimmed_fasta, "fasta"))

    short_trimmed = []
    for r in trimmed_records:
        short_id = reverse_map.get(r.id, r.id)
        short_trimmed.append(
            SeqRecord(r.seq, id=short_id, description="")
        )

    SeqIO.write(short_trimmed, short_trimmed_fasta, "fasta")
    AlignIO.write(
        AlignIO.read(short_trimmed_fasta, "fasta"),
        trimmed_nexus_short,
        "nexus"
    )

    return {
        "combined_fasta":       results["fasta"]["full"],
        "combined_nexus":       results["nexus"]["full"],
        "short_combined_nexus": results["nexus"]["short"],
        "trimmed_fasta":        trimmed_fasta,
        "trimmed_nexus_short":  trimmed_nexus_short,
        "id_map":               id_map,
        "id_map_path":          map_path,
        "trim_stats":           trim_stats,
        "n_scaffold":           n_scaffold,
        "n_queries":            len(query_records),
    }


# ─────────────────────────────────────────────────────────────────────────────
# BUILD QUERY SUBSAMPLES
# ─────────────────────────────────────────────────────────────────────────────

def build_query_subsamples(query_seq: str,
                            query_taxon: str,
                            query_gene: str,
                            query_accession: str,
                            query_orientation: str,
                            insert_fraction: float = 0.30,
                            target_detection_prob: float = 0.90) -> dict:
    """
    Draw n_samples random subsamples from a query sequence, where
    subsample length and count are chosen to give target_detection_prob
    probability of overlapping a chimeric insert of size insert_fraction.

    Parameters
    ----------
    query_seq             : str   full query sequence string
    query_taxon           : str   e.g. "Solenodon_paradoxus"
    query_gene            : str   e.g. "12S_rRNA"
    query_accession       : str   e.g. "AY530070"
    query_orientation     : str   "+" or "-"
    insert_fraction       : float expected chimeric insert as fraction
                                  of total sequence length
    target_detection_prob : float desired probability of detecting
                                  a chimeric junction

    Returns
    -------
    dict with:
        sample_records   : list   SeqRecord objects ready for alignment
        sample_metadata  : list   dicts with start, end, length per sample
        n_samples        : int    number of subsamples drawn
        subsample_length : int    length of each subsample in bp
        seq_length       : int    length of query sequence
    """
    seq_length = len(query_seq)

    n = recommended_subsample_length(
        seq_length, insert_fraction=insert_fraction
    )
    n_samples = recommended_n_samples(
        seq_length            = seq_length,
        subsample_length      = n,
        target_detection_prob = target_detection_prob,
        insert_fraction       = insert_fraction,
    )

    print(f"\nQuery: {query_taxon} | {query_gene} | {query_accession}")
    print(f"Sequence: {seq_length}bp | Subsample: {n}bp | "
          f"Samples: {n_samples}")

    sample_records  = []
    sample_metadata = []

    for i in range(1, n_samples + 1):
        start  = random.randint(0, seq_length - n)
        sample = query_seq[start:start + n]

        header = (f"{query_taxon}|{query_gene}|orien:{query_orientation}|"
                  f"accession:{query_accession}|Sample{i}")

        sample_records.append(
            SeqRecord(Seq(sample), id=header, description="")
        )
        sample_metadata.append({
            "sample": i,
            "start":  start,
            "end":    start + n,
            "length": n,
        })
        print(f"  Sample {i}/{n_samples} | "
              f"position {start}-{start + n}")

    return {
        "sample_records":   sample_records,
        "sample_metadata":  sample_metadata,
        "n_samples":        n_samples,
        "subsample_length": n,
        "seq_length":       seq_length,
    }


# ─────────────────────────────────────────────────────────────────────────────
# INTERNAL TREE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _get_parent(tree, clade):
    """Return the parent node of clade in tree, or None if root."""
    for node in tree.find_clades(order="level"):
        if clade in node.clades:
            return node
    return None


def _in_same_clade(tree, target, reference_terminals) -> bool:
    """
    Return True if target terminal falls within the most recent
    common ancestor clade of reference_terminals.
    """
    try:
        mrca = tree.common_ancestor(reference_terminals)
        return target in mrca.get_terminals()
    except Exception:
        return False


# ─────────────────────────────────────────────────────────────────────────────
# PHYLOGENETIC PLACEMENT SCORING
# ─────────────────────────────────────────────────────────────────────────────

def score_phylogenetic_placement(treefile: str,
                                  id_map: dict,
                                  query_accession: str,
                                  query_taxon: str,
                                  taxonomy: dict,
                                  output_dir: str,
                                  save_path: str = None,
                                  bootstrap_threshold: float = 70.0) -> dict:
    """
    Score the phylogenetic placement of query subsamples relative to
    their expected taxonomic ingroup.

    For each subsample three metrics are computed:

    placement_ratio
        Distance to nearest outgroup divided by distance to nearest
        ingroup. Values > 1 indicate ingroup placement; values < 1
        indicate the subsample is closer to the outgroup than to its
        expected clade, which is a chimera signal.

    parent_bootstrap
        Bootstrap support of the node directly containing the query
        subsample. Low bootstrap (< bootstrap_threshold) means the
        placement is uncertain regardless of placement_ratio — this
        is important because low-support mixed clades (e.g. Mus inside
        Eulipotyphla at 30-40% bootstrap) are artefacts, not genuine
        signal.

    clade_ingroup_fraction
        Fraction of the query's immediate clade (parent node) that are
        ingroup taxa, excluding other query subsamples. A subsample in
        a clade that is 90% ingroup with 85% bootstrap is a confident
        ingroup placement. A subsample in a 50% ingroup clade with 40%
        bootstrap is not.

    bootstrap_weighted_ratio
        placement_ratio multiplied by parent_bootstrap / 100. This
        down-weights placements that happen to have a good ratio but
        land in poorly supported nodes — their contribution to the
        sequence-level summary is correspondingly reduced.

    confident_ingroup
        True only when all three conditions are met:
          - placement_ratio > 1 (closer to ingroup)
          - parent_bootstrap >= bootstrap_threshold
          - clade_ingroup_fraction >= 0.8

    is_outlier
        True when confident_ingroup is False. This replaces the simpler
        placement_ratio < 1 threshold from the previous version.

    Parameters
    ----------
    treefile            : str    path to IQ-TREE .treefile (Newick)
    id_map              : dict   short_id -> full_id mapping from
                                 add_queries_to_alignment()
    query_accession     : str    accession to identify query subsamples
    query_taxon         : str    e.g. "Solenodon_paradoxus"
    taxonomy            : dict   from load_taxonomy()
    output_dir          : str    directory to save figure
    save_path           : str    figure save directory (default:
                                 output_dir)
    bootstrap_threshold : float  minimum bootstrap to consider a
                                 placement confident (default 70)

    Returns
    -------
    dict with per-sample scores and sequence-level summary including:
        chimera_flag               : bool
        prop_confident_ingroup     : float
        mean_bootstrap_weighted_ratio : float
        mean_parent_bootstrap      : float
        sample_scores              : list of per-sample dicts
    """
    save_path = save_path or output_dir
    os.makedirs(save_path, exist_ok=True)

    # ── Determine expected ingroup from taxonomy ───────────────────────────────
    expected_order = get_expected_order(query_taxon, taxonomy)
    ingroup_genera = (get_ingroup_genera(expected_order, taxonomy)
                      if expected_order else set())

    print(f"\nExpected order : {expected_order}")
    print(f"Ingroup genera : {', '.join(sorted(ingroup_genera))}")

    # ── Load tree and restore full headers ────────────────────────────────────
    tree = Phylo.read(treefile, "newick")
    for clade in tree.get_terminals():
        if clade.name and clade.name in id_map:
            clade.name = id_map[clade.name]

    # ── Classify terminals as query / ingroup / outgroup ─────────────────────
    query_genus   = query_taxon.split("_")[0]
    all_terminals = tree.get_terminals()

    def is_query(name):
        return ("Sample" in (name or "")) and (
            query_accession in (name or "") or
            query_genus     in (name or "")
        )

    def get_genus(name):
        return (name or "").split("|")[0].split("_")[0]

    query_terminals    = [c for c in all_terminals if is_query(c.name)]
    ingroup_terminals  = [c for c in all_terminals
                          if not is_query(c.name)
                          and get_genus(c.name) in ingroup_genera]
    outgroup_terminals = [c for c in all_terminals
                          if not is_query(c.name)
                          and get_genus(c.name) not in ingroup_genera]

    print(f"\nTerminals identified:")
    print(f"  Query subsamples : {len(query_terminals)}")
    print(f"  Ingroup tips     : {len(ingroup_terminals)}")
    print(f"  Outgroup tips    : {len(outgroup_terminals)}")

    if not ingroup_terminals:
        raise ValueError(
            f"No ingroup terminals found for order '{expected_order}'. "
            f"Check that scaffold sequences include taxa from your "
            f"taxonomy CSV."
        )

    # ── Score each query subsample ────────────────────────────────────────────
    sample_scores = []

    for qt in query_terminals:

        def min_dist_to(terminals):
            dists = []
            for t in terminals:
                try:
                    dists.append((tree.distance(qt, t), t.name))
                except Exception:
                    continue
            return (min(dists, key=lambda x: x[0])
                    if dists else (None, None))

        ig_dist, ig_name = min_dist_to(ingroup_terminals)
        og_dist, og_name = min_dist_to(outgroup_terminals)

        if ig_dist is None:
            continue

        placement_ratio = (
            og_dist / ig_dist
            if og_dist is not None and ig_dist > 0
            else float("inf")
        )

        in_expected = _in_same_clade(tree, qt, ingroup_terminals)

        # Bootstrap of the node directly containing this subsample
        parent    = _get_parent(tree, qt)
        bootstrap = (float(parent.confidence)
                     if parent and parent.confidence is not None
                     else None)

        # Composition of parent clade — what fraction are ingroup?
        clade_ingroup_fraction = None
        if parent:
            clade_members = parent.get_terminals()
            non_query     = [c for c in clade_members
                             if not is_query(c.name)]
            if non_query:
                n_ig = sum(1 for c in non_query
                           if get_genus(c.name) in ingroup_genera)
                clade_ingroup_fraction = round(n_ig / len(non_query), 4)

        # Bootstrap-weighted placement ratio —
        # down-weights placements in poorly supported nodes
        if bootstrap is not None and placement_ratio != float("inf"):
            bootstrap_weighted_ratio = round(
                placement_ratio * (bootstrap / 100.0), 4
            )
        else:
            bootstrap_weighted_ratio = None

        # Confident ingroup requires all three conditions:
        # good ratio, sufficient bootstrap, predominantly ingroup clade
        confident_ingroup = (
            placement_ratio > 1.0
            and bootstrap is not None
            and bootstrap >= bootstrap_threshold
            and clade_ingroup_fraction is not None
            and clade_ingroup_fraction >= 0.8
        )

        score = {
            "sample_name":               qt.name,
            "dist_to_nearest_ingroup":   round(ig_dist, 6),
            "dist_to_nearest_outgroup":  (round(og_dist, 6)
                                          if og_dist is not None
                                          else None),
            "placement_ratio":           round(placement_ratio, 4),
            "nearest_ingroup_taxon":     ig_name,
            "nearest_outgroup_taxon":    og_name,
            "in_expected_clade":         in_expected,
            "parent_bootstrap":          bootstrap,
            "clade_ingroup_fraction":    clade_ingroup_fraction,
            "bootstrap_weighted_ratio":  bootstrap_weighted_ratio,
            "confident_ingroup":         confident_ingroup,
            "is_outlier":                not confident_ingroup,
        }
        sample_scores.append(score)

        flag = ("✓ CONFIDENT" if confident_ingroup
                else "~ UNCERTAIN" if placement_ratio > 1.0
                else "⚠ OUTLIER")
        print(f"\n  {flag}  {qt.name}")
        print(f"    nearest ingroup:          {ig_name} "
              f"(d={ig_dist:.4f})")
        if og_dist is not None:
            print(f"    nearest outgroup:         {og_name} "
                  f"(d={og_dist:.4f})")
        print(f"    placement ratio:          {placement_ratio:.4f}")
        print(f"    parent bootstrap:         {bootstrap}")
        print(f"    clade ingroup fraction:   {clade_ingroup_fraction}")
        print(f"    bootstrap-weighted ratio: {bootstrap_weighted_ratio}")
        print(f"    confident ingroup:        {confident_ingroup}")

    # ── Sequence-level summary ────────────────────────────────────────────────
    ratios     = [s["placement_ratio"] for s in sample_scores
                  if s["placement_ratio"] != float("inf")]
    bw_ratios  = [s["bootstrap_weighted_ratio"] for s in sample_scores
                  if s["bootstrap_weighted_ratio"] is not None]
    bs_values  = [s["parent_bootstrap"] for s in sample_scores
                  if s["parent_bootstrap"] is not None]
    confident  = [s for s in sample_scores if s["confident_ingroup"]]
    outliers   = [s for s in sample_scores if s["is_outlier"]]

    n_scored   = len(sample_scores)

    summary = {
        "query_accession":               query_accession,
        "query_taxon":                   query_taxon,
        "expected_order":                expected_order,
        "n_ingroup_tips":                len(ingroup_terminals),
        "n_outgroup_tips":               len(outgroup_terminals),
        "n_samples_scored":              n_scored,
        "n_outliers":                    len(outliers),
        "n_confident_ingroup":           len(confident),
        "prop_confident_ingroup":        (round(len(confident) /
                                               n_scored, 4)
                                          if n_scored else None),
        "prop_outlier_samples":          (round(len(outliers) /
                                               n_scored, 4)
                                          if n_scored else None),
        "mean_placement_ratio":          (round(float(np.mean(ratios)),
                                               4) if ratios else None),
        "min_placement_ratio":           (round(float(np.min(ratios)),
                                               4) if ratios else None),
        "max_placement_ratio":           (round(float(np.max(ratios)),
                                               4) if ratios else None),
        "placement_ratio_variance":      (round(float(np.var(ratios)),
                                               4) if ratios else None),
        "mean_bootstrap_weighted_ratio": (round(float(np.mean(bw_ratios)),
                                               4) if bw_ratios else None),
        "mean_parent_bootstrap":         (round(float(np.mean(bs_values)),
                                               4) if bs_values else None),
        # chimera_flag is True if any samples are outliers OR if no
        # samples achieved confident ingroup placement
        "chimera_flag":                  (len(outliers) > 0
                                          or len(confident) == 0),
        "sample_scores":                 sample_scores,
    }

    print(f"\n── Placement summary "
          f"────────────────────────────────────────────")
    print(f"  Expected order              : {expected_order}")
    print(f"  Confident ingroup placements: "
          f"{len(confident)}/{n_scored}")
    print(f"  Outlier samples             : "
          f"{len(outliers)}/{n_scored}")
    print(f"  Mean placement ratio        : "
          f"{summary['mean_placement_ratio']}")
    print(f"  Bootstrap-weighted ratio    : "
          f"{summary['mean_bootstrap_weighted_ratio']}")
    print(f"  Mean parent bootstrap       : "
          f"{summary['mean_parent_bootstrap']}")
    print(f"  Chimera flag                : {summary['chimera_flag']}")

    _plot_phylogeny(
        tree            = tree,
        query_taxon     = query_taxon,
        query_accession = query_accession,
        expected_order  = expected_order,
        outlier_names   = [s["sample_name"] for s in outliers],
        ingroup_genera  = ingroup_genera,
        output_dir      = save_path,
        base_name       = f"{query_taxon}_{query_accession}",
    )

    return summary


# ─────────────────────────────────────────────────────────────────────────────
# PLOT
# ─────────────────────────────────────────────────────────────────────────────

def _plot_phylogeny(tree, query_taxon, query_accession,
                    expected_order, outlier_names,
                    ingroup_genera, output_dir, base_name):
    """
    Plot the phylogeny with colour-coded labels:
        Green     = query subsample with confident ingroup placement
        Red       = query subsample that is an outlier
        Steel blue = scaffold ingroup taxon
        Dark grey  = scaffold outgroup taxon
    """
    query_genus = query_taxon.split("_")[0]

    def get_genus(name):
        return (name or "").split("|")[0].split("_")[0]

    def label_colour(name):
        is_q = "Sample" in (name or "") and (
            query_accession in (name or "") or
            query_genus     in (name or "")
        )
        if is_q:
            return "firebrick" if name in outlier_names else "seagreen"
        if get_genus(name) in ingroup_genera:
            return "steelblue"
        return "dimgrey"

    n_terminals = len(tree.get_terminals())
    fig, ax     = plt.subplots(
        figsize=(14, max(8, n_terminals * 0.35))
    )

    Phylo.draw(tree, axes=ax, do_show=False, label_colors=label_colour)

    ax.set_title(
        f"{query_taxon.replace('_', ' ')} ({query_accession}) — "
        f"phylogenetic placement\n"
        f"Expected order: {expected_order}   |   "
        f"Green = confident ingroup   "
        f"Red = outlier   "
        f"Blue = scaffold ingroup   "
        f"Grey = scaffold outgroup",
        fontsize=10,
    )
    ax.set_xlabel("Branch length (substitutions per site)")
    plt.tight_layout()

    fig_path = os.path.join(output_dir, f"{base_name}_tree.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Figure saved: {fig_path}")
    return fig_path


# ─────────────────────────────────────────────────────────────────────────────
# SCAFFOLD TREE VISUALISATION (separate from query verification plots)
# ─────────────────────────────────────────────────────────────────────────────

def plot_scaffold_tree(treefile: str,
                        id_map_path: str = None,
                        save_path: str = None):
    """
    Plot a scaffold tree with readable labels. Restores full taxon
    names from id_map if provided.

    Parameters
    ----------
    treefile    : str   path to IQ-TREE .treefile
    id_map_path : str   path to _id_map.json. If provided, short IDs
                        are replaced with Genus_species_accession labels.
    save_path   : str   path to save PNG (default: treefile + .png)
    """
    tree = Phylo.read(treefile, "newick")

    if id_map_path and Path(id_map_path).exists():
        with open(id_map_path, "r") as f:
            id_map = json.load(f)
        for clade in tree.get_terminals():
            if clade.name and clade.name in id_map:
                full_id = id_map[clade.name]
                parts   = full_id.split("|")
                taxon   = parts[0]
                acc     = next(
                    (p.split(":", 1)[1] for p in parts
                     if p.startswith("accession:")), ""
                )
                clade.name = f"{taxon}_{acc}" if acc else taxon

    n_terminals = len(tree.get_terminals())
    fig, ax     = plt.subplots(
        figsize=(14, max(8, n_terminals * 0.4))
    )

    Phylo.draw(tree, axes=ax, do_show=False)

    gene_name = Path(treefile).stem
    ax.set_title(f"Scaffold tree — {gene_name}", fontsize=12)
    ax.set_xlabel("Branch length (substitutions per site)")
    plt.tight_layout()

    if save_path is None:
        save_path = str(treefile).replace(".treefile", "_tree.png")

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figure saved: {save_path}")
    return save_path


# ─────────────────────────────────────────────────────────────────────────────
# SINGLE SEQUENCE PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def verify_by_phylogeny(query_seq: str,
                         query_taxon: str,
                         query_gene: str,
                         query_accession: str,
                         query_orientation: str,
                         scaffold_fasta_aligned: str,
                         output_dir: str,
                         taxonomy: dict,
                         iqtree_bin: str = "iqtree2",
                         model: str = "GTR+G",
                         bootstrap: int = 1000,
                         insert_fraction: float = 0.30,
                         target_detection_prob: float = 0.90,
                         n_threads: int = 4,
                         seed: int = 12345,
                         email: str = "om380@cam.ac.uk",
                         min_coverage: float = 0.5,
                         bootstrap_threshold: float = 70.0,
                         save_path: str = None) -> dict:
    """
    Verify a single sequence by adding subsamples to an existing
    scaffold alignment and scoring their phylogenetic placement.

    Does NOT re-align the scaffold — only adds query subsamples using
    profile alignment, preserving scaffold alignment quality.

    Output saved to: {output_dir}/{query_taxon}_{query_accession}_phylo/

    Parameters
    ----------
    query_seq              : str    full query sequence
    query_taxon            : str    e.g. "Solenodon_paradoxus"
    query_gene             : str    e.g. "12S_rRNA"
    query_accession        : str    e.g. "AY530070"
    query_orientation      : str    "+" or "-"
    scaffold_fasta_aligned : str    path to full-header aligned scaffold
                                    FASTA (_aln.fasta, not _short_aln)
    output_dir             : str    parent directory for outputs
    taxonomy               : dict   from load_taxonomy()
    iqtree_bin             : str    IQ-TREE binary path
    model                  : str    substitution model
    bootstrap              : int    ultrafast bootstrap replicates
    insert_fraction        : float  expected chimeric insert fraction
    target_detection_prob  : float  desired P(detect chimera)
    n_threads              : int    CPU threads
    seed                   : int    random seed
    email                  : str    EBI ClustalOmega email
    min_coverage           : float  end-trimming coverage threshold
    bootstrap_threshold    : float  minimum bootstrap for confident
                                    ingroup placement (default 70)
    save_path              : str    figure directory (default: query dir)

    Returns
    -------
    dict combining outputs from all four pipeline steps
    """
    query_dir = os.path.join(
        output_dir,
        f"{query_taxon}_{query_accession}_phylo"
    )
    os.makedirs(query_dir, exist_ok=True)

    # 1. Generate query subsamples
    subsample_result = build_query_subsamples(
        query_seq             = query_seq,
        query_taxon           = query_taxon,
        query_gene            = query_gene,
        query_accession       = query_accession,
        query_orientation     = query_orientation,
        insert_fraction       = insert_fraction,
        target_detection_prob = target_detection_prob,
    )

    # 2. Profile alignment: add subsamples to scaffold
    align_result = add_queries_to_alignment(
        scaffold_fasta_aligned = scaffold_fasta_aligned,
        query_records          = subsample_result["sample_records"],
        output_dir             = query_dir,
        query_taxon            = query_taxon,
        query_accession        = query_accession,
        email                  = email,
        min_coverage           = min_coverage,
    )

    # 3. IQ-TREE on trimmed short-ID Nexus
    iqtree_result = run_iqtree(
        nexus_alignment = align_result["trimmed_nexus_short"],
        output_dir      = query_dir,
        iqtree_bin      = iqtree_bin,
        model           = model,
        bootstrap       = bootstrap,
        n_threads       = n_threads,
        seed            = seed,
    )

    # 4. Score placement with bootstrap-weighted metrics
    placement_result = score_phylogenetic_placement(
        treefile            = iqtree_result["treefile"],
        id_map              = align_result["id_map"],
        query_accession     = query_accession,
        query_taxon         = query_taxon,
        taxonomy            = taxonomy,
        output_dir          = query_dir,
        save_path           = save_path or query_dir,
        bootstrap_threshold = bootstrap_threshold,
    )

    return {
        **subsample_result,
        **align_result,
        **iqtree_result,
        **placement_result,
        "query_dir": query_dir,
    }


# ─────────────────────────────────────────────────────────────────────────────
# BATCH PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def verify_batch_by_phylogeny(query_fasta: str,
                               scaffold_fasta_aligned: str,
                               output_dir: str,
                               taxonomy: dict,
                               iqtree_bin: str = "iqtree2",
                               model: str = "GTR+G",
                               bootstrap: int = 1000,
                               insert_fraction: float = 0.30,
                               target_detection_prob: float = 0.90,
                               n_threads: int = 4,
                               seed: int = 12345,
                               email: str = "om380@cam.ac.uk",
                               min_coverage: float = 0.5,
                               bootstrap_threshold: float = 70.0,
                               interval: int = 30) -> list:
    """
    Run verify_by_phylogeny() on every sequence in a FASTA file.

    FASTA header format expected:
        >Species_name|gene|orien:+/-|accession:XXXX

    Parameters
    ----------
    query_fasta            : str   FASTA of sequences to verify
    scaffold_fasta_aligned : str   full-header aligned scaffold FASTA
    output_dir             : str   parent directory for all outputs
    taxonomy               : dict  from load_taxonomy()
    iqtree_bin             : str   IQ-TREE binary
    model                  : str   substitution model
    bootstrap              : int   ultrafast bootstrap replicates
    insert_fraction        : float expected chimeric insert fraction
    target_detection_prob  : float desired P(detect chimera)
    n_threads              : int   CPU threads
    seed                   : int   random seed
    email                  : str   EBI ClustalOmega email
    min_coverage           : float end-trimming coverage threshold
    bootstrap_threshold    : float minimum bootstrap for confident
                                   ingroup placement
    interval               : int   seconds to wait between submissions

    Returns
    -------
    list of result dicts, one per sequence
    """
    records = list(SeqIO.parse(query_fasta, "fasta"))
    total   = len(records)
    results = []

    print(f"\nBatch phylogenetic verification: {total} sequences")
    print(f"Scaffold alignment : {scaffold_fasta_aligned}")
    print(f"Output directory   : {output_dir}\n")

    for i, record in enumerate(records, 1):
        parts       = record.description.split("|")
        taxon       = parts[0].strip()
        gene        = parts[1].strip() if len(parts) > 1 else "unknown"
        orientation = "+"
        accession   = None

        for p in parts:
            if p.startswith("orien:"):
                orientation = p.split(":", 1)[1].strip()
            if p.startswith("accession:"):
                accession = p.split(":", 1)[1].strip()

        if not accession:
            print(f"  [{i}/{total}] Skipping {taxon} — "
                  f"no accession in header")
            continue

        print(f"\n[{i}/{total}] {taxon} | {gene} | {accession}")

        try:
            result = verify_by_phylogeny(
                query_seq              = str(record.seq),
                query_taxon            = taxon,
                query_gene             = gene,
                query_accession        = accession,
                query_orientation      = orientation,
                scaffold_fasta_aligned = scaffold_fasta_aligned,
                output_dir             = output_dir,
                taxonomy               = taxonomy,
                iqtree_bin             = iqtree_bin,
                model                  = model,
                bootstrap              = bootstrap,
                insert_fraction        = insert_fraction,
                target_detection_prob  = target_detection_prob,
                n_threads              = n_threads,
                seed                   = seed,
                email                  = email,
                min_coverage           = min_coverage,
                bootstrap_threshold    = bootstrap_threshold,
            )
            result["_header"] = record.description
            results.append(result)

            print(f"  ✓ Done | chimera_flag="
                  f"{result['chimera_flag']} | "
                  f"confident_ingroup="
                  f"{result['n_confident_ingroup']}/"
                  f"{result['n_samples_scored']} | "
                  f"mean_bs_weighted_ratio="
                  f"{result['mean_bootstrap_weighted_ratio']}")

        except Exception as e:
            print(f"  ✗ Failed: {e}")
            import traceback
            traceback.print_exc()

        if i < total:
            print(f"  Waiting {interval}s before next submission...")
            time.sleep(interval)

    print(f"\nBatch complete: {len(results)}/{total} succeeded")
    return results