# _build_reference_fastas.py

import os
import glob
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def build_combined_fasta(
        fasta_paths: list,
        output_path: str,
) -> dict:
    """
    Combine multiple FASTA files into a single file, deduplicating
    by (taxon, gene) keeping the longest sequence per pair.

    All input FASTAs must use pipe-delimited headers:
        >Taxon_name|GENE|orien:+|accession:XXXX

    Parameters
    ----------
    fasta_paths : list   list of FASTA file paths to combine
    output_path : str    path to write combined output FASTA

    Returns
    -------
    dict with n_input, n_written, n_deduplicated, output_path
    """
    found   = {}
    n_input = 0
    n_dedup = 0

    for fasta_path in fasta_paths:
        if not fasta_path or not Path(fasta_path).exists():
            print(f"  Warning: not found: {fasta_path}")
            continue

        print(f"  Scanning: {Path(fasta_path).name}")
        n_file = 0

        for record in SeqIO.parse(fasta_path, "fasta"):
            n_input += 1
            taxon, gene = _parse_pipe_header(record.id)
            if taxon is None or gene is None:
                continue

            key      = (taxon, gene)
            existing = found.get(key)
            if existing is None or len(record.seq) > len(existing.seq):
                if existing is not None:
                    n_dedup += 1
                found[key] = record
                n_file    += 1

        print(f"    {n_file} sequences loaded")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(list(found.values()), output_path, "fasta")

    print(f"  Combined FASTA written: {output_path} "
          f"({len(found)} sequences, {n_dedup} deduplicated)")

    return {
        "n_input":        n_input,
        "n_written":      len(found),
        "n_deduplicated": n_dedup,
        "output_path":    output_path,
    }


def _parse_pipe_header(record_id: str) -> tuple:
    """
    Parse (taxon, gene) from a pipe-delimited FASTA header.
    Handles:
        Taxon_name|GENE|orien:+|accession:XXXX
        Taxon_name|GENE
    Returns (None, None) if not parseable.
    """
    if '|' not in record_id:
        return None, None

    parts = record_id.split('|')
    if len(parts) < 2:
        return None, None

    taxon = parts[0].strip()
    gene  = parts[1].strip().upper()

    if (not taxon
            or '_' not in taxon
            or not taxon[0].isupper()):
        return None, None

    if not gene:
        return None, None

    rrna_map = {
        "12S_RRNA": "12S_rRNA",
        "16S_RRNA": "16S_rRNA",
        "12S":      "12S_rRNA",
        "16S":      "16S_rRNA",
    }
    gene = rrna_map.get(gene, gene)

    return taxon, gene


def build_reference_fastas(
        cleared_by_taxon: dict,
        combined_fasta: str,
        output_dir: str,
        calibration_fasta: str = None,
        min_taxa_per_gene: int = 3,
        genes_to_include: list = None,
) -> dict:
    """
    Extract full-length sequences for taxa that passed BLAST
    verification (cleared sequences only — not flagged, not review)
    from a combined FASTA file, and write per-gene FASTA files ready
    for alignment with ClustalOmega.

    The cleared_by_taxon dict should contain ONLY sequences that were
    not flagged by the BLAST model. This ensures the reference
    alignment for family distance verification is built from sequences
    believed to be genuine, giving a reliable baseline for within-
    family distance calibration.

    Flagged sequences (LIKELY_TP, REVIEW, LIKELY_FP, ORPHAN_FP) are
    NOT included in the reference alignment — they are tested against
    it.

    Parameters
    ----------
    cleared_by_taxon  : dict   taxon -> list of cleared genes.
                               Use lists["cleared_by_taxon"] from
                               extract_verification_lists(), which
                               contains only sequences that passed
                               BLAST verification.
    combined_fasta    : str    path to combined FASTA containing all
                               sequences (from build_combined_fasta())
    output_dir        : str    directory to write [gene].fasta files
    calibration_fasta : str    optional separate calibration FASTA
                               (calibration sequences are also genuine
                               and should be included in reference)
    min_taxa_per_gene : int    skip genes with fewer taxa
    genes_to_include  : list   optional gene filter

    Returns
    -------
    dict with written_genes, skipped_genes, not_found, n_files_written
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    targets = set()
    for taxon, genes in cleared_by_taxon.items():
        for gene in genes:
            gene_norm = _normalise_gene(gene)
            if genes_to_include is None or gene_norm in [
                    _normalise_gene(g) for g in genes_to_include]:
                targets.add((taxon, gene_norm))

    print(f"\n{'='*60}")
    print(f"BUILD REFERENCE FASTAS")
    print(f"{'='*60}")
    print(f"  Combined FASTA : {combined_fasta}")
    print(f"  Output dir     : {output_dir}")
    print(f"  Target pairs   : {len(targets)}")
    print(f"  (Using cleared sequences only — flagged sequences "
          f"excluded from reference alignment)")

    if len(targets) == 0:
        print(f"\n  WARNING: No targets — cleared_by_taxon is empty.")
        return {
            "written_genes":   {},
            "skipped_genes":   {},
            "not_found":       [],
            "n_files_written": 0,
        }

    found      = {}
    scan_paths = [combined_fasta]
    if calibration_fasta:
        scan_paths.append(calibration_fasta)

    for fasta_path in scan_paths:
        if not fasta_path or not Path(fasta_path).exists():
            print(f"  Warning: FASTA not found: {fasta_path}")
            continue

        print(f"\n  Scanning: {Path(fasta_path).name}")
        n_scanned = 0
        n_matched = 0

        for record in SeqIO.parse(fasta_path, "fasta"):
            n_scanned += 1
            taxon, gene = _parse_pipe_header(record.id)

            if taxon is None or gene is None:
                continue

            key = (taxon, gene)
            if key not in targets:
                continue

            existing = found.get(key)
            if existing is None or len(record.seq) > len(existing.seq):
                found[key]  = record
                n_matched  += 1

        print(f"    Scanned {n_scanned}, matched {n_matched}")

    by_gene       = defaultdict(list)
    for (taxon, gene), record in found.items():
        by_gene[gene].append((taxon, record))

    not_found     = sorted(targets - set(found.keys()))
    written_genes = {}
    skipped_genes = {}
    n_files       = 0

    if not_found:
        by_missing_gene = defaultdict(list)
        for taxon, gene in not_found:
            by_missing_gene[gene].append(taxon)
        print(f"\n  Not found ({len(not_found)} pairs):")
        for gene, taxa in sorted(by_missing_gene.items()):
            print(f"    {gene:<15} {', '.join(sorted(taxa)[:5])}"
                  f"{'...' if len(taxa) > 5 else ''}")

    print(f"\n  Writing per-gene FASTA files...")
    for gene, taxon_records in sorted(by_gene.items()):
        if len(taxon_records) < min_taxa_per_gene:
            skipped_genes[gene] = (
                f"Only {len(taxon_records)} taxa "
                f"(minimum {min_taxa_per_gene})"
            )
            continue

        out_path       = Path(output_dir) / f"{gene}.fasta"
        output_records = []
        taxa_written   = []

        for taxon, record in sorted(taxon_records, key=lambda x: x[0]):
            output_records.append(SeqRecord(
                seq         = record.seq,
                id          = taxon,
                name        = "",
                description = "",
            ))
            taxa_written.append(taxon)

        SeqIO.write(output_records, str(out_path), "fasta")
        written_genes[gene] = taxa_written
        n_files += 1
        print(f"    {gene:<15} {len(taxa_written):>3} taxa "
              f"→ {out_path.name}")

    print(f"\n{'='*60}")
    print(f"BUILD REFERENCE FASTAS COMPLETE")
    print(f"{'='*60}")
    print(f"  Files written : {n_files}")
    print(f"  Genes skipped : {len(skipped_genes)}")
    print(f"  Pairs missing : {len(not_found)}")

    return {
        "written_genes":   written_genes,
        "skipped_genes":   skipped_genes,
        "not_found":       not_found,
        "n_files_written": n_files,
        "output_dir":      output_dir,
    }


def _normalise_gene(gene: str) -> str:
    rrna_map = {
        "12S_RRNA": "12S_rRNA",
        "16S_RRNA": "16S_rRNA",
        "12S":      "12S_rRNA",
        "16S":      "16S_rRNA",
    }
    upper = gene.strip().upper()
    return rrna_map.get(upper, upper)


def add_calibration_to_qualified(
        qualified: dict,
        calibration_set: set,
        gbk_index: dict,
        calibration_fasta: str,
) -> dict:
    """
    Add calibration set sequences to the qualified dict so they are
    included in the reference FASTA build.

    Calibration sequences were never scored by score_from_csv() so
    they never appear in qualify_cleared_sequences() output. This
    function finds which taxa the calibration accessions belong to
    (via the GBK index) and adds all their genes to qualified.

    Parameters
    ----------
    qualified         : dict   existing qualified dict
    calibration_set   : set    accession strings from calibration
    gbk_index         : dict   from build_gbk_index()
    calibration_fasta : str    path to calibration FASTA

    Returns
    -------
    dict   updated qualified dict with calibration sequences added
    """
    import copy

    updated = copy.deepcopy(qualified)

    acc_to_taxon = {}
    for taxon, accessions in gbk_index.items():
        for acc in accessions:
            acc_to_taxon[acc] = taxon

    cal_taxa_genes = defaultdict(set)

    if calibration_fasta and Path(calibration_fasta).exists():
        for record in SeqIO.parse(calibration_fasta, "fasta"):
            taxon, gene = _parse_pipe_header(record.id)
            if taxon is None or gene is None:
                continue
            taxon_accessions = gbk_index.get(taxon, [])
            if any(acc in calibration_set
                   for acc in taxon_accessions):
                cal_taxa_genes[taxon].add(gene)

    n_added = 0
    for taxon, genes in cal_taxa_genes.items():
        if taxon not in updated:
            updated[taxon] = list(genes)
            n_added += len(genes)
        else:
            existing       = set(updated[taxon])
            updated[taxon] = list(existing | genes)
            n_added       += len(genes - existing)

    print(f"  Added {n_added} calibration (taxon, gene) pairs "
          f"from {len(cal_taxa_genes)} taxa to qualified set")

    return updated


def resolve_accession_keys(cleared_by_taxon: dict,
                            results_csv: str) -> dict:
    """
    Convert a cleared_by_taxon dict keyed by accession strings
    to one keyed by actual Genus_species taxon names, using the
    taxon column from blast_results.csv for the mapping.

    Parameters
    ----------
    cleared_by_taxon : dict   accession or taxon -> list of genes
    results_csv      : str    path to blast_results.csv

    Returns
    -------
    dict   taxon_name -> list of genes, with accession keys replaced
           by Genus_species names where resolvable. Entries that
           already look like taxon names are kept as-is.
    """
    import csv

    # Build accession -> taxon map from CSV
    acc_to_taxon = {}
    with open(results_csv, "r", newline="",
              encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc   = row.get("accession", "").strip()
            taxon = row.get("taxon", "").strip()
            if acc and taxon:
                acc_to_taxon[acc] = taxon

    resolved    = {}
    n_resolved  = 0
    n_kept      = 0
    n_failed    = 0

    for key, genes in cleared_by_taxon.items():
        # Check if key already looks like a taxon name
        if ("_" in key
                and key[0].isupper()
                and not any(c.isdigit() for c in key)):
            resolved[key] = genes
            n_kept += 1
            continue

        # Try to resolve accession to taxon name
        taxon = acc_to_taxon.get(key)

        # Try without version suffix (LC124901.1 -> LC124901)
        if taxon is None and "." in key:
            taxon = acc_to_taxon.get(key.split(".")[0])

        if taxon and "_" in taxon and taxon[0].isupper():
            if taxon not in resolved:
                resolved[taxon] = list(genes)
            else:
                # Merge genes if taxon already present
                existing = set(resolved[taxon])
                resolved[taxon] = list(existing | set(genes))
            n_resolved += 1
        else:
            # Cannot resolve — keep original key as fallback
            resolved[key] = genes
            n_failed += 1

    print(f"  Resolved accession keys to taxon names:")
    print(f"    Already taxon names : {n_kept}")
    print(f"    Resolved from CSV   : {n_resolved}")
    print(f"    Could not resolve   : {n_failed}")
    print(f"    Final taxon count   : {len(resolved)}")

    return resolved