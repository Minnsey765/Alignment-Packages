from collections import defaultdict
import csv
import json
import math
import numpy as np
from pathlib import Path
from collections import defaultdict
import copy
from bioinf_packages.alignment_funcs._read_gbk import read_gbk

def extract_verification_lists(
        result: dict,
        output_dir: str = None,
) -> dict:
    """
    Extract flagged and cleared sequence lists from score_from_csv()
    output. Groups by taxon (accession field) and gene.

    Parameters
    ----------
    result     : dict   return value of score_from_csv()
    output_dir : str    optional directory to write text files

    Returns
    -------
    dict with:
        flagged_by_taxon  : dict   taxon -> list of flagged genes
        cleared_by_taxon  : dict   taxon -> list of cleared genes
        flagged_list      : list   (taxon, gene) tuples
        cleared_list      : list   (taxon, gene) tuples
        flagged_taxa      : set    unique taxa with any chimeric gene
        cleared_taxa      : set    taxa with ALL genes cleared
        mixed_taxa        : set    taxa with both flagged and cleared genes
    """
    flagged_by_taxon  = defaultdict(list)
    cleared_by_taxon  = defaultdict(list)

    for (acc, gene), call in result["chimera_calls"].items():
        if call["chimera"]:
            flagged_by_taxon[acc].append(gene)
        else:
            cleared_by_taxon[acc].append(gene)

    flagged_taxa = set(flagged_by_taxon.keys())
    cleared_taxa = {t for t in cleared_by_taxon
                    if t not in flagged_taxa}
    mixed_taxa   = {t for t in flagged_taxa
                    if t in cleared_by_taxon}

    print(f"\n{'='*60}")
    print(f"VERIFICATION SUMMARY")
    print(f"{'='*60}")
    print(f"  Total sequences scored  : {result['n_query_seqs']}")
    print(f"  Flagged chimeric        : {result['n_chimeras']}")
    print(f"  Cleared genuine         : {result['n_genuine']}")
    print(f"\n  Unique taxa with any flag : {len(flagged_taxa)}")
    print(f"  Unique taxa fully cleared : {len(cleared_taxa)}")
    print(f"  Taxa with mixed results   : {len(mixed_taxa)}")

    if mixed_taxa:
        print(f"\n  Mixed taxa (some genes flagged, some cleared):")
        for taxon in sorted(mixed_taxa):
            flagged = sorted(flagged_by_taxon[taxon])
            cleared = sorted(cleared_by_taxon[taxon])
            print(f"    {taxon}")
            print(f"      Flagged : {', '.join(flagged)}")
            print(f"      Cleared : {', '.join(cleared)}")

    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Flagged sequences — one line per (taxon, gene)
        flagged_path = Path(output_dir) / "flagged_sequences.txt"
        with open(flagged_path, "w") as f:
            f.write("taxon\tgene\ttrigger\tprop_neg\n")
            for (acc, gene) in sorted(result["flagged"]):
                c = result["chimera_calls"][(acc, gene)]
                f.write(f"{acc}\t{gene}\t"
                        f"{c['trigger']}\t{c['prop_neg']:.4f}\n")
        print(f"\n  Flagged list saved: {flagged_path}")

        # Cleared sequences
        cleared_path = Path(output_dir) / "cleared_sequences.txt"
        with open(cleared_path, "w") as f:
            f.write("taxon\tgene\tprop_neg\tmean_posterior\n")
            for (acc, gene), c in sorted(
                    result["chimera_calls"].items()):
                if not c["chimera"]:
                    f.write(f"{acc}\t{gene}\t"
                            f"{c['prop_neg']:.4f}\t"
                            f"{c['mean_posterior']:.4f}\n")
        print(f"  Cleared list saved: {cleared_path}")

        # Taxa summary
        taxa_path = Path(output_dir) / "taxa_summary.txt"
        with open(taxa_path, "w") as f:
            f.write("taxon\tstatus\tflagged_genes\tcleared_genes\n")
            all_taxa = flagged_taxa | cleared_taxa | mixed_taxa
            for taxon in sorted(all_taxa):
                if taxon in mixed_taxa:
                    status = "MIXED"
                elif taxon in flagged_taxa:
                    status = "ALL_FLAGGED"
                else:
                    status = "CLEARED"
                flagged = ','.join(sorted(
                    flagged_by_taxon.get(taxon, [])))
                cleared = ','.join(sorted(
                    cleared_by_taxon.get(taxon, [])))
                f.write(f"{taxon}\t{status}\t"
                        f"{flagged}\t{cleared}\n")
        print(f"  Taxa summary saved: {taxa_path}")

    return {
        "flagged_by_taxon": dict(flagged_by_taxon),
        "cleared_by_taxon": dict(cleared_by_taxon),
        "flagged_list":     list(result["flagged"]),
        "cleared_list":     [(acc, gene)
                             for (acc, gene), c
                             in result["chimera_calls"].items()
                             if not c["chimera"]],
        "flagged_taxa":     flagged_taxa,
        "cleared_taxa":     cleared_taxa,
        "mixed_taxa":       mixed_taxa,
    }

def build_gbk_index(gbk_folder: str, glossary: str) -> dict:
    from pathlib import Path
    from collections import defaultdict

    index     = defaultdict(list)
    gbk_files = list(Path(gbk_folder).glob("*.gbk"))
    gbk_files += list(Path(gbk_folder).glob("*.gb"))

    print(f"Indexing {len(gbk_files)} GBK files...")
    n_failed = 0
    for gbk_path in gbk_files:
        accession = gbk_path.stem
        try:
            genes = read_gbk(accession, gbk_folder, glossary)
            if genes:
                taxon = genes[0]["species"]
                index[taxon].append(accession)
        except Exception as e:
            n_failed += 1
            if n_failed <= 5:   # print first 5 failures only
                print(f"  Failed {accession}: {e}")

    print(f"  {len(index)} unique taxa indexed "
          f"({n_failed} files failed)")
    return dict(index)


def qualify_cleared_sequences(
        cleared_by_taxon: dict,
        gbk_folder: str,
        glossary: str,
        gbk_index: dict = None,
        calibration_accessions: set = None,
) -> dict:
    """
    Filter cleared sequences to those with at least one PubMed
    reference in their GenBank record, plus calibration sequences.

    Parameters
    ----------
    cleared_by_taxon       : dict   from extract_verification_lists()
    gbk_folder             : str    directory containing .gbk files
    glossary               : str    path to glossary CSV for read_gbk()
    gbk_index              : dict   optional pre-built taxon->accessions
                                    index from build_gbk_index(). If
                                    None, built automatically (slow).
    calibration_accessions : set    accessions from original genuine
                                    calibration set — included
                                    unconditionally.

    Returns
    -------
    dict with qualified, unqualified, pubmed_found, not_found,
              n_qualified, n_unqualified, n_calibration
    """
    from Bio import SeqIO
    from pathlib import Path

    calibration_accessions = calibration_accessions or set()

    if gbk_index is None:
        gbk_index = build_gbk_index(gbk_folder, glossary)

    qualified    = defaultdict(list)
    unqualified  = defaultdict(list)
    pubmed_found = defaultdict(list)
    not_found    = []
    n_calibration = 0

    print(f"\n{'='*60}")
    print(f"QUALIFY CLEARED SEQUENCES")
    print(f"{'='*60}")
    print(f"  Taxa to check     : {len(cleared_by_taxon)}")
    print(f"  Calibration set   : {len(calibration_accessions)}")

    for taxon, genes in sorted(cleared_by_taxon.items()):

        accessions = gbk_index.get(taxon, [])

        if not accessions:
            not_found.append(taxon)
            for gene in genes:
                if taxon in calibration_accessions:
                    qualified[taxon].append(gene)
                    n_calibration += 1
                else:
                    unqualified[taxon].append(gene)
            continue

        # Check PubMed references across all GBK files for this taxon
        pubmed_ids = []
        for acc in accessions:
            gbk_path = Path(gbk_folder) / f"{acc}.gbk"
            if not gbk_path.exists():
                gbk_path = Path(gbk_folder) / f"{acc}.gb"
            if not gbk_path.exists():
                continue
            try:
                for record in SeqIO.parse(str(gbk_path), "genbank"):
                    for ref in record.annotations.get("references", []):
                        pid = getattr(ref, "pubmed_id", "").strip()
                        if pid:
                            pubmed_ids.append(pid)
            except Exception:
                continue

        pubmed_ids = list(set(pubmed_ids))
        has_pubmed = len(pubmed_ids) > 0

        if has_pubmed:
            pubmed_found[taxon] = pubmed_ids

        for gene in genes:
            in_calibration = any(
                acc in calibration_accessions for acc in accessions
            )
            if has_pubmed or in_calibration:
                qualified[taxon].append(gene)
                if in_calibration and not has_pubmed:
                    n_calibration += 1
            else:
                unqualified[taxon].append(gene)

    n_qualified   = sum(len(g) for g in qualified.values())
    n_unqualified = sum(len(g) for g in unqualified.values())

    print(f"\n  Qualified (taxon, gene) pairs  : {n_qualified}")
    print(f"    via PubMed                   : "
          f"{n_qualified - n_calibration}")
    print(f"    via calibration set          : {n_calibration}")
    print(f"  Unqualified pairs              : {n_unqualified}")
    print(f"  Taxa with no GBK found         : {len(not_found)}")

    if not_found:
        print(f"\n  Taxa without GBK files:")
        for t in sorted(not_found):
            print(f"    {t}")

    if unqualified:
        print(f"\n  Unqualified taxa:")
        for taxon, genes in sorted(unqualified.items()):
            print(f"    {taxon:<40} {', '.join(sorted(genes))}")

    return {
        "qualified":     dict(qualified),
        "unqualified":   dict(unqualified),
        "pubmed_found":  dict(pubmed_found),
        "not_found":     not_found,
        "n_qualified":   n_qualified,
        "n_unqualified": n_unqualified,
        "n_calibration": n_calibration,
    }

def assess_flagged_sequences(
        result: dict,
        qual: dict,
        gbk_folder: str,
        glossary: str,
        gbk_index: dict = None,
) -> dict:
    """
    Assess flagged sequences using the qualified cleared set as
    context. Produces a recommended action for each flagged sequence.

    Recommendations
    ---------------
    LIKELY_FP    — flag is isolated (single gene flagged for this
                   taxon, all others cleared and qualified), trigger
                   is post_var (noisy BLAST run), and prop_neg is
                   below 0.5. Likely a false positive.
    REVIEW       — ambiguous — multiple genes flagged or trigger is
                   prop_neg with high prop_neg value.
    LIKELY_TP    — most genes flagged for this taxon, or prop_neg
                   very high (>0.7), or sequence has no qualified
                   cleared genes to provide context.
    ORPHAN_FP    — taxon is in known orphan taxa list and all other
                   genes for this taxon are rescued/cleared — orphan
                   taxa routinely produce false positives.

    Parameters
    ----------
    result     : dict   return value of score_from_csv()
    qual       : dict   return value of qualify_cleared_sequences()
    gbk_folder : str    directory containing .gbk files
    glossary   : str    path to glossary CSV
    gbk_index  : dict   pre-built index from build_gbk_index()

    Returns
    -------
    dict with:
        assessments    : dict  (taxon, gene) -> assessment dict
        likely_fp      : list  (taxon, gene) pairs likely false positive
        review         : list  (taxon, gene) pairs needing review
        likely_tp      : list  (taxon, gene) pairs likely true positive
        orphan_fp      : list  (taxon, gene) pairs orphan false positives
        exclude_list   : list  (taxon, gene) pairs to exclude from
                               final dataset (LIKELY_TP only)
        retain_list    : list  (taxon, gene) pairs to retain despite
                               flag (LIKELY_FP + ORPHAN_FP)
    """
    from collections import defaultdict

    orphan_taxa = set(
        result.get("orphan_taxa", [])
        # fall back to empty if not stored in result
    )

    # Load orphan taxa from calibration JSON if available
    # (score_from_csv doesn't currently return these — add if needed)

    flagged_by_taxon  = defaultdict(list)
    cleared_by_taxon  = defaultdict(list)
    qualified_by_taxon = qual.get("qualified", {})

    for (taxon, gene), call in result["chimera_calls"].items():
        if call["chimera"]:
            flagged_by_taxon[taxon].append((gene, call))
        else:
            cleared_by_taxon[taxon].append(gene)

    assessments = {}
    likely_fp   = []
    review      = []
    likely_tp   = []
    orphan_fp   = []

    for taxon, flagged_genes in sorted(flagged_by_taxon.items()):

        n_flagged   = len(flagged_genes)
        n_cleared   = len(cleared_by_taxon.get(taxon, []))
        n_qualified = len(qualified_by_taxon.get(taxon, []))
        is_orphan   = any(taxon.startswith(g)
                          for g in orphan_taxa)

        for gene, call in flagged_genes:
            trigger  = call.get("trigger", "")
            prop_neg = call.get("prop_neg", 0.0)
            post_var = call.get("posterior_variance", 0.0)

            # Orphan taxon — orphan rescue fires in calibration but
            # not always in verification. Check if other genes cleared.
            if is_orphan and n_cleared > 0:
                recommendation = "ORPHAN_FP"
                orphan_fp.append((taxon, gene))

            # Single gene flagged, most genes cleared and qualified,
            # trigger is variance-based with moderate prop_neg
            elif (n_flagged == 1
                  and n_qualified >= 3
                  and trigger == "posterior_variance"
                  and prop_neg < 0.5):
                recommendation = "LIKELY_FP"
                likely_fp.append((taxon, gene))

            # Most genes flagged — strong chimera signal
            elif (n_flagged > n_cleared
                  or prop_neg > 0.7
                  or (n_qualified == 0 and n_flagged > 2)):
                recommendation = "LIKELY_TP"
                likely_tp.append((taxon, gene))

            # Ambiguous — needs review
            else:
                recommendation = "REVIEW"
                review.append((taxon, gene))

            assessments[(taxon, gene)] = {
                "taxon":          taxon,
                "gene":           gene,
                "trigger":        trigger,
                "prop_neg":       prop_neg,
                "posterior_var":  post_var,
                "n_flagged_genes_this_taxon": n_flagged,
                "n_cleared_genes_this_taxon": n_cleared,
                "n_qualified_cleared":        n_qualified,
                "is_orphan":                  is_orphan,
                "recommendation": recommendation,
            }

    exclude_list = likely_tp
    retain_list  = likely_fp + orphan_fp

    print(f"\n{'='*60}")
    print(f"FLAGGED SEQUENCE ASSESSMENT")
    print(f"{'='*60}")
    print(f"  Total flagged     : {len(assessments)}")
    print(f"  LIKELY_FP         : {len(likely_fp)}")
    print(f"  ORPHAN_FP         : {len(orphan_fp)}")
    print(f"  REVIEW            : {len(review)}")
    print(f"  LIKELY_TP         : {len(likely_tp)}")
    print(f"\n  Exclude from dataset : {len(exclude_list)}")
    print(f"  Retain despite flag  : {len(retain_list)}")

    if review:
        print(f"\n  Sequences requiring review:")
        for taxon, gene in sorted(review):
            a = assessments[(taxon, gene)]
            print(f"    {taxon:<40} {gene:<10} "
                  f"trigger={a['trigger']:<22} "
                  f"prop_neg={a['prop_neg']:.3f} "
                  f"flagged={a['n_flagged_genes_this_taxon']}/"
                  f"{a['n_flagged_genes_this_taxon']+a['n_cleared_genes_this_taxon']} genes")

    if likely_tp:
        print(f"\n  Likely true positives (recommend exclude):")
        for taxon, gene in sorted(likely_tp):
            a = assessments[(taxon, gene)]
            print(f"    {taxon:<40} {gene:<10} "
                  f"prop_neg={a['prop_neg']:.3f}")

    return {
        "assessments":  assessments,
        "likely_fp":    likely_fp,
        "orphan_fp":    orphan_fp,
        "review":       review,
        "likely_tp":    likely_tp,
        "exclude_list": exclude_list,
        "retain_list":  retain_list,
    }

from bioinf_packages.verify_funcs.batch_blast_verify import score_from_csv
from bioinf_packages.alignment_funcs._extract_accessions import extract_accessions

"""
results = score_from_csv(
    results_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/blast_results.csv",
    calibration_json = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/model_calibration.json",
    output_path = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/chimera_calls.csv",
    label_filter = "GENUINE",
)

lists = extract_verification_lists(result=results)

calibration_set = extract_accessions("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/callibration_accessions.csv")
print(calibration_set)
index = build_gbk_index(
    gbk_folder="C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
    glossary ="C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv",
)

qual = qualify_cleared_sequences(
    cleared_by_taxon       = lists["cleared_by_taxon"],
    gbk_folder             = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
    calibration_accessions = calibration_set,
    gbk_index = index,
    glossary = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"
)

assessment = assess_flagged_sequences(
    result     = results,
    qual       = qual,
    gbk_folder = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info",
    glossary   = "C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv",
    gbk_index  = index,
)

# Final dataset: qualified cleared + retained flagged
final_dataset = {
    taxon: genes
    for taxon, genes in qual["qualified"].items()
}
for taxon, gene in assessment["retain_list"]:
    final_dataset.setdefault(taxon, []).append(gene)

"""