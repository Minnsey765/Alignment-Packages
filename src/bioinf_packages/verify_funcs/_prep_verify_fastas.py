import csv
import glob
import os
from pathlib import Path
from collections import defaultdict

def prepare_verification_fasta(
        gene_fasta_dir: str,
        calibration_csv: str,
        output_fasta: str,
        glossary: dict,
        gene_symbol_map: dict = None,
        min_length: int = 100,
        exclude_mode: str = "accession_and_gene",
) -> dict:
    """
    Combine per-gene FASTA files into a single FASTA for verification,
    excluding any sequences that were used in calibration.

    Reads all [gene_symbol].fasta files from gene_fasta_dir, looks up
    which (accession, gene) pairs appear in the calibration CSV, and
    writes all remaining sequences to a single output FASTA with
    pipe-delimited headers matching the format expected by
    prepare_subsample_fastas().

    Expected input FASTA header format (flexible):
        >accession [any description]
        >accession.version [any description]
        >genus_species_accession [any description]
    The accession is extracted from the first whitespace-delimited token.

    Output FASTA header format (for prepare_subsample_fastas):
        >Genus_species|gene|orien:+|accession:XXXX

    Parameters
    ----------
    gene_fasta_dir  : str   directory containing per-gene FASTA files
                            named [gene_symbol].fasta or
                            [gene_symbol].fa. Gene symbol is taken
                            from the filename stem.
    calibration_csv : str   path to blast_results.csv used for
                            calibration. Any (accession, gene) pair
                            appearing here is excluded from output.
    output_fasta    : str   path to write the combined output FASTA.
    glossary        : dict  from load_glossary(). Used to normalise
                            gene symbols from filenames to canonical
                            names (e.g. "cytb" -> "CYTB").
    gene_symbol_map : dict  optional manual overrides for gene symbol
                            normalisation, e.g. {"cox1": "COX1"}.
                            Applied before glossary lookup.
    min_length      : int   minimum sequence length to include.
                            Default 100bp.
    exclude_mode    : str   "accession_and_gene" — exclude only the
                            specific (accession, gene) combination
                            that appears in the calibration CSV.
                            "accession_only" — exclude any sequence
                            whose accession appears in the calibration
                            CSV regardless of gene.
                            Default "accession_and_gene" since you may
                            want to verify non-calibrated genes from
                            the same accession.

    Returns
    -------
    dict with:
        n_input_seqs     : int   total sequences read across all files
        n_excluded       : int   sequences excluded (in calibration set)
        n_too_short      : int   sequences excluded (too short)
        n_written        : int   sequences written to output FASTA
        excluded_by_gene : dict  gene -> count excluded
        written_by_gene  : dict  gene -> count written
        missing_taxon    : list  accessions where taxon could not be
                                 parsed from header
        output_fasta     : str
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    # ── Load calibration exclusion set ────────────────────────────────────────
    calibration_pairs    = set()   # (accession, gene) pairs to exclude
    calibration_accessions = set() # accessions to exclude (for accession_only)

    with open(calibration_csv, "r", newline="",
              encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc  = row.get("accession", "").strip()
            gene = row.get("gene", "").strip()
            if acc:
                calibration_accessions.add(acc)
            if acc and gene:
                calibration_pairs.add((acc, gene))

    print(f"\n{'='*60}")
    print(f"PREPARE VERIFICATION FASTA")
    print(f"{'='*60}")
    print(f"  Gene FASTA dir    : {gene_fasta_dir}")
    print(f"  Calibration CSV   : {calibration_csv}")
    print(f"  Output FASTA      : {output_fasta}")
    print(f"  Exclude mode      : {exclude_mode}")
    print(f"  Calibration pairs : {len(calibration_pairs)} "
          f"({len(calibration_accessions)} unique accessions)")

    # ── Find all gene FASTA files ─────────────────────────────────────────────
    fasta_files = []
    for ext in ("*.fasta", "*.fa", "*.fas", "*.fna"):
        fasta_files.extend(
            glob.glob(os.path.join(gene_fasta_dir, ext))
        )
    fasta_files = sorted(set(fasta_files))

    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files found in {gene_fasta_dir}"
        )

    print(f"  FASTA files found : {len(fasta_files)}")

    # ── Helper: normalise gene symbol from filename ───────────────────────────
    def normalise_gene_symbol(raw_symbol: str) -> str:
        """
        Normalise a gene symbol from a FASTA filename to canonical
        form using the glossary, e.g. "cytb" -> "CYTB",
        "cox1" -> "COX1", "12s_rrna" -> "12S_rRNA".
        """
        # Apply manual overrides first
        if gene_symbol_map:
            override = gene_symbol_map.get(raw_symbol.lower())
            if override:
                return override

        # Try glossary lookup — glossary maps raw -> canonical
        upper = raw_symbol.upper()
        for raw, canon in glossary.items():
            if isinstance(canon, str):
                if raw.upper() == upper:
                    return canon
            elif isinstance(canon, dict):
                if raw.upper() == upper:
                    return canon.get("canonical", raw_symbol)

        # Fall back to uppercase
        return raw_symbol.upper()

    # ── Helper: extract accession from FASTA header ───────────────────────────
    def extract_accession(header: str) -> str:
        """
        Extract accession from a FASTA record ID.
        Handles formats:
            KF783046        → KF783046
            KF783046.1      → KF783046   (strip version)
            gb|KF783046.1|  → KF783046
        """
        token = header.strip().split()[0]
        # Strip pipe-delimited accession format
        if "|" in token:
            parts = token.split("|")
            # Find the part that looks like an accession
            for p in parts:
                p = p.strip()
                if p and not p.lower() in ("gb", "ref", "emb",
                                            "dbj", "pir", "sp"):
                    token = p
                    break
        # Strip version number
        if "." in token:
            token = token.split(".")[0]
        return token.strip()

    # ── Helper: extract taxon from FASTA header ───────────────────────────────
    def extract_taxon(record_id: str, description: str) -> str:
        """
        Extract taxon (Genus_species) from a FASTA header.

        Tries several common formats:
        1. Pipe-delimited: Genus_species|gene|...
        2. Description contains "[organism=Genus species]"
        3. Description starts with Genus species
        4. Falls back to first two words of description capitalised
        """
        # Format 1: already pipe-delimited (our own format)
        if "|" in record_id:
            parts = record_id.split("|")
            if parts[0] and "_" in parts[0]:
                return parts[0].strip()

        # Format 2: [organism=Genus species]
        import re
        m = re.search(r'\[organism=([^\]]+)\]', description,
                      re.IGNORECASE)
        if m:
            org = m.group(1).strip().replace(" ", "_")
            return org

        # Format 3: description starts with two capitalised words
        desc_words = description.strip().split()
        if (len(desc_words) >= 2
                and desc_words[0][0].isupper()
                and desc_words[1][0].islower()):
            return f"{desc_words[0]}_{desc_words[1]}"

        return "Unknown_sp"

    # ── Process each gene FASTA ───────────────────────────────────────────────
    n_input      = 0
    n_excluded   = 0
    n_too_short  = 0
    n_written    = 0
    missing_taxon = []
    excluded_by_gene = defaultdict(int)
    written_by_gene  = defaultdict(int)

    output_records = []

    for fasta_file in fasta_files:
        filename    = Path(fasta_file).stem
        gene_symbol = normalise_gene_symbol(filename)

        records = list(SeqIO.parse(fasta_file, "fasta"))
        print(f"\n  {filename}.fasta → gene={gene_symbol} "
              f"({len(records)} sequences)")

        for record in records:
            n_input += 1

            accession = extract_accession(record.id)
            seq       = str(record.seq).replace("-", "").replace(
                "N", "").replace("n", "")

            # Length filter
            if len(seq) < min_length:
                n_too_short += 1
                continue

            # Exclusion check
            pair = (accession, gene_symbol)
            if exclude_mode == "accession_and_gene":
                if pair in calibration_pairs:
                    n_excluded += 1
                    excluded_by_gene[gene_symbol] += 1
                    continue
            elif exclude_mode == "accession_only":
                if accession in calibration_accessions:
                    n_excluded += 1
                    excluded_by_gene[gene_symbol] += 1
                    continue

            # Extract taxon
            taxon = extract_taxon(
                record.id, record.description
            )
            if taxon == "Unknown_sp":
                missing_taxon.append(
                    f"{accession} ({gene_symbol})"
                )

            # Build output header in prepare_subsample_fastas format
            header = (
                f"{taxon}|{gene_symbol}|orien:+|"
                f"accession:{accession}"
            )

            output_records.append(
                SeqRecord(
                    Seq(str(record.seq).replace("-", "")),
                    id          = header,
                    description = "",
                )
            )
            n_written += 1
            written_by_gene[gene_symbol] += 1

    # ── Write combined output FASTA ───────────────────────────────────────────
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(output_records, output_fasta, "fasta")

    # ── Print summary ─────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"  Total input sequences : {n_input}")
    print(f"  Excluded (calibration): {n_excluded}")
    print(f"  Excluded (too short)  : {n_too_short}")
    print(f"  Written to output     : {n_written}")
    print(f"  Output FASTA          : {output_fasta}")

    print(f"\n  Written per gene:")
    for gene, count in sorted(written_by_gene.items()):
        excl = excluded_by_gene.get(gene, 0)
        print(f"    {gene:<12} written={count:<5} "
              f"excluded={excl}")

    if missing_taxon:
        print(f"\n  WARNING: {len(missing_taxon)} sequences where "
              f"taxon could not be parsed from header:")
        for m in missing_taxon[:10]:
            print(f"    {m}")
        if len(missing_taxon) > 10:
            print(f"    ... and {len(missing_taxon)-10} more")
        print(f"  These will use taxon=Unknown_sp which will score "
              f"poorly in species_parser(). Consider fixing headers "
              f"manually or providing an accession->taxon lookup.")

    return {
        "n_input_seqs":     n_input,
        "n_excluded":       n_excluded,
        "n_too_short":      n_too_short,
        "n_written":        n_written,
        "excluded_by_gene": dict(excluded_by_gene),
        "written_by_gene":  dict(written_by_gene),
        "missing_taxon":    missing_taxon,
        "output_fasta":     output_fasta,
    }

from bioinf_packages.verify_funcs._gene_parser import load_glossary
prepare_verification_fasta(
    gene_fasta_dir = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data",
    calibration_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/blast_results.csv",
    output_fasta = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/all_queries.fasta",
    glossary = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"),
)