import csv
import glob
import os
from pathlib import Path
from collections import defaultdict




import csv
import glob
import os
import re
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
    writes all remaining sequences to a single output FASTA.

    Expected input FASTA header formats (all handled):
        >Genus_species|GENE|orien:+|accession:GU981015
        >GU981015 Genus species gene description
        >GU981015.1 [organism=Genus species]
        >gb|GU981015.1| description

    Output FASTA header format:
        >Genus_species|GENE|orien:+|accession:GU981015

    Parameters
    ----------
    gene_fasta_dir  : str   directory containing per-gene FASTA files
                            named [gene_symbol].fasta or .fa
    calibration_csv : str   path to blast_results.csv used for
                            calibration
    output_fasta    : str   path to write the combined output FASTA
    glossary        : dict  from load_glossary()
    gene_symbol_map : dict  optional manual gene symbol overrides
    min_length      : int   minimum sequence length to include
    exclude_mode    : str   "accession_and_gene" or "accession_only"

    Returns
    -------
    dict with n_input_seqs, n_excluded, n_too_short, n_written,
              excluded_by_gene, written_by_gene, missing_taxon,
              output_fasta
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    # ── Load calibration exclusion set ────────────────────────────────────────
    calibration_pairs      = set()
    calibration_accessions = set()

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
        # Handle rRNA genes explicitly before any case conversion —
        # these must preserve mixed case to match the CSV and headers.
        # Filenames may be "12s_rrna", "12S_rRNA", "12s_rRNA" etc.
        rrna_map = {
            "12s_rrna": "12S_rRNA",
            "16s_rrna": "16S_rRNA",
            "12s":      "12S_rRNA",
            "16s":      "16S_rRNA",
            "12s_rrna": "12S_rRNA",
            "16s_rrna": "16S_rRNA",
        }
        lower = raw_symbol.lower()
        if lower in rrna_map:
            return rrna_map[lower]

        # Apply manual overrides
        if gene_symbol_map:
            override = gene_symbol_map.get(lower)
            if override:
                return override

        # Try glossary lookup
        upper = raw_symbol.upper()
        for raw, canon in glossary.items():
            if isinstance(canon, str):
                if raw.upper() == upper:
                    return canon
            elif isinstance(canon, dict):
                if raw.upper() == upper:
                    return canon.get("canonical", raw_symbol)

        # Fall back to uppercase — note this will incorrectly case
        # any rRNA genes not caught above, so the rrna_map above
        # should be exhaustive for your gene set.
        return upper


    # ── Helper: parse header fields ───────────────────────────────────────────
    def parse_header(record_id: str,
                     description: str) -> tuple:
        """
        Parse (taxon, accession) from a FASTA record ID and
        description.

        Handles two primary formats:

        Format A — pipe-delimited (our own format from previous runs):
            Genus_species|GENE|orien:+|accession:GU981015
            Fields are split on '|' and the accession is read from
            the field starting with 'accession:'.

        Format B — GenBank download format:
            GU981015.1 Genus species mitochondrion gene description
            [organism=Genus species]
            Accession is the first token (before any dot for version),
            taxon is extracted from [organism=...] or from the first
            two words of the description.

        Returns
        -------
        (taxon, accession) both as strings, or (None, None) if the
        header cannot be parsed.
        """
        record_id   = record_id.strip()
        description = description.strip()

        # ── Format A: pipe-delimited ──────────────────────────────────────────
        if '|' in record_id:
            fields  = record_id.split('|')
            taxon   = fields[0].strip()   # Genus_species
            acc     = None

            for field in fields:
                field = field.strip()
                if field.lower().startswith('accession:'):
                    acc = field.split(':', 1)[1].strip()
                    break

            # Validate: taxon should look like Genus_species
            # (starts with capital, has underscore, no digits)
            if (taxon
                    and '_' in taxon
                    and taxon[0].isupper()
                    and not any(c.isdigit() for c in taxon)):
                # taxon looks valid
                pass
            else:
                taxon = None

            if taxon and acc:
                return taxon, acc

            # If accession field is missing, taxon may itself be
            # an accession — fall through to Format B
            if not acc and taxon:
                # Check if first field is actually an accession
                # (contains digits, typical accession pattern)
                if any(c.isdigit() for c in taxon):
                    # It's an accession in the first field position
                    acc   = taxon.split('.')[0]
                    taxon = _extract_taxon_from_description(
                        description
                    )
                    return taxon, acc

        # ── Format B: GenBank / standard accession-first format ───────────────
        # Accession is the first whitespace-delimited token
        first_token = record_id.split()[0]

        # Handle pipe-wrapped accessions: gb|GU981015.1|
        if '|' in first_token:
            parts = first_token.split('|')
            for p in parts:
                p = p.strip()
                if (p
                        and p.lower() not in
                        ('gb', 'ref', 'emb', 'dbj', 'pir', 'sp',
                         'tr', 'up')
                        and p):
                    first_token = p
                    break

        # Strip version number (GU981015.1 -> GU981015)
        acc = first_token.split('.')[0].strip()

        # Taxon from description
        taxon = _extract_taxon_from_description(description)

        return taxon, acc

    def _extract_taxon_from_description(description: str) -> str:
        """
        Extract Genus_species from a sequence description string.

        Tries in order:
        1. [organism=Genus species] qualifier
        2. First two words where word 1 is capitalised and
           word 2 is lowercase
        3. Returns 'Unknown_sp' as fallback
        """
        # [organism=Genus species]
        m = re.search(r'\[organism=([^\]]+)\]',
                      description, re.IGNORECASE)
        if m:
            org   = m.group(1).strip()
            parts = org.split()
            if len(parts) >= 2:
                return f"{parts[0]}_{parts[1]}"
            return org.replace(' ', '_')

        # First two capitalised/lowercase words
        words = description.strip().split()
        if (len(words) >= 2
                and words[0]
                and words[0][0].isupper()
                and words[1]
                and words[1][0].islower()
                and words[1].isalpha()):
            return f"{words[0]}_{words[1]}"

        return "Unknown_sp"

    # ── Process each gene FASTA ───────────────────────────────────────────────
    n_input           = 0
    n_excluded        = 0
    n_too_short       = 0
    n_written         = 0
    missing_taxon     = []
    excluded_by_gene  = defaultdict(int)
    written_by_gene   = defaultdict(int)
    output_records    = []

    for fasta_file in fasta_files:
        filename    = Path(fasta_file).stem
        gene_symbol = normalise_gene_symbol(filename)

        records = list(SeqIO.parse(fasta_file, "fasta"))
        print(f"\n  {filename}.fasta → gene={gene_symbol} "
              f"({len(records)} sequences)")

        n_gene_excluded = 0
        n_gene_written  = 0

        for record in records:
            n_input += 1

            # Parse taxon and accession from header
            taxon, accession = parse_header(
                record.id, record.description
            )

            if accession is None:
                missing_taxon.append(
                    f"{record.id} ({gene_symbol}) — "
                    f"could not parse accession"
                )
                n_too_short += 1
                continue

            if taxon is None or taxon == "Unknown_sp":
                missing_taxon.append(
                    f"{accession} ({gene_symbol}) — "
                    f"could not parse taxon"
                )
                taxon = "Unknown_sp"

            # Strip gaps from sequence for length check
            seq_clean = (str(record.seq)
                         .replace("-", "")
                         .replace("N", "")
                         .replace("n", ""))

            if len(seq_clean) < min_length:
                n_too_short += 1
                continue

            # Exclusion check
            if exclude_mode == "accession_and_gene":
                if (accession, gene_symbol) in calibration_pairs:
                    n_excluded += 1
                    n_gene_excluded += 1
                    excluded_by_gene[gene_symbol] += 1
                    continue
            elif exclude_mode == "accession_only":
                if accession in calibration_accessions:
                    n_excluded += 1
                    n_gene_excluded += 1
                    excluded_by_gene[gene_symbol] += 1
                    continue

            # Build output header
            header = (
                f"{taxon}|{gene_symbol}|orien:+|"
                f"accession:{accession}"
            )

            output_records.append(
                SeqRecord(
                    Seq(str(record.seq).replace("-", "")),
                    id          = header,
                    name        = "",
                    description = "",
                )
            )
            n_written       += 1
            n_gene_written  += 1
            written_by_gene[gene_symbol] += 1

        print(f"    written={n_gene_written}  "
              f"excluded={n_gene_excluded}")

    # ── Write output FASTA ────────────────────────────────────────────────────
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
        print(f"    {gene:<15} written={count:<5} excluded={excl}")

    if missing_taxon:
        print(f"\n  WARNING: {len(missing_taxon)} sequences with "
              f"parsing issues:")
        for m in missing_taxon[:10]:
            print(f"    {m}")
        if len(missing_taxon) > 10:
            print(f"    ... and {len(missing_taxon)-10} more")

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
    gene_fasta_dir  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data",
    calibration_csv = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/callibration/blast_results.csv",
    output_fasta    = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/verify_data/query_verify/all_queries.fasta",
    glossary        = load_glossary("C:/Users/ojmin/OneDrive/Documents/UNI/Python_Packages/src/bioinf_packages/dictionary_funcs/glossary.csv"),
)