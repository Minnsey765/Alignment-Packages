"""
slice_nesophontes.py
--------------------
Three functions for extracting Nesophontes gene sequences from a Yuan et al
NEXUS alignment, splitting the combined 12S/16S rRNA block via NCBI BLAST,
and appending the sequences into an existing per-gene alignment directory.
 
FUNCTIONS (call from __main__ block at the bottom, or import into your own script):
 
    slice_nexus(nexus_file, output_dir, taxon="nesophontes")
        Slices all 21 genes from the NEXUS alignment and writes per-gene FASTAs
        to output_dir, including a combined 12S_16S file for BLAST splitting.
 
    blast_split_rrna(combined_fasta, s12_fasta, s16_fasta, output_dir)
        Calls NCBI BLAST to locate 12S and 16S regions within the combined
        Nesophontes rRNA sequence and writes separate 12S_rRNA.fasta and
        16S_rRNA.fasta files to output_dir.
 
    append_to_alignment_dir(sliced_dir, alignment_dir, taxon_label="Nesophontes_sp",
                            accession="BraceEtAl")
        Appends the Nesophontes sequence from each per-gene FASTA in sliced_dir
        into the matching [gene].fasta file in alignment_dir, using the standard
        header format:  >Nesophontes_sp|GENE|orien:+|accession:BraceEtAl
 
Requirements: Python 3.6+, requests (pip install requests)
"""
 
import re
import sys
import time
import requests
from pathlib import Path
 
# ---------------------------------------------------------------------------
# Gene boundaries — 1-based inclusive — from Yuan et al MrBayes charset block
# ---------------------------------------------------------------------------
GENE_REGIONS = {
    "APOB":               [(1,     1179)],
    "BDNF":               [(2065,  2634)],
    "BRCA1":              [(2635,  3744)],
    "GHR":                [(7897,  8817)],
    "RAG1":               [(15202, 16269)],
    "RAG2":               [(16270, 16998)],
    "ATP6":               [(1180,  1860)],
    "ATP8":               [(1861,  2064)],
    "COX1":               [(3745,  5289)],
    "COX2":               [(5290,  5973)],
    "COX3":               [(5974,  6756)],
    "CYTB":               [(6757,  7896)],
    "ND1":                [(8818,  9771)],
    "ND2":                [(9772,  10818)],
    "ND3":                [(10819, 11163)],
    "ND4":                [(11164, 12543)],
    "ND4L":               [(12544, 12840)],
    "ND5":                [(12841, 14664)],
    "ND6":                [(14665, 15201)],
    # Combined — split by blast_split_rrna() below
    "12S_16S_rRNA_combined": [(16999, 18762)],
}
 
# Maps the combined key and per-gene keys to the filenames expected in the
# alignment directory.  Keys here must match filenames minus ".fasta".
GENE_TO_ALIGNMENT_FILENAME = {
    "APOB":    "APOB",
    "BDNF":    "BDNF",
    "BRCA1":   "BRCA1",
    "GHR":     "GHR",
    "RAG1":    "RAG1",
    "RAG2":    "RAG2",
    "ATP6":    "ATP6",
    "ATP8":    "ATP8",
    "COX1":    "COX1",
    "COX2":    "COX2",
    "COX3":    "COX3",
    "CYTB":    "CYTB",
    "ND1":     "ND1",
    "ND2":     "ND2",
    "ND3":     "ND3",
    "ND4":     "ND4",
    "ND4L":    "ND4L",
    "ND5":     "ND5",
    "ND6":     "ND6",
    "12S_rRNA": "12S_rRNA",
    "16S_rRNA": "16S_rRNA",
}
 
 
# ===========================================================================
# Internal helpers
# ===========================================================================
 
def _parse_nexus(filepath):
    """
    Parse a NEXUS file → dict {taxon_name: full_aligned_sequence}.
    Handles both interleaved and non-interleaved formats.
    """
    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
 
    matrix_match = re.search(
        r"matrix\s*(.*?)\s*;", text, re.IGNORECASE | re.DOTALL
    )
    if not matrix_match:
        raise ValueError(f"No MATRIX block found in {filepath}")
 
    sequences = {}
    for line in matrix_match.group(1).splitlines():
        line = line.strip()
        if not line or line.startswith("["):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        name = parts[0].strip("'\"")
        fragment = "".join(parts[1:])
        sequences[name] = sequences.get(name, "") + fragment
    return sequences
 
 
def _find_taxon(sequences, target):
    """Return the key in sequences that contains target (case-insensitive)."""
    matches = [k for k in sequences if target.lower() in k.lower()]
    if not matches:
        available = "\n  ".join(sorted(sequences.keys())[:30])
        raise KeyError(
            f"Taxon '{target}' not found in alignment.\n"
            f"First 30 taxa present:\n  {available}"
        )
    if len(matches) > 1:
        print(f"  [warning] Multiple taxa match '{target}': {matches}")
        print(f"  [warning] Using first: {matches[0]}")
    return matches[0]
 
 
def _slice(sequence, regions):
    """Extract and concatenate columns (1-based inclusive) from sequence."""
    return "".join(sequence[s - 1 : e] for s, e in regions)
 
 
def _gap_pct(seq):
    gaps = seq.count("-") + seq.count("?") + seq.count("N") + seq.count("n")
    return 100 * (len(seq) - gaps) / len(seq) if seq else 0.0
 
 
def _write_fasta(path, header, sequence, line_width=80):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(sequence), line_width):
            fh.write(sequence[i : i + line_width] + "\n")
 
 
def _read_fasta_first(path):
    """Return (header_without_gt, sequence) for the first record in a FASTA."""
    header, fragments = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    break           # only want the first record
                header = line[1:]
            elif header is not None:
                fragments.append(line)
    if header is None:
        raise ValueError(f"No FASTA records found in {path}")
    return header, "".join(fragments)
 
 
# ===========================================================================
# Function 1 — slice all genes from the NEXUS alignment
# ===========================================================================
 
def slice_nexus(nexus_file, output_dir, taxon="nesophontes"):
    """
    Extract per-gene sequences for `taxon` from a NEXUS alignment and write
    one FASTA per gene into output_dir.
 
    Parameters
    ----------
    nexus_file  : str or Path — path to the Yuan et al .nex file
    output_dir  : str or Path — directory to write per-gene FASTAs into
                  (created if it does not exist)
    taxon       : str — substring to match against taxon names (case-insensitive)
 
    Returns
    -------
    dict  {gene_name: sequence_string}  — the sliced sequences (gaps included)
    """
    nexus_file = Path(nexus_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
 
    print(f"\n[slice_nexus] Parsing {nexus_file.name} ...")
    sequences = _parse_nexus(nexus_file)
    print(f"  Found {len(sequences)} taxa in alignment.")
 
    taxon_key = _find_taxon(sequences, taxon)
    full_seq   = sequences[taxon_key]
    print(f"  Matched taxon: '{taxon_key}'  (alignment length: {len(full_seq)} bp)")
 
    if len(full_seq) < 18762:
        print(f"  [warning] Sequence shorter than expected 18762 cols — "
              f"check gene boundaries.")
 
    print(f"\n  {'Gene':<30} {'Columns':<22} {'bp':>6}  {'Coverage':>9}")
    print(f"  {'-'*30} {'-'*22} {'-'*6}  {'-'*9}")
 
    sliced = {}
    for gene, regions in GENE_REGIONS.items():
        seq  = _slice(full_seq, regions)
        pct  = _gap_pct(seq)
        cols = ", ".join(f"{s}-{e}" for s, e in regions)
        flag = "  *** mostly gaps" if pct < 20 else ""
        print(f"  {gene:<30} {cols:<22} {len(seq):>6}  {pct:>8.1f}%{flag}")
 
        out_path = output_dir / f"{gene}.fasta"
        _write_fasta(out_path, f"{taxon_key}_{gene}", seq)
        sliced[gene] = seq
 
    print(f"\n  Written {len(sliced)} FASTA files to '{output_dir}/'")
    print(f"  NOTE: '12S_16S_rRNA_combined.fasta' must be split — run blast_split_rrna() next.")
    return sliced
 
 
# ===========================================================================
# Function 2 — BLAST to split the combined 12S/16S block
# ===========================================================================
 
def blast_split_rrna(combined_fasta, s12_fasta, s16_fasta, output_dir):
    """
    Use NCBI BLAST (blastn, pairwise mode via the NCBI CGI API) to locate
    the 12S and 16S rRNA regions within the Nesophontes combined rRNA block,
    then write separate 12S_rRNA.fasta and 16S_rRNA.fasta files.
 
    Each Solenodon sequence is used as the QUERY; the Nesophontes combined
    block is the SUBJECT.  Hit coordinates on the subject tell us where each
    gene sits in the Nesophontes sequence.
 
    Parameters
    ----------
    combined_fasta : str or Path  FASTA with the Nesophontes 12S+16S block
    s12_fasta      : str or Path  FASTA with your Solenodon 12S rRNA sequence
    s16_fasta      : str or Path  FASTA with your Solenodon 16S rRNA sequence
    output_dir     : str or Path  where to write 12S_rRNA.fasta / 16S_rRNA.fasta
 
    Returns
    -------
    tuple (qs_12, qe_12, qs_16, qe_16)  1-based inclusive hit coords on the
    Nesophontes subject; None if no hit was found.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
 
    _, subject_seq = _read_fasta_first(combined_fasta)
    subject_clean  = re.sub(r"[-?]", "", subject_seq)
    print(f"\n[blast_split_rrna] Nesophontes combined rRNA: {len(subject_clean)} bp (gaps stripped)")
 
    BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
 
    def _submit(query_path, label):
        _, qseq = _read_fasta_first(query_path)
        qclean  = re.sub(r"[-?]", "", qseq)
        print(f"\n  [{label}] Query: {len(qclean)} bp  Subject: {len(subject_clean)} bp")
        print(f"  Submitting to NCBI BLAST ...")
 
        data = {
            "CMD":            "Put",
            "PROGRAM":        "blastn",
            "BLAST_PROGRAMS": "blastn",
            "QUERY":          qclean,
            "SUBJECTS":       subject_clean,
            "DATABASE":       "",
            "HITLIST_SIZE":   "5",
            "EXPECT":         "1e-3",
            "FORMAT_TYPE":    "XML",
            "WORD_SIZE":      "7",        # smaller word size = more sensitive
        }
        resp = requests.post(BLAST_URL, data=data, timeout=120)
        resp.raise_for_status()
 
        # ----------------------------------------------------------------
        # The NCBI CGI returns a page containing exactly one line of the form:
        #     RID = <11-char-alphanum>
        # followed later by:
        #     RTOE = <seconds>
        # We must match RID but NOT RTOE.
        # ----------------------------------------------------------------
        rid_match = re.search(r"^\s*RID = ([A-Z0-9]{6,})\s*$", resp.text, re.MULTILINE)
 
        if not rid_match:
            # Print full response for diagnosis
            print("\n  ===== FULL NCBI RESPONSE (for diagnosis) =====")
            print(resp.text[:3000])
            print("  ===== END RESPONSE =====\n")
            raise RuntimeError(
                f"[{label}] Could not extract RID from NCBI response.\n"
                f"  HTTP status : {resp.status_code}\n"
                f"  Query len   : {len(qclean)} bp\n"
                f"  Subject len : {len(subject_clean)} bp\n"
                f"  See response dump above for the actual error from NCBI."
            )
 
        rid  = rid_match.group(1)
        rtoe_match = re.search(r"^\s*RTOE = (\d+)\s*$", resp.text, re.MULTILINE)
        rtoe = int(rtoe_match.group(1)) if rtoe_match else 20
        print(f"  RID  = {rid}")
        print(f"  RTOE = {rtoe} s  (estimated wait time from NCBI)")
        return rid, rtoe, qclean
 
    def _poll(rid, rtoe, label):
        # Wait the suggested time before first poll
        wait = max(rtoe, 10)
        print(f"  Waiting {wait} s before first poll ...")
        time.sleep(wait)
 
        for attempt in range(40):
            poll = requests.get(BLAST_URL, params={
                "CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
            }, timeout=120)
            poll.raise_for_status()
            raw = poll.text
 
            # The status line is inside an HTML comment: <!-- QBlastInfoBegin ... -->
            status_match = re.search(r"Status=(\w+)", raw)
            status = status_match.group(1) if status_match else "UNKNOWN"
            print(f"  [{attempt+1:02d}] Status={status}")
 
            if status == "READY":
                if "<Hit>" not in raw:
                    print(f"  [warning] No hits found for {label}.")
                return raw
            if status == "FAILED":
                raise RuntimeError(f"BLAST FAILED for {label} (RID={rid})")
            if status == "UNKNOWN":
                # Dump a snippet to help diagnose
                snippet = raw[:400].replace("\n", " ")
                print(f"  [debug] Response snippet: {snippet}")
 
            time.sleep(15)
 
        raise TimeoutError(
            f"BLAST timed out for {label} (RID={rid}).\n"
            f"Check manually: https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&RID={rid}"
        )
 
    def _parse_hit_coords(xml, label):
        """Pull the best HSP hit coords on the SUBJECT (Nesophontes)."""
        hs = re.search(r"<Hsp_hit-from>(\d+)</Hsp_hit-from>",   xml)
        he = re.search(r"<Hsp_hit-to>(\d+)</Hsp_hit-to>",       xml)
        qs = re.search(r"<Hsp_query-from>(\d+)</Hsp_query-from>", xml)
        qe = re.search(r"<Hsp_query-to>(\d+)</Hsp_query-to>",   xml)
        if not hs or not he:
            print(f"  [warning] No hit coordinates in XML for {label}.")
            return None, None
        start = min(int(hs.group(1)), int(he.group(1)))
        end   = max(int(hs.group(1)), int(he.group(1)))
        qspan = f"{qs.group(1)}–{qe.group(1)}" if qs and qe else "?"
        print(f"  {label}: Nesophontes subject hit = {start}–{end}  "
              f"(Solenodon query covered: {qspan})")
        return start, end
 
    def _write_gene(seq, start, end, gene_label, out_path):
        if start is None:
            print(f"  [skip] {gene_label} — no hit, file not written.")
            return
        sliced = seq[start - 1 : end]
        _write_fasta(out_path, f"Nesophontes_sp_{gene_label}", sliced)
        print(f"  Wrote {gene_label}: {len(sliced)} bp → {out_path.name}")
 
    # ------------------------------------------------------------------
    # 12S
    # ------------------------------------------------------------------
    rid12, rtoe12, _ = _submit(s12_fasta, "12S rRNA")
    xml12 = _poll(rid12, rtoe12, "12S rRNA")
    qs_12, qe_12 = _parse_hit_coords(xml12, "12S rRNA")
 
    # ------------------------------------------------------------------
    # 16S
    # ------------------------------------------------------------------
    rid16, rtoe16, _ = _submit(s16_fasta, "16S rRNA")
    xml16 = _poll(rid16, rtoe16, "16S rRNA")
    qs_16, qe_16 = _parse_hit_coords(xml16, "16S rRNA")
 
    # ------------------------------------------------------------------
    # Write output FASTAs
    # ------------------------------------------------------------------
    print("\n  Writing output FASTAs ...")
    _write_gene(subject_clean, qs_12, qe_12, "12S_rRNA", output_dir / "12S_rRNA.fasta")
    _write_gene(subject_clean, qs_16, qe_16, "16S_rRNA", output_dir / "16S_rRNA.fasta")
 
    return qs_12, qe_12, qs_16, qe_16
 
 
# ===========================================================================
# Function 3 — append Nesophontes sequences into the alignment directory
# ===========================================================================
 
def append_to_alignment_dir(
    sliced_dir,
    alignment_dir,
    taxon_label="Nesophontes_sp",
    accession="BraceEtAl",
    orientation="+",
):
    """
    For each of the 21 genes, append the Nesophontes sequence into the
    corresponding [gene].fasta file in alignment_dir.
 
    Header format written:
        >Nesophontes_sp|GENE|orien:+|accession:BraceEtAl
 
    Parameters
    ----------
    sliced_dir    : str or Path — directory produced by slice_nexus() (and
                    blast_split_rrna()), containing files named e.g. CYTB.fasta,
                    12S_rRNA.fasta, 16S_rRNA.fasta, etc.
    alignment_dir : str or Path — directory containing your existing per-gene
                    FASTAs named exactly [gene].fasta
    taxon_label   : str — taxon part of the header (default: 'Nesophontes_sp')
    accession     : str — accession part of the header (default: 'BraceEtAl')
    orientation   : str — strand orientation for header (default: '+')
 
    Returns
    -------
    list of gene names successfully appended
    """
    sliced_dir    = Path(sliced_dir)
    alignment_dir = Path(alignment_dir)
 
    print(f"\n[append_to_alignment_dir]")
    print(f"  Sliced dir   : {sliced_dir}")
    print(f"  Alignment dir: {alignment_dir}")
 
    appended, skipped, missing = [], [], []
 
    for gene, aln_stem in GENE_TO_ALIGNMENT_FILENAME.items():
        # Where we expect to find the sliced sequence
        sliced_path = sliced_dir / f"{gene}.fasta"
        if not sliced_path.exists():
            # blast_split_rrna writes 12S_rRNA.fasta / 16S_rRNA.fasta directly
            sliced_path = sliced_dir / f"{aln_stem}.fasta"
        if not sliced_path.exists():
            print(f"  [skip]   {gene:<20} — sliced file not found in {sliced_dir}")
            missing.append(gene)
            continue
 
        # Target alignment file
        aln_path = alignment_dir / f"{aln_stem}.fasta"
        if not aln_path.exists():
            print(f"  [skip]   {gene:<20} — {aln_path.name} not found in alignment_dir")
            skipped.append(gene)
            continue
 
        # Read the sliced sequence
        try:
            _, seq = _read_fasta_first(sliced_path)
        except ValueError as e:
            print(f"  [error]  {gene:<20} — {e}")
            skipped.append(gene)
            continue
 
        # Strip alignment gaps — ClustalOmega requires unaligned input;
        # gaps will be reintroduced during alignment
        import re as _re
        seq_clean = _re.sub(r"[-?]", "", seq)
 
        # Skip if the sequence is empty or entirely gaps (no data for this gene)
        if not seq_clean:
            print(f"  [no data]{gene:<20} — sequence is entirely gaps, skipping "
                  f"(Nesophontes has no {gene} data in this alignment)")
            missing.append(gene)
            continue
 
        real_pct = 100 * len(seq_clean) / len(seq) if seq else 0
        if real_pct < 10:
            print(f"  [sparse] {gene:<20} — only {real_pct:.1f}% real bases after gap removal "
                  f"({len(seq_clean)} bp). Writing anyway — review before alignment.")
 
        # Report N content — Ns are valid IUPAC and kept, but flag high counts
        n_count = seq_clean.upper().count("N")
        n_pct   = 100 * n_count / len(seq_clean) if seq_clean else 0
        if n_pct > 30:
            print(f"  [Ns]     {gene:<20} — {n_pct:.1f}% ambiguous bases "
                  f"({n_count}/{len(seq_clean)} bp are N). Kept as-is for alignment.")
 
        # Build the standardised header
        header = f"{taxon_label}|{aln_stem}|orien:{orientation}|accession:{accession}"
 
        # Check the sequence isn't already in the file
        with open(aln_path) as fh:
            existing = fh.read()
        if taxon_label in existing:
            print(f"  [exists] {gene:<20} — {taxon_label} already present in {aln_path.name}, skipping.")
            skipped.append(gene)
            continue
 
        # Detect whether existing file uses continuous or wrapped sequences
        non_header_lines = [l.strip() for l in existing.splitlines()
                            if l.strip() and not l.startswith(">")]
        wrapped = bool(non_header_lines) and len(non_header_lines[0]) < 200
 
        # Append matching the existing format
        with open(aln_path, "a") as fh:
            fh.write(f"\n>{header}\n")
            if wrapped:
                for i in range(0, len(seq_clean), 80):
                    fh.write(seq_clean[i : i + 80] + "\n")
            else:
                fh.write(seq_clean + "\n")
 
        n_str = f", {n_pct:.1f}% N" if n_count > 0 else ""
        print(f"  [ok]     {gene:<20} → {aln_path.name}  "
              f"({len(seq_clean)} bp after gap removal{n_str})")
        appended.append(gene)
 
    print(f"\n  Appended : {len(appended)} genes")
    print(f"  Skipped  : {len(skipped)} genes (file not found or already present)")
    print(f"  Missing  : {len(missing)} genes (not in sliced_dir — run blast_split_rrna for rRNA)")
    return appended

def manual_split_rrna(combined_fasta, s12_start, s12_end, s16_start, s16_end, output_dir):
    """
    Manually split the combined rRNA FASTA using coordinates from BLAST.
    
    s12_start, s12_end : hit coords for 12S on the Nesophontes subject (1-based)
    s16_start, s16_end : hit coords for 16S on the Nesophontes subject (1-based)
    """
    import re
    from pathlib import Path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    _, seq = _read_fasta_first(combined_fasta)
    clean  = re.sub(r"[-?]", "", seq)
    print(f"Combined rRNA: {len(clean)} bp")

    for label, start, end, fname in [
        ("12S_rRNA", s12_start, s12_end, "12S_rRNA.fasta"),
        ("16S_rRNA", s16_start, s16_end, "16S_rRNA.fasta"),
    ]:
        sliced = clean[start - 1 : end]   # 1-based → Python slice
        out    = output_dir / fname
        _write_fasta(out, f"Nesophontes_sp_{label}", sliced)
        print(f"  {label}: {len(sliced)} bp → {out.name}")


# ===========================================================================
# __main__ — edit the paths here and run:  python slice_nesophontes.py
# ===========================================================================

if __name__ == "__main__":

    # ------------------------------------------------------------------
    # STEP 1 — Slice all genes from the Yuan et al NEXUS alignment
    # ------------------------------------------------------------------
    NEXUS_FILE  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/yuanetal_appendix/nexus/adunatorN_5Jcal_R4C8.nex"      # <-- path to your NEXUS file
    SLICED_DIR  = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/nesophontes"            # <-- where per-gene FASTAs are saved

#    slice_nexus(
#        nexus_file  = NEXUS_FILE,
#        output_dir  = SLICED_DIR,
#        taxon       = "nesophontes",        # case-insensitive substring match
#    )
#
#    # ------------------------------------------------------------------
#    # STEP 2 — Split the combined 12S/16S rRNA block using NCBI BLAST
#    #          (needs an internet connection; takes ~2 minutes)
#    # ------------------------------------------------------------------
#    blast_split_rrna(
#        combined_fasta = f"{SLICED_DIR}/12S_16S_rRNA_combined.fasta",
#        s12_fasta      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/nesophontes/Solenodon_12S_rRNA.fasta",   # <-- your Solenodon 12S FASTA
#        s16_fasta      = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/nesophontes/Solenodon_16S_rRNA.fasta",   # <-- your Solenodon 16S FASTA
#        output_dir     = SLICED_DIR,              # writes 12S_rRNA.fasta + 16S_rRNA.fasta here
#    )
#    manual_split_rrna(
#    combined_fasta = f"{SLICED_DIR}/12S_16S_rRNA_combined.fasta",
#    s12_start = 89,
#    s12_end   = 711,
#    s16_start = 1162,
#    s16_end   = 1740,
#    output_dir = SLICED_DIR,
#    )


    # ------------------------------------------------------------------
    # STEP 3 — Append Nesophontes sequences into your alignment directory
    # ------------------------------------------------------------------
    ALIGNMENT_DIR = "C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_data_clean/relaxed"         # <-- directory with your [gene].fasta files

    append_to_alignment_dir(
        sliced_dir    = SLICED_DIR,
        alignment_dir = ALIGNMENT_DIR,
        taxon_label   = "Nesophontes_sp",
        accession     = "BraceEtAl",
        orientation   = "+",
    )

