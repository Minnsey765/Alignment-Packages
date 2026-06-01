from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def assess_peer_review_status(
    gbk_file: str | Path,
    require_pubmed: bool = False,
) -> Dict:
    """
    Heuristically assess whether a GenBank sequence appears to be
    supported by peer-reviewed literature.

    Parameters
    ----------
    gbk_file : str | Path
        Path to a .gbk / .gb / .genbank file.
    require_pubmed : bool
        If True, only counts PubMed-linked references as peer-reviewed.

    Returns
    -------
    dict
        {
            "accession": str,
            "peer_reviewed": bool,
            "confidence": str,
            "evidence": list[str],
            "references": list[dict],
        }
    """

    record: SeqRecord = SeqIO.read(str(gbk_file), "genbank")

    accession = record.id
    references = []
    evidence = []

    peer_reviewed = False
    confidence = "low"

    # GenBank stores references in annotations["references"]
    refs = record.annotations.get("references", [])

    if not refs:
        evidence.append("No references found in GenBank record.")
    else:
        for ref in refs:
            ref_info = {
                "title": getattr(ref, "title", None),
                "journal": getattr(ref, "journal", None),
                "pubmed_id": getattr(ref, "pubmed_id", None),
                "authors": getattr(ref, "authors", None),
            }

            references.append(ref_info)

            journal = (ref_info["journal"] or "").lower()
            pubmed_id = ref_info["pubmed_id"]

            # Detect direct submissions
            if "direct submission" in journal:
                evidence.append(
                    "Contains a direct submission reference "
                    "(not independently peer reviewed)."
                )
                continue

            # Detect unpublished sequences
            if "unpublished" in journal:
                evidence.append("Contains unpublished reference.")
                continue

            # PubMed-linked publication
            if pubmed_id:
                evidence.append(
                    f"Found PubMed-linked publication: PMID {pubmed_id}"
                )

                if require_pubmed:
                    peer_reviewed = True

            # Journal heuristic
            likely_journal = (
                journal
                and "submitted" not in journal
                and "patent" not in journal
                and "thesis" not in journal
            )

            if likely_journal:
                evidence.append(
                    f"Appears associated with journal publication: "
                    f"{ref_info['journal']}"
                )

                peer_reviewed = True

        # Confidence scoring
        pubmed_count = sum(
            1 for r in references if r.get("pubmed_id")
        )

        if peer_reviewed:
            if pubmed_count >= 1:
                confidence = "high"
            else:
                confidence = "medium"
        else:
            confidence = "low"

    return {
        "accession": accession,
        "peer_reviewed": peer_reviewed,
        "confidence": confidence,
        "evidence": evidence,
        "references": references,
    }
#print(assess_peer_review_status("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info/AB175120.gbk"))
#print(assess_peer_review_status("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info/NC_042734.gbk")["confidence"])


def loop_peer_review(
    gbk_directory: str | Path,
    require_pubmed: bool = False,
    recursive: bool = False,
) -> List[Dict]:
    """
    Return only GenBank entries with HIGH confidence peer-review support.

    Parameters
    ----------
    gbk_directory : str | Path
        Directory containing .gbk files.
    require_pubmed : bool
        Passed through to assess_peer_review_status().
    recursive : bool
        If True, search subdirectories recursively.

    Returns
    -------
    list[dict]
        [
            {
                "file": str,
                "accession": str,
                "peer_reviewed": bool,
                "confidence": "high",
            },
            ...
        ]
    """

    gbk_directory = Path(gbk_directory)

    if recursive:
        gbk_files = list(gbk_directory.rglob("*.gbk"))
    else:
        gbk_files = list(gbk_directory.glob("*.gbk"))

    high_confidence_results = []

    for gbk_file in gbk_files:
        try:
            result = assess_peer_review_status(
                gbk_file,
                require_pubmed=require_pubmed,
            )

            if result["confidence"] == "high":
                high_confidence_results.append({
                    #"file": str(gbk_file),
                    "accession": result["accession"],
                    "peer_reviewed": result["peer_reviewed"],
                    "confidence": result["confidence"],
                })

        except Exception:
            # Skip unreadable or malformed files
            continue

    return high_confidence_results

#print((loop_peer_review("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/aligment/code/fasta_info")))