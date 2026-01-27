import pandas as pd

def extract_accessions(csv_path: str) -> list[str]:
    """
    Reads a CSV containing accession numbers under species (rows) and genes (columns),
    and returns a unique list of accession numbers.

    Args:
        csv_path (str): Path to the input CSV file.

    Returns:
        list[str]: Unique accession numbers.
    """
    # Load the table
    df = pd.read_csv(csv_path)

    # Drop the first column if it’s species names
    df_no_species = df.drop(columns=df.columns[0], errors="ignore")

    # Flatten table into a single list
    accessions = df_no_species.values.flatten()

    # Remove NaN/empty cells, entries with spaces, and strip whitespace
    accessions = [
        str(acc).strip()
        for acc in accessions
        if pd.notna(acc) and str(acc).strip() != "" and " " not in str(acc)
    ]

    # Remove duplicates while preserving order
    unique_accessions = list(dict.fromkeys(accessions))

    return unique_accessions

#print((extract_accessions("C:/Users/ojmin/OneDrive/Documents/UNI/MPhil/Project/newGenBank.csv")))