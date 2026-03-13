import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
SFARI_PATH = BASE_DIR / "data" / "SFARI.csv"
SFARI_COLUMNS = [
    "gene_symbol",
    "gene_name",
    "ensembl_id",
    "chromosome",
    "genetic_category",
    "gene_score",
    "syndromic",
    "eagle",
    "number_of_reports"]

def load_sfari() -> pd.DataFrame:
    df = pd.read_csv(
        SFARI_PATH,
        sep=",",
        header=0,
        names=SFARI_COLUMNS
    )
    df["gene_symbol_norm"] = (
        df["gene_symbol"]
        .astype(str)
        .str.upper()
    )
    return df

def get_sfari_genes_df(df_sfari: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    genes_df = pd.DataFrame({"gene_symbol_norm": [g.upper() for g in genes]})
    return df_sfari.merge(
        genes_df,
        on="gene_symbol_norm",
        how="inner"
    )