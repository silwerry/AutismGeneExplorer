import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
TRRUST_PATH = BASE_DIR / "data" / "trrust_rawdata.human.tsv"

COLUMNS = ["TF", "Target", "Mode", "Source"]


def load_trrust() -> pd.DataFrame:
    df = pd.read_csv(
        TRRUST_PATH,
        sep="\t",
        header=None,
        names=COLUMNS
    )
    df["TF_norm"] = df["TF"].astype(str).str.upper()
    df["Target_norm"] = df["Target"].astype(str).str.upper()
    return df

def get_regulators_df(df_trrust: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    genes_df = pd.DataFrame({"Target_norm": [g.upper() for g in genes]})
    return df_trrust.merge(
        genes_df,
        on="Target_norm",
        how="inner"
    )


def get_targets_df(df_trrust: pd.DataFrame, tfs: list[str]) -> pd.DataFrame:
    tfs_df = pd.DataFrame({"TF_norm": [t.upper() for t in tfs]})
    return df_trrust.merge(
        tfs_df,
        on="TF_norm",
        how="inner"
    )

