import pandas as pd
import re


GENE_REGEX = re.compile(r"^[A-Z0-9-]{2,10}$")


def extract_gene_list(df: pd.DataFrame, gene_column: str) -> list[str]:
    if gene_column not in df.columns:
        raise ValueError(f"Column '{gene_column}' not found")
    return (
        df[gene_column]
        .dropna()
        .astype(str)
        .str.upper()
        .unique()
        .tolist()
    )


def guess_gene_column(df_input: pd.DataFrame) -> list[dict]:
    """Rank columns likely to contain gene symbols (DB-agnostic)."""
    results = []
    for col in df_input.columns:
        if df_input[col].dtype != object:
            continue
        vals = (
            df_input[col]
            .dropna()
            .astype(str)
            .str.upper()
        )
        if vals.empty:
            continue
        gene_like = vals.str.match(GENE_REGEX)
        score = gene_like.mean()  # fraction gene-like
        if score == 0:
            continue
        results.append({
            "column": col,
            "gene_like_fraction": round(score, 3),
            "n_unique": vals.nunique()
        })
    return sorted(results, key=lambda x: x["gene_like_fraction"], reverse=True)

  
def genes_in_trrust(
    df_genes: pd.DataFrame,
    df_trrust: pd.DataFrame,
    gene_column: str
) -> list[str]:
    genes = extract_gene_list(df_genes, gene_column)
    trrust_set = set(
        df_trrust["TF"].str.upper()
    ).union(
        df_trrust["Target"].str.upper()
    )
    return sorted(set(genes) & trrust_set)
