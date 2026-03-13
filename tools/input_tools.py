import pandas as pd
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent
INPUT_DIR = DATA_DIR / "data" / "input"


def list_input_files() -> list[str]:
    if not INPUT_DIR.exists():
        return []
    return sorted(p.name for p in INPUT_DIR.iterdir() if p.is_file())


def load_table(filename: str) -> pd.DataFrame:
    path = INPUT_DIR / filename
    if not path.exists():
        raise ValueError(f"{filename} not found in input directory")

    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in [".xls", ".xlsx"]:
        return pd.read_excel(path, engine="openpyxl")

    raise ValueError("Unsupported file format")
