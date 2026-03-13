import re
import requests

BASE_URL = "https://api.geneontology.org/api"
TIMEOUT = 10
HUMAN_TAXON = "NCBITaxon:9606"


def _clean_gene_symbol(gene: str) -> str:
    return re.sub(r"\s+", "", gene.strip())


def map_gene_to_uniprot(gene_symbol: str) -> str | None:
    
    #Map human gene symbol (e.g., HOXD13) to UniProtKB ID using GO autocomplete.
    #Returns UniProtKB ID (e.g., UniProtKB:P31271) or None.
    
    gene_symbol = _clean_gene_symbol(gene_symbol)

    url = f"{BASE_URL}/search/entity/autocomplete/{gene_symbol}"
    params = {"category": "gene"}

    r = requests.get(url, params=params, timeout=TIMEOUT)
    r.raise_for_status()
    data = r.json()

    for hit in data.get("docs", []):
        if (
            hit.get("taxon") == HUMAN_TAXON
            and hit.get("id", "").startswith("UniProtKB:")
            and hit.get("label", "").upper() == gene_symbol.upper()
        ):
            return hit["id"]

    return None


def fetch_go_function(uniprot_id: str) -> list[dict]:
    
    #Fetch GO molecular function annotations for given UniProtKB ID.
    #Returns list of dicts with GO ID and label.
    
    url = f"{BASE_URL}/bioentity/gene/{uniprot_id}/function"
    r = requests.get(url, timeout=TIMEOUT)
    r.raise_for_status()
    data = r.json()

    annotations = []
    for assoc in data.get("associations", []):
        obj = assoc.get("object", {})
        if "molecular_activity" in obj.get("category", []):
            annotations.append(
                {
                    "go_id": obj.get("id"),
                    "label": obj.get("label"),
                }
            )

    return annotations

def fetch_go_biological_process(uniprot_id: str) -> list[dict]:
    #Fetch GO biological processes annotations for given UniProtKB ID.
    url = f"{BASE_URL}/bioentity/gene/{uniprot_id}/function"
    r = requests.get(url, timeout=TIMEOUT)
    r.raise_for_status()
    data = r.json()

    annotations = []

    for assoc in data.get("associations", []):
        obj = assoc.get("object", {})
        if "biological_process" in obj.get("category", []):
            annotations.append({
                "go_id": obj.get("id"),
                "label": obj.get("label"),
            })

    return annotations
    

def fetch_go_cellular_component(uniprot_id: str) -> list[dict]:
    #Fetch GO cellular component annotations for given UniProtKB ID.
    url = f"{BASE_URL}/bioentity/gene/{uniprot_id}/function"
    r = requests.get(url, timeout=TIMEOUT)
    r.raise_for_status()
    data = r.json()

    annotations = []

    for assoc in data.get("associations", []):
        obj = assoc.get("object", {})
        if "cellular_component" in obj.get("category", []):
            annotations.append(
                {
                    "go_id": obj.get("id"),
                    "label": obj.get("label"),
                }
            )

    return annotations


def format_go_annotations(annotations: list[dict]) -> str:
    if not annotations:
        return "No GO molecular function annotations found."
    return "\n".join(f"{a['go_id']} | {a['label']}" for a in annotations)
