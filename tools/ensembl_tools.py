import re
import requests
import matplotlib.pyplot as plt
from Bio import Phylo
from io import StringIO
import os

BASE_URL = "https://rest.ensembl.org"
TIMEOUT = 60


def _clean_gene_symbol(gene: str) -> str:
    return re.sub(r"\s+", "", gene.strip())


def fetch_orthologs(species: str, gene_symbol: str) -> list[dict]:
    gene_symbol = _clean_gene_symbol(gene_symbol)
    url = f"{BASE_URL}/homology/symbol/{species}/{gene_symbol}"
    headers = {"Content-Type": "application/json"}

    r = requests.get(url, headers=headers, timeout=TIMEOUT)
    r.raise_for_status()
    data = r.json()

    results = []
    for entry in data.get("data", []):
        for homology in entry.get("homologies", []):
            target = homology.get("target", {})
            results.append(
                {
                    "species": target.get("species"),
                    "ensembl_id": target.get("id"),
                    "type": homology.get("type"),
                    "percent_identity": homology.get("target_percent_id"),
                }
            )
    return results

def save_tree_visualization(gene_symbol: str, species: str = "homo_sapiens", dpi: int = 200) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from Bio import Phylo
    from io import StringIO
    import re

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(BASE_DIR, "..", "plots")
    os.makedirs(out_dir, exist_ok=True)

    url = f"{BASE_URL}/genetree/member/symbol/{species}/{gene_symbol}"

    r_json = requests.get(url, params={
        "content-type": "application/json",
        "sequence": "none",
        "aligned": "0"
    }, timeout=TIMEOUT)
    r_json.raise_for_status()
    tree_json = r_json.json()

    id_to_label = {}

    def extract_labels(node):
        if not isinstance(node, dict):
            return
        node_id = node.get("id", "")
        if node_id and node_id.startswith("ENS"):
            seq = node.get("sequence", {})
            mol = seq.get("mol_seq", {}) if isinstance(seq, dict) else {}
            tax = node.get("taxonomy", {})
            sci = tax.get("scientific_name", "") if isinstance(tax, dict) else ""
            gene_name = ""
            for g in node.get("genes", []):
                if isinstance(g, dict):
                    gi = g.get("gene", {})
                    gene_name = gi.get("name", "") or gi.get("id", "")
                    sp_info = gi.get("species", {})
                    if isinstance(sp_info, dict):
                        sci = sp_info.get("scientific_name", sci)
                    break
            parts = sci.split()
            short_sp = f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) >= 2 else sci
            label = f"{gene_name} [{short_sp}]" if gene_name else f"[{short_sp}]"
            if label.strip("[ ]"):
                id_to_label[node_id] = label

        for child in node.get("children", []):
            extract_labels(child)

    extract_labels(tree_json)

    r_nwk = requests.get(url, params={
        "content-type": "text/x-nh",
        "sequence": "none",
        "aligned": "0"
    }, timeout=TIMEOUT)
    r_nwk.raise_for_status()
    newick_raw = r_nwk.text


    print(f"[DEBUG] id_to_label має {len(id_to_label)} записів")
    if id_to_label:
        sample = list(id_to_label.items())[:3]
        print(f"[DEBUG] Приклади: {sample}")

    ids_in_newick = re.findall(r'ENS[A-Z0-9]+', newick_raw)
    print(f"[DEBUG] ENS* у Newick: {ids_in_newick[:5]}")

    def replace_id(m):
        eid = m.group(0)
        label = id_to_label.get(eid, "")
        if not label:
            return eid 
        label = re.sub(r'[():,;]', '_', label)
        return label

    newick_labeled = re.sub(r'ENS[A-Z0-9]+', replace_id, newick_raw)

    tree = Phylo.read(StringIO(newick_labeled), "newick")

    CLADE_COLORS = {
        "homo":       "#9333ea", "pan":        "#a855f7",
        "primates":   "#6366f1", "macaca":     "#4f46e5",
        "mus":        "#16a34a", "rattus":     "#15803d",
        "rodentia":   "#166534", "cricet":     "#14532d",
        "canis":      "#b45309", "felidae":    "#92400e",
        "carnivora":  "#d97706", "bos":        "#ea580c",
        "sus":        "#c2410c", "pecora":     "#7c2d12",
        "cetacea":    "#0369a1", "salmo":      "#0284c7",
        "oryzias":    "#0ea5e9", "cichlidae":  "#38bdf8",
        "cyprin":     "#1d4ed8", "percomorph": "#2563eb",
        "gallus":     "#dc2626", "aves":       "#b91c1c",
        "xenopus":    "#ca8a04", "anura":      "#a16207",
        "drosophila": "#c026d3",
    }

    def get_color(name: str) -> str:
        low = (name or "").lower()
        for key, color in CLADE_COLORS.items():
            if key in low:
                return color
        return "#1e293b"

    def label_node(clade) -> str:
        name = clade.name or ""
        name = re.sub(r'_\d+$', '', name).replace("_", " ")
        return name[:40] + "…" if len(name) > 40 else name

    terminals = tree.count_terminals()
    fig_h = max(20, terminals * 0.22)
    fig, ax = plt.subplots(figsize=(26, fig_h))
    fig.patch.set_facecolor("#ffffff")
    ax.set_facecolor("#ffffff")

    Phylo.draw(
        tree, axes=ax, do_show=False,
        label_func=label_node,
        label_colors=lambda name: (
            "#9333ea" if "h. sapiens" in (name or "").lower()
                      or ("homo" in (name or "").lower() and "sapiens" in (name or "").lower())
            else get_color(name)
        ),
    )

    for line in ax.get_lines():
        line.set_color("#334155")
        line.set_linewidth(0.7)
        line.set_alpha(1.0)

    for text in ax.texts:
        txt = text.get_text()
        is_human = "h. sapiens" in txt.lower() or (
            "homo" in txt.lower() and "sapiens" in txt.lower()
        )
        color = "#9333ea" if is_human else get_color(txt)
        text.set_color(color)
        text.set_fontsize(9.5 if is_human else 7)
        text.set_fontfamily("monospace")
        text.set_fontweight("bold" if is_human else "normal")
        if is_human:
            text.set_bbox(dict(
                boxstyle="round,pad=0.25",
                facecolor="#f3e8ff",
                edgecolor="#9333ea",
                linewidth=0.8,
            ))

    ax.set_title(
        f"Філогенетичне дерево  ·  {gene_symbol}  ·  {species}  ·  Ensembl",
        color="#0f172a", fontsize=14, fontweight="bold",
        fontfamily="monospace", pad=16,
    )
    ax.axis("off")

    patches = [
        mpatches.Patch(color="#9333ea", label="Homo sapiens"),
        mpatches.Patch(color="#6366f1", label="Primates"),
        mpatches.Patch(color="#16a34a", label="Rodentia"),
        mpatches.Patch(color="#d97706", label="Carnivora"),
        mpatches.Patch(color="#ea580c", label="Artiodactyla"),
        mpatches.Patch(color="#0284c7", label="Риби"),
        mpatches.Patch(color="#dc2626", label="Aves"),
        mpatches.Patch(color="#ca8a04", label="Amphibia"),
        mpatches.Patch(color="#1e293b", label="Інші"),
    ]
    legend = ax.legend(
        handles=patches, loc="lower left", framealpha=0.9,
        facecolor="#f8fafc", edgecolor="#cbd5e1",
        labelcolor="#0f172a", fontsize=9,
        title="Таксономічна група", title_fontsize=9,
    )
    legend.get_title().set_color("#475569")

    out_path = os.path.join(out_dir, f"{gene_symbol}_tree.png")
    try:
        plt.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="#ffffff")
        return os.path.abspath(out_path)
    finally:
        plt.close()

def fetch_gene_tree(
    species: str,
    gene_symbol: str,
    target_species: list[str] | None = None,
) -> dict | None:
    gene_symbol = _clean_gene_symbol(gene_symbol)
    params = "?aligned=0&sequence=none"

    if target_species:
        species_str = ",".join(target_species)
        params += f"&prune_species={species_str}"

    url = f"{BASE_URL}/genetree/member/symbol/{species}/{gene_symbol}{params}"
    headers = {"Content-Type": "application/json"}

    r = requests.get(url, headers=headers, timeout=TIMEOUT)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


def format_orthologs(orthologs: list[dict]) -> str:
    if not orthologs:
        return "No orthologs found."

    lines = []
    for o in orthologs:
        lines.append(
            f"{o['species']} | {o['ensembl_id']} | "
            f"{o['type']} | identity={o['percent_identity']}%"
        )
    return "\n".join(lines)


def fetch_gene_coords(species: str, gene_symbol: str) -> dict | None:
    """Отримує координати гена для пошуку регуляторних зон."""
    gene_symbol = _clean_gene_symbol(gene_symbol)
    url = f"{BASE_URL}/lookup/symbol/{species}/{gene_symbol}"
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, headers=headers, timeout=TIMEOUT)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


def fetch_regulatory_features(species: str, chrom: str, start: int, end: int) -> list[dict]:
    """Шукає промотори, енхансери та сайти зв'язування у вказаному регіоні."""
    region = f"{chrom}:{max(1, start - 5000)}-{end + 5000}"
    url = f"{BASE_URL}/overlap/region/{species}/{region}?feature=regulatory"
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, headers=headers, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()

def find_regulatory_hits(input_df, reg_elements: list[dict]) -> list[dict]:
    """
    Порівнює мутації з файлу з регуляторними елементами Ensembl.

    """
    if "position" in input_df.columns:
        pos_col = "position"
    elif "start" in input_df.columns:
        pos_col = "start"
    else:
        raise ValueError(
            f"Колонка позиції не знайдена. Доступні колонки: {list(input_df.columns)}"

        )
    hits = []
    for _, row in input_df.iterrows():
        mutation_pos = int(row[pos_col])
        for re_elem in reg_elements:
            if re_elem["start"] <= mutation_pos <= re_elem["end"]:
                hits.append(
                    {
                        "mutation_pos": mutation_pos,
                        "ref": row.get("ref", "N/A"),
                        "alt": row.get("alt", "N/A"),
                        "reg_id": re_elem.get("id", "N/A"),
                        "reg_type": re_elem.get("feature_type", "Unknown"),
                        "reg_start": re_elem.get("start"),
                        "reg_end": re_elem.get("end"),
                        "distance_to_reg_start": mutation_pos - re_elem["start"],
                    }
                )
    return hits