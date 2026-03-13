import os
import sys
import argparse
import requests
from io import StringIO

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
except ImportError:
    print("❌ Встановіть matplotlib: pip install matplotlib")
    sys.exit(1)

try:
    from Bio import Phylo
except ImportError:
    print("❌ Встановіть biopython: pip install biopython")
    sys.exit(1)

BASE_URL = "https://rest.ensembl.org"
TIMEOUT  = 30

CLADE_COLORS = {
    "homo":          "#e879f9",
    "pan":           "#c084fc",
    "hominidae":     "#a78bfa",
    "primates":      "#818cf8",
    "simiiformes":   "#6366f1",
    "catarrhini":    "#7c3aed",
    "rodentia":      "#34d399",
    "mus":           "#6ee7b7",
    "carnivora":     "#f59e0b",
    "canis":         "#fbbf24",
    "felidae":       "#fde68a",
    "salmo":         "#38bdf8",
    "oncorhynchus":  "#7dd3fc",
    "danio":         "#bae6fd",
    "oryzias":       "#93c5fd",
    "cichlidae":     "#60a5fa",
    "percomorph":    "#3b82f6",
    "cyprin":        "#2563eb",
    "drosophila":    "#fb923c",
    "caenorhabditis":"#a3e635",
    "arabidopsis":   "#86efac",
    "xenopus":       "#fcd34d",
    "gallus":        "#fca5a5",
}
DEFAULT_COLOR    = "#94a3b8"
HIGHLIGHT_TARGET = "Homo sapiens"


# ── Допоміжні функції ─────────────────────────────────────────

def get_color(name: str) -> str:
    low = name.lower()
    for key, color in CLADE_COLORS.items():
        if key in low:
            return color
    return DEFAULT_COLOR


def fetch_newick(gene: str, species: str) -> str | None:
    url = f"{BASE_URL}/genetree/member/symbol/{species}/{gene}"
    params = {"content-type": "text/x-nh", "sequence": "none", "aligned": "0"}
    try:
        r = requests.get(url, params=params, timeout=TIMEOUT)
        r.raise_for_status()
        return r.text
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            print(f"❌ Ген '{gene}' не знайдено у виді '{species}'.")
            print("   Перевір назву гену або спробуй інший вид (--species).")
        else:
            print(f"❌ HTTP помилка: {e}")
        return None
    except Exception as e:
        print(f"❌ Помилка запиту: {e}")
        return None


def label_node(clade) -> str:
    name = clade.name or ""
    return name[:28] + "…" if len(name) > 30 else name


def draw_tree(tree, gene: str, species: str, out_path: str, dpi: int):
    terminals = tree.count_terminals()
    fig_h = max(14, terminals * 0.18)
    fig_w = 22

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("#04080f")
    ax.set_facecolor("#04080f")

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=label_node,
        label_colors=lambda name: (
            "#e879f9" if HIGHLIGHT_TARGET.lower() in (name or "").lower()
            else get_color(name or "")
        ),
    )

    for line in ax.get_lines():
        line.set_color("#1e3a5c")
        line.set_linewidth(0.7)

    for text in ax.texts:
        txt = text.get_text()
        is_target = HIGHLIGHT_TARGET.lower() in txt.lower()
        text.set_color("#e879f9" if is_target else get_color(txt))
        text.set_fontsize(9 if is_target else 6.5)
        text.set_fontfamily("monospace")
        text.set_fontweight("bold" if is_target else "normal")
        if is_target:
            text.set_bbox(dict(
                boxstyle="round,pad=0.3",
                facecolor="#2d0a3e",
                edgecolor="#e879f9",
                linewidth=0.8,
            ))

    ax.set_title(
        f"Філогенетичне дерево гену  {gene}  ·  Ensembl REST API  ·  {species}",
        color="#e2e8f0", fontsize=13, fontweight="bold",
        fontfamily="monospace", pad=14,
    )
    ax.axis("off")

    patches = [
        mpatches.Patch(color="#e879f9", label="Homo sapiens"),
        mpatches.Patch(color="#818cf8", label="Primates"),
        mpatches.Patch(color="#34d399", label="Rodentia"),
        mpatches.Patch(color="#f59e0b", label="Carnivora"),
        mpatches.Patch(color="#38bdf8", label="Риби (Osteichthyes)"),
        mpatches.Patch(color="#fb923c", label="Drosophila"),
        mpatches.Patch(color="#fcd34d", label="Xenopus"),
        mpatches.Patch(color="#fca5a5", label="Gallus"),
        mpatches.Patch(color="#94a3b8", label="Інші"),
    ]
    legend = ax.legend(
        handles=patches, loc="lower left",
        framealpha=0.15, facecolor="#0d1b2e", edgecolor="#1e3a5c",
        labelcolor="#e2e8f0", fontsize=8,
        title="Таксономічна група", title_fontsize=8,
    )
    legend.get_title().set_color("#94a3b8")

    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="#04080f")
    plt.close()
    print(f"✅ Збережено: {out_path}")
    print(f"   Листків у дереві: {terminals}")



def main():
    parser = argparse.ArgumentParser(
        description="Філогенетичне дерево гену з Ensembl",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Приклади:
  python gene_tree_viz.py SHANK3
  python gene_tree_viz.py TP53
  python gene_tree_viz.py BRCA1 --dpi 300
  python gene_tree_viz.py MECP2 --species mus_musculus
  python gene_tree_viz.py PTEN --out /home/user/results
        """,
    )
    parser.add_argument("gene",
        help="Символ гену (наприклад: SHANK3, TP53, BRCA1)")
    parser.add_argument("--species", default="homo_sapiens",
        help="Вид (default: homo_sapiens)")
    parser.add_argument("--dpi", type=int, default=200,
        help="Роздільна здатність PNG (default: 200)")
    parser.add_argument("--out",
        help="Папка для збереження (default: plots/ поруч зі скриптом)")

    args = parser.parse_args()

    gene    = args.gene.strip().upper()
    species = args.species.strip().lower()
    out_dir = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{gene}_tree.png")

    print(f"⏳ Завантажую дерево для {gene} ({species})…")
    newick = fetch_newick(gene, species)
    if not newick:
        sys.exit(1)

    print("🌿 Парсинг Newick…")
    try:
        tree = Phylo.read(StringIO(newick), "newick")
    except Exception as e:
        print(f"❌ Помилка парсингу: {e}")
        sys.exit(1)

    print("🎨 Малюю дерево…")
    draw_tree(tree, gene, species, out_path, args.dpi)


if __name__ == "__main__":
    main()