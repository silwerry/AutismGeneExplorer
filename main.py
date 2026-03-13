from mcp.server.fastmcp import FastMCP
from tools.trrust_tools import load_trrust, get_regulators_df, get_targets_df
from tools.input_tools import load_table, list_input_files as list_input
from tools.integration_tools import genes_in_trrust, guess_gene_column
from tools.sfari_tools import load_sfari, get_sfari_genes_df
from tools.go_tools import (
    map_gene_to_uniprot, fetch_go_biological_process, fetch_go_function,
    format_go_annotations, fetch_go_cellular_component
)
from tools.ensembl_tools import (
    fetch_orthologs, format_orthologs, fetch_gene_tree,
    fetch_gene_coords, fetch_regulatory_features, find_regulatory_hits, save_tree_visualization
)
from tools.priorety_tools import calculate_consensus_score, discover_connections
from tools.druggability_analysis import (
    analyze_druggability_v2,
    get_opentargets_profile,
    get_ot_tractability,
    get_ot_genetic_constraint,
    get_ot_known_drugs,
    get_ot_safety_liabilities,
    _symbol_to_ensembl,
)


df_trrust = load_trrust()
df_sfari  = load_sfari()
mcp = FastMCP("asd-isef")


# ── Input / Table tools ───────────────────────────────────────────────────────

@mcp.tool()
def list_input_files() -> str:
    files = list_input()
    return "No input files found." if not files else "\n".join(files)


@mcp.tool()
def describe_table(filename: str) -> str:
    df = load_table(filename)
    return "\n".join(f"{c}: {t}" for c, t in df.dtypes.items())


@mcp.tool()
def suggest_gene_column(filename: str) -> str:
    df = load_table(filename)
    ranked = guess_gene_column(df)
    if not ranked:
        return "No candidate gene columns detected."
    return "\n".join(
        f"{r['column']} | gene_like={r['gene_like_fraction']} | n={r['n_unique']}"
        for r in ranked
    )


@mcp.tool()
def data_in_trrust(filename: str, gene_column: str) -> str:
    df_input = load_table(filename)
    hits = genes_in_trrust(df_input, df_trrust, gene_column)
    return "No genes found in TRRUST." if not hits else "\n".join(hits)


# ── TRRUST tools ──────────────────────────────────────────────────────────────

@mcp.tool()
def trrust_regulators(genes: list[str]) -> str:
    hits = get_regulators_df(df_trrust, genes)
    if hits.empty:
        return f"No regulators found for {', '.join(genes)}."
    lines = []
    for _, r in hits.iterrows():
        lines.append(f"{r.TF} | {r.Target} | {r.Mode} | {r.Source}")
    return "\n".join(lines)


@mcp.tool()
def trrust_targets(tfs: list[str]) -> str:
    hits = get_targets_df(df_trrust, tfs)
    if hits.empty:
        return f"No targets found for {', '.join(tfs)}."
    lines = []
    for _, r in hits.iterrows():
        lines.append(f"{r.TF} | {r.Target} | {r.Mode} | {r.Source}")
    return "\n".join(lines)


# ── SFARI tools ───────────────────────────────────────────────────────────────

@mcp.tool()
def sfari_genes(genes: list[str]) -> str:
    hits = get_sfari_genes_df(df_sfari, genes)
    if hits.empty:
        return "No genes found in SFARI."
    lines = []
    for _, row in hits.iterrows():
        lines.append(
            f"{row.gene_symbol} | "
            f"score={row.gene_score} | "
            f"category={row.genetic_category}"
        )
    return "\n".join(lines)


# ── GO tools ──────────────────────────────────────────────────────────────────

@mcp.tool()
def go_molecular_function(gene_symbol: str) -> str:
    """Map human gene symbol → UniProtKB → extract GO molecular functions."""
    uniprot_id = map_gene_to_uniprot(gene_symbol)
    if not uniprot_id:
        return f"No human UniProtKB ID found for {gene_symbol}."
    annotations = fetch_go_function(uniprot_id)
    formatted = format_go_annotations(annotations)
    return f"{gene_symbol} → {uniprot_id}\n{formatted}"


@mcp.tool()
def biological_process(gene_symbol: str) -> str:
    """Map human gene symbol → UniProtKB → extract biological functions."""
    uniprot_id = map_gene_to_uniprot(gene_symbol)
    if not uniprot_id:
        return f"No human UniProtKB ID found for {gene_symbol}."
    annotations = fetch_go_biological_process(uniprot_id)
    formatted = format_go_annotations(annotations)
    return f"{gene_symbol} → {uniprot_id}\n{formatted}"


@mcp.tool()
def cellular_component(gene_symbol: str) -> str:
    """Map human gene symbol → UniProtKB → extract cellular component."""
    uniprot_id = map_gene_to_uniprot(gene_symbol)
    if not uniprot_id:
        return f"No human UniProtKB ID found for {gene_symbol}."
    annotations = fetch_go_cellular_component(uniprot_id)
    formatted = format_go_annotations(annotations)
    return f"{gene_symbol} → {uniprot_id}\n{formatted}"


# ── Ensembl tools ─────────────────────────────────────────────────────────────

@mcp.tool()
def gene_orthologs(gene_symbol: str, species: str = "homo_sapiens") -> str:
    """Fetch orthologs for a gene using Ensembl Comparative Genomics API."""
    orthologs = fetch_orthologs(species, gene_symbol)
    if not orthologs:
        return f"No orthologs found for {gene_symbol} in {species}."
    formatted = format_orthologs(orthologs)
    return f"Orthologs for {gene_symbol} ({species}):\n{formatted}"


@mcp.tool()
def gene_tree(gene_symbol: str, species: str = "homo_sapiens") -> str:
    """Fetch gene tree for a gene using Ensembl REST API and generate visualization."""
    tree = fetch_gene_tree(species, gene_symbol)
    if not tree:
        return f"No gene tree found for {gene_symbol} in {species}."

    try:
        image_path = save_tree_visualization(gene_symbol, species)
        viz_status = f"\n✅ Візуалізацію збережено: {image_path}"
    except Exception as e:
        viz_status = f"\n❌ Помилка візуалізації: {str(e)}"

    tree_data = tree.get("tree", {}) if isinstance(tree, dict) else {}
    tree_id = (
        tree.get("id") or tree.get("stable_id") or tree.get("tree_id")
        or tree_data.get("id") or tree_data.get("stable_id")
        or tree_data.get("tree_id") or "Unknown tree ID"
    )

    def parse_node(node: dict, depth: int = 0) -> list[str]:
        if not isinstance(node, dict):
            return []
        lines = []
        indent = "  " * depth
        node_type = node.get("node_type", "")
        if isinstance(node_type, dict):
            node_type = node_type.get("value", "")
        tax = node.get("taxonomy", {})
        tax_name = tax.get("scientific_name", "") if isinstance(tax, dict) else str(tax or "")
        events = node.get("events", {})
        duplication = ""
        if isinstance(events, dict):
            ev_type = events.get("type", {})
            duplication = ev_type.get("value", "") if isinstance(ev_type, dict) else str(ev_type or "")
        children = node.get("children", [])
        if not children:
            genes = node.get("genes", [])
            for gene in genes:
                if not isinstance(gene, dict): continue
                gene_info = gene.get("gene", {})
                if not isinstance(gene_info, dict): continue
                gene_id   = gene_info.get("id", "?")
                gene_name = gene_info.get("name", "")
                sp_info   = gene_info.get("species", {})
                sp = sp_info.get("scientific_name", tax_name) if isinstance(sp_info, dict) else tax_name
                lines.append(f"{indent}🍃 {gene_name or gene_id} [{sp}]")
        else:
            label     = tax_name or str(node_type) or "node"
            event_str = f" ({duplication})" if duplication else ""
            lines.append(f"{indent}🌿 {label}{event_str}")
            for child in children:
                lines.extend(parse_node(child, depth + 1))
        return lines

    tree_lines = parse_node(tree_data)
    MAX_LINES  = 200
    truncated  = ""
    if len(tree_lines) > MAX_LINES:
        truncated  = f"\n... (показано {MAX_LINES} з {len(tree_lines)} вузлів)"
        tree_lines = tree_lines[:MAX_LINES]

    tree_str = "\n".join(tree_lines) if tree_lines else "Tree structure is empty."
    return (
        f"Gene tree for {gene_symbol} ({species})\n"
        f"Tree ID: {tree_id}\n"
        f"{'='*50}\n"
        f"{tree_str}"
        f"{truncated}\n"
        f"{'='*50}"
        f"{viz_status}"
    )


@mcp.tool()
def gene_regulation(gene_symbol: str, species: str = "homo_sapiens") -> str:
    coords = fetch_gene_coords(species, gene_symbol)
    if not coords:
        return f"Ген {gene_symbol} не знайдено."
    chrom       = coords.get("seq_region_name")
    start, end  = coords.get("start"), coords.get("end")
    regs        = fetch_regulatory_features(species, chrom, start, end)
    if not regs:
        return "Регуляторних елементів не знайдено."
    lines = [f"Результати для гена {gene_symbol} (chr{chrom}:{start}-{end}):\n"]
    lines.append(f"{'ID':<20} | {'Тип':<20} | {'Координати'}")
    lines.append("-" * 70)
    for r in regs:
        f_id   = r.get("id", "N/A")
        f_type = r.get("feature_type", "Unknown")
        loc    = f"{r.get('start')}-{r.get('end')}"
        lines.append(f"{f_id:<20} | {f_type:<20} | {loc}")
    return "\n".join(lines)


# ── Priority / Report tools ───────────────────────────────────────────────────

@mcp.tool()
def calculate_gene_priority(gene_symbol: str) -> str:
    """
    Швидкий розрахунок пріоритетності гена на основі SFARI, TRRUST, GO та Ensembl.
    """
    sfari_hits  = get_sfari_genes_df(df_sfari, [gene_symbol])
    sfari_data  = {"category": sfari_hits.iloc[0].gene_score} if not sfari_hits.empty else None
    trrust_hits = get_regulators_df(df_trrust, [gene_symbol])
    trrust_list = trrust_hits.to_dict("records") if not trrust_hits.empty else []
    uniprot_id  = map_gene_to_uniprot(gene_symbol)
    go_data     = fetch_go_biological_process(uniprot_id) if uniprot_id else []
    homology    = fetch_orthologs("homo_sapiens", gene_symbol)

    all_data = {
        "sfari":   sfari_data,
        "trrust":  trrust_list,
        "go":      go_data,
        "ensembl": homology,
    }

    consensus_score = calculate_consensus_score(all_data)
    connections     = discover_connections(gene_symbol, all_data, df_sfari)
    percentage      = consensus_score * 100
    level = "High 🔴" if consensus_score >= 0.7 else "Medium 🟡" if consensus_score >= 0.4 else "Low 🟢"

    report = [
        f"## Gene Priority Score: {gene_symbol}",
        f"**Consensus Score:** {consensus_score} ({percentage:.1f}%)",
        f"**Risk Level:** {level}",
        "---",
        "### Score Components:",
        f"* **SFARI Category:** {sfari_data['category'] if sfari_data else 'Not Listed'}",
        f"* **Regulatory Links:** {len(trrust_list)} TFs",
        f"* **GO Terms Found:** {len(go_data)}",
        f"* **Homology Depth:** {len(homology)} orthologs/paralogs",
        "",
    ]
    if connections:
        report.append("### Key Insights:")
        for conn in connections:
            report.append(f"* {conn}")

    return "\n".join(report)


@mcp.tool()
def get_gene_report(gene_symbol: str) -> dict:
    """
    Повний звіт по гену: пріоритет, регуляція та фармакологічний профіль.
    Використовує Open Targets замість DepMap для оцінки безпеки.
    """
    sfari_hits  = get_sfari_genes_df(df_sfari, [gene_symbol])
    sfari_data  = {"category": sfari_hits.iloc[0].gene_score} if not sfari_hits.empty else None
    trrust_hits = get_regulators_df(df_trrust, [gene_symbol])
    trrust_list = trrust_hits.to_dict("records") if not trrust_hits.empty else []
    uniprot_id  = map_gene_to_uniprot(gene_symbol)
    go_bp       = fetch_go_biological_process(uniprot_id) if uniprot_id else []
    go_mf       = fetch_go_function(uniprot_id) if uniprot_id else []
    go_cc       = fetch_go_cellular_component(uniprot_id) if uniprot_id else []
    all_go      = go_bp + go_mf + go_cc
    homology    = fetch_orthologs("homo_sapiens", gene_symbol)
    coords      = fetch_gene_coords("homo_sapiens", gene_symbol)
    reg_data    = []
    if coords:
        reg_data = fetch_regulatory_features(
            "homo_sapiens", coords["seq_region_name"], coords["start"], coords["end"]
        )

    all_data = {
        "sfari":      sfari_data,
        "trrust":     trrust_list,
        "go":         all_go,
        "ensembl":    homology,
        "regulation": reg_data,
    }

    score        = calculate_consensus_score(all_data)
    insights     = discover_connections(gene_symbol, all_data, df_sfari)
    drug_analysis = analyze_druggability_v2(gene_symbol, all_data)

    return {
        "gene":           gene_symbol,
        "priority_score": score,
        "pharmacogenomics_profile": {
            "tier":                      drug_analysis.get("tier"),
            "strategy":                  drug_analysis.get("strategy"),
            "safety_status":             drug_analysis.get("safety_status"),
            "pli_score":                 drug_analysis.get("pli_score"),
            "tractable_modalities":      drug_analysis.get("tractable_modalities"),
            "safety_liabilities_count":  drug_analysis.get("safety_liabilities_count"),
            "validated_drugs":           drug_analysis.get("top_drugs", []),
        },
        "knowledge_graph": insights,
        "regulation_summary": {
            "total_elements": len(reg_data),
            "promoters":  len([r for r in reg_data if "Promoter" in r.get("feature_type", "")]),
            "enhancers":  len([r for r in reg_data if "Enhancer" in r.get("feature_type", "")]),
        },
        "sfari_status": sfari_data.get("category") if sfari_data else "Not in SFARI",
    }


# ── Open Targets tools ────────────────────────────────────────────────────────

@mcp.tool()
def get_ot_tractability(gene_symbol: str) -> str:
    """
    Оцінює придатність гена для таргетування (Open Targets tractability).
    Показує доступні модальності: Small Molecule, Antibody, PROTAC.
    """
    ensembl_id = _symbol_to_ensembl(gene_symbol)
    if not ensembl_id:
        return f"Could not resolve Ensembl ID for {gene_symbol}."

    from tools.druggability_analysis import get_ot_tractability as _fetch
    items = _fetch(ensembl_id)
    if not items:
        return f"No tractability data found for {gene_symbol} ({ensembl_id})."

    lines = [f"Tractability for {gene_symbol} ({ensembl_id}):"]
    for item in items:
        status = "✅" if item.get("value") else "❌"
        lines.append(f"  {status} [{item.get('modality')}] {item.get('label')}")
    return "\n".join(lines)


@mcp.tool()
def get_genetic_constraint(gene_symbol: str) -> str:
    """
    Отримує genetic constraint (pLI, oe ratio) з Open Targets.
    Високий pLI → ген нетерпимий до мутацій, критичний для розвитку.
    """
    ensembl_id = _symbol_to_ensembl(gene_symbol)
    if not ensembl_id:
        return f"Could not resolve Ensembl ID for {gene_symbol}."

    from tools.druggability_analysis import get_ot_genetic_constraint as _fetch
    items = _fetch(ensembl_id)
    if not items:
        return f"No genetic constraint data found for {gene_symbol} ({ensembl_id})."

    lines = [f"Genetic Constraint for {gene_symbol} ({ensembl_id}):"]
    for c in items:
        c_type = c.get("constraintType", "?")
        score  = c.get("score")
        oe     = c.get("oe")
        obs    = c.get("obs")
        exp    = c.get("exp")
        risk   = ""
        if c_type == "lof" and score is not None:
            if score >= 0.9:
                risk = " 🔴 High constraint (pLI≥0.9)"
            elif score >= 0.5:
                risk = " 🟡 Moderate constraint"
            else:
                risk = " 🟢 Low constraint"
        lines.append(
            f"  [{c_type}] pLI/score={score} | oe={oe} "
            f"(obs={obs}, exp={exp}){risk}"
        )
    return "\n".join(lines)


@mcp.tool()
def get_known_drugs(gene_symbol: str) -> str:
    """
    Повертає відомі ліки для гена-мішені з Open Targets / ChEMBL.
    Включає фазу клінічних випробувань, механізм дії та показання.
    """
    ensembl_id = _symbol_to_ensembl(gene_symbol)
    if not ensembl_id:
        return f"Could not resolve Ensembl ID for {gene_symbol}."

    from tools.druggability_analysis import get_ot_known_drugs as _fetch
    rows = _fetch(ensembl_id)
    if not rows:
        return f"No known drugs found for {gene_symbol} ({ensembl_id})."

    lines = [f"Known Drugs for {gene_symbol} ({ensembl_id}):"]
    for row in rows:
        drug     = row.get("drug", {})
        name     = drug.get("name", "Unknown")
        d_type   = drug.get("drugType", "")
        phase    = row.get("phase", "?")
        status   = row.get("status", "")
        moa      = row.get("mechanismOfAction", "")
        disease  = row.get("disease", {}).get("name", "")
        lines.append(
            f"     {name} ({d_type}) | Phase {phase} | {status}\n"
            f"     MOA: {moa}\n"
            f"     Indication: {disease}"
        )
    return "\n".join(lines)


@mcp.tool()
def get_safety_liabilities(gene_symbol: str) -> str:
    """
    Отримує safety liabilities гена з Open Targets.
    Показує відомі побічні ефекти, токсичні події та фармаковіgilance-сигнали.
    """
    ensembl_id = _symbol_to_ensembl(gene_symbol)
    if not ensembl_id:
        return f"❌ Could not resolve Ensembl ID for {gene_symbol}."

    from tools.druggability_analysis import get_ot_safety_liabilities as _fetch
    items = _fetch(ensembl_id)
    if not items:
        return f"No safety liabilities found for {gene_symbol} ({ensembl_id}). Gene appears safe to target."

    lines = [f"Safety Liabilities for {gene_symbol} ({ensembl_id}):"]
    for item in items:
        event      = item.get("event", "Unknown event")
        datasource = item.get("datasource", "")
        effects    = item.get("effects", [])
        effect_str = ", ".join(
            f"{e.get('direction')} ({e.get('dosing')})" for e in effects if e
        ) if effects else "N/A"
        studies    = item.get("studies", [])
        study_str  = "; ".join(
            s.get("description", "") for s in studies if s.get("description")
        ) if studies else ""
        lines.append(f"  ⚠️  {event} [{datasource}]")
        lines.append(f"      Effects: {effect_str}")
        if study_str:
            lines.append(f"      Studies: {study_str[:200]}")
    return "\n".join(lines)


@mcp.tool()
def get_therapeutic_insight(gene_symbol: str) -> dict:
    """
    терапевтичний профіль гена:
    DGIdb (ліки) + Open Targets (tractability, safety, known drugs) + стратегія.
    """
    uniprot_id = map_gene_to_uniprot(gene_symbol)
    go_data    = []
    if uniprot_id:
        mf      = fetch_go_function(uniprot_id)
        cc      = fetch_go_cellular_component(uniprot_id)
        go_data = mf + cc

    result = analyze_druggability_v2(gene_symbol, {"go": go_data})

    return {
        "gene":       gene_symbol,
        "assessment": result,
    }


if __name__ == "__main__":
    mcp.run(transport="stdio")