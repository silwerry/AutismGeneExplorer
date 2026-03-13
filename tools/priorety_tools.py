def calculate_go_score(go_terms):
    if not isinstance(go_terms, list) or len(go_terms) == 0:
        return 0.0

    relevant_keywords = [
        "synapse", "synaptic", "postsynaptic", "presynaptic",
        "neuron", "neuro", "axon", "dendrite", "brain", 
        "behavior", "cognition", "learning", "neurogenesis", "asd", "autism"
    ]

    relevant_count = 0
    for term in go_terms:
        label = ""
        if isinstance(term, dict):
            label = term.get("label", "").lower()
        elif isinstance(term, str):
            label = term.lower()
            
        if any(keyword in label for keyword in relevant_keywords):
            relevant_count += 1

    ratio = relevant_count / len(go_terms)
    return round(min(0.2, ratio * 0.2), 4)

def discover_connections(gene_symbol, all_data, df_sfari):
    insights = []
    trrust_list = all_data.get("trrust") or []
    sfari_genes_list = df_sfari['gene_symbol'].tolist() if df_sfari is not None else []

    for item in trrust_list:
        tf = item.get("TF")
        if tf in sfari_genes_list:
            insights.append(f"Regulatory Link: {gene_symbol} is regulated by {tf} (SFARI Risk Gene).")

    ensembl_data = all_data.get("ensembl") or []
    
    orthologs_count = 0
    paralogs_count = 0
    
    for entry in ensembl_data:
        h_type = entry.get("type", "").lower()
        if "ortholog" in h_type:
            orthologs_count += 1
        elif "paralog" in h_type:
            paralogs_count += 1

    if orthologs_count > 0:
        insights.append(f"Evolutionary Insight: Gene is highly conserved with {orthologs_count} orthologs found.")
    
    if paralogs_count > 0:
        insights.append(f"Duplication Alert: Found {paralogs_count} paralogs. This might indicate functional redundancy.")
    elif orthologs_count > 0 and paralogs_count == 0:
        insights.append(f"High Vulnerability: No paralogs found (single-copy gene). Mutations might have higher impact.")

    return insights

def calculate_consensus_score(gene_data):
    score = 0.0

    sfari = gene_data.get("sfari") or {}
    raw_cat = sfari.get("category", "")
    try:
        cat_int = int(float(str(raw_cat).strip()))
    except (ValueError, TypeError):
        cat_int = None

    if cat_int == 1:
        score += 0.6
    elif cat_int == 2:
        score += 0.4
    elif cat_int == 3:
        score += 0.2

    trrust = gene_data.get("trrust") or []
    if isinstance(trrust, list):
        score += min(0.2, len(trrust) * 0.05)

    go_terms = gene_data.get("go") or []
    score += calculate_go_score(go_terms)

    ensembl_data = gene_data.get("ensembl") or []
    if len(ensembl_data) > 100:
        score += 0.05

    return round(min(score, 1.0), 3)