import requests
import re
DGIDB_URL           = "https://dgidb.org/api/v2/interactions.json"
OPENTARGETS_URL     = "https://api.platform.opentargets.org/api/v4/graphql"
ENSEMBL_LOOKUP_URL  = "https://rest.ensembl.org/xrefs/symbol/homo_sapiens"
TIMEOUT             = 15


def _clean_symbol(gene: str) -> str:
    """Приведення назви гена до стандартного вигляду."""
    return re.sub(r"\s+", "", gene.strip()).upper()


def _symbol_to_ensembl(gene_symbol: str) -> str | None:
    """Конвертація символу гена у Ensembl ID через Ensembl REST API."""
    try:
        r = requests.get(
            f"{ENSEMBL_LOOKUP_URL}/{gene_symbol}",
            params={"content-type": "application/json"},
            timeout=TIMEOUT,
        )
        r.raise_for_status()
        data = r.json()
        for entry in data:
            if entry.get("type") == "gene":
                return entry.get("id")
        return None
    except Exception:
        return None



def fetch_dgidb_interactions(gene_symbol: str) -> list[dict]:
    """Отримує дані про ліки з DGIdb."""
    gene_symbol = _clean_symbol(gene_symbol)
    try:
        r = requests.get(DGIDB_URL, params={"genes": gene_symbol}, timeout=TIMEOUT)
        r.raise_for_status()
        data = r.json()

        interactions_found = []
        for term in data.get("matchedTerms", []):
            if term.get("searchTerm", "").upper() != gene_symbol:
                continue
            for interaction in term.get("interactions", []):
                types = interaction.get("interactionTypes", [])
                action_labels = (
                    [t.get("type", "unknown") for t in types] if types else ["unknown"]
                )
                sources = interaction.get("sources", [])
                interactions_found.append({
                    "drug": interaction.get("drugName"),
                    "actions": action_labels,
                    "score": interaction.get("score", 0),
                    "sources_count": len(sources),
                    "is_high_confidence": interaction.get("score", 0) > 5,
                })

        return sorted(interactions_found, key=lambda x: x["score"], reverse=True)

    except Exception as e:
        print(f"DGIdb Error: {e}")
        return []



def _run_ot_query(query: str, variables: dict) -> dict | None:
    """Виконує GraphQL-запит до Open Targets API."""
    try:
        r = requests.post(
            OPENTARGETS_URL,
            json={"query": query, "variables": variables},
            timeout=TIMEOUT,
        )
        r.raise_for_status()
        data = r.json()
        return data.get("data", {}).get("target")
    except Exception as e:
        print(f"Open Targets Error: {e}")
        return None


def get_ot_tractability(ensembl_id: str) -> list[dict]:
    """
    Отримує tractability assessments — придатність білка для таргетування
    по модальностях: Small Molecule, Antibody, PROTAC, Other.
    """
    query = """
    query tractability($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        tractability {
          label
          modality
          value
        }
      }
    }
    """
    result = _run_ot_query(query, {"ensemblId": ensembl_id})
    if not result:
        return []
    return result.get("tractability") or []


def get_ot_genetic_constraint(ensembl_id: str) -> list[dict]:
    """
    Отримує genetic constraint (pLI, oe ratio) — нетерпимість гена до мутацій.
    Високий pLI / низький oe → ген критичний для розвитку (важливо для ASD).
    """
    query = """
    query geneticConstraint($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        geneticConstraint {
          constraintType
          exp
          obs
          score
          oe
          oeLower
          oeUpper
        }
      }
    }
    """
    result = _run_ot_query(query, {"ensemblId": ensembl_id})
    if not result:
        return []
    return result.get("geneticConstraint") or []



def get_ot_known_drugs(ensembl_id: str, size: int = 10) -> list[dict]:
    """
    Отримує відомі ліки для мішені з ChEMBL через Open Targets.
    Включає фази клінічних випробувань та механізми дії.
    """
    query = """
    query knownDrugs($ensemblId: String!, $size: Int!) {
      target(ensemblId: $ensemblId) {
        knownDrugs(size: $size) {
          rows {
            drug {
              id
              name
              drugType
            }
            mechanismOfAction
            phase
            status
            disease {
              id
              name
            }
          }
        }
      }
    }
    """
    result = _run_ot_query(query, {"ensemblId": ensembl_id, "size": size})
    if not result:
        return []
    kd = result.get("knownDrugs") or {}
    return kd.get("rows") or []


def get_ot_safety_liabilities(ensembl_id: str) -> list[dict]:
    """
    Отримує safety liabilities — профіль безпеки гена з клінічних баз даних.
    Включає побічні ефекти, токсичність та фармаконаглядові-сигнали.
    """
    query = """
    query safetyLiabilities($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        safetyLiabilities {
          event
          eventId
          effects {
            direction
            dosing
          }
          datasource
          literature
          url
          studies {
            description
            name
            type
          }
        }
      }
    }
    """
    result = _run_ot_query(query, {"ensemblId": ensembl_id})
    if not result:
        return []
    return result.get("safetyLiabilities") or []


def get_opentargets_profile(gene_symbol: str) -> dict:
    """
    Головна функція: збирає повний Open Targets профіль гена.
    Конвертує символ → Ensembl ID → виконує всі запити.
    
    Повертає:
        ensembl_id, tractability, genetic_constraint,
        associated_diseases, known_drugs, safety_liabilities
    """
    gene_symbol = _clean_symbol(gene_symbol)

    ensembl_id = _symbol_to_ensembl(gene_symbol)
    if not ensembl_id:
        return {
            "ensembl_id": None,
            "error": f"Could not resolve Ensembl ID for {gene_symbol}",
            "tractability": [],
            "genetic_constraint": [],
            "associated_diseases": [],
            "known_drugs": [],
            "safety_liabilities": [],
        }

    return {
        "ensembl_id": ensembl_id,
        "error": None,
        "tractability":          get_ot_tractability(ensembl_id),
        "genetic_constraint":    get_ot_genetic_constraint(ensembl_id),
        "known_drugs":           get_ot_known_drugs(ensembl_id),
        "safety_liabilities":    get_ot_safety_liabilities(ensembl_id),
    }


def analyze_druggability_v2(gene_symbol: str, all_data: dict) -> dict:
    """
    Об'єднує дані з DGIdb (ліки), Open Targets (безпека + tractability) та GO-анотацій.

    Тієри:
        Gold   — є валідовані ліки у DGIdb або Open Targets knownDrugs (phase >= 3)
        Silver — druggable клас білка (кіназа / рецептор / іонний канал)
                 або tractability: Small Molecule / Antibody = True
        Bronze — недостатньо доказів
    """
    gene_symbol = _clean_symbol(gene_symbol)

    dgidb_interactions = fetch_dgidb_interactions(gene_symbol)
    ot_profile         = get_opentargets_profile(gene_symbol)

    ot_drugs           = ot_profile.get("known_drugs", [])
    tractability       = ot_profile.get("tractability", [])
    safety_liabilities = ot_profile.get("safety_liabilities", [])
    genetic_constraint = ot_profile.get("genetic_constraint", [])

    go_items = all_data.get("go", [])
    go_text  = " ".join(
        str(i.get("label", "")).lower() for i in go_items if isinstance(i, dict)
    )

    pli_score = None
    for c in genetic_constraint:
        if c.get("constraintType") == "lof":
            pli_score = c.get("score")
            break

    if pli_score is not None:
        if pli_score >= 0.9:
            safety_status = "High Constraint (pLI≥0.9) — caution"
        elif pli_score >= 0.5:
            safety_status = "Moderate Constraint"
        else:
            safety_status = "Low Constraint — Safe Target"
    else:
        safety_status = "No constraint data"

    tractable_modalities = [
        t["modality"] for t in tractability if t.get("value") is True
    ]

    ot_top_drugs = [
        {
            "drug": row.get("drug", {}).get("name"),
            "phase": row.get("phase"),
            "mechanism": row.get("mechanismOfAction"),
            "disease": row.get("disease", {}).get("name"),
            "source": "OpenTargets",
        }
        for row in sorted(ot_drugs, key=lambda x: x.get("phase") or 0, reverse=True)[:5]
    ]

    dgidb_top_drugs = [
        {**d, "source": "DGIdb"}
        for d in dgidb_interactions[:5]
    ]

    all_top_drugs = ot_top_drugs + dgidb_top_drugs

    report = {
        "gene":              gene_symbol,
        "tier":              "Bronze",
        "score":             0.1,
        "safety_status":     safety_status,
        "pli_score":         pli_score,
        "tractable_modalities": tractable_modalities,
        "safety_liabilities_count": len(safety_liabilities),
        "verdict":           "Low therapeutic evidence.",
        "top_drugs":         all_top_drugs,
        "strategy":          "Unknown",
        "ot_ensembl_id":     ot_profile.get("ensembl_id"),
    }

    has_phase3_drugs = any((row.get("phase") or 0) >= 3 for row in ot_drugs)
    has_any_drugs    = bool(dgidb_interactions) or bool(ot_drugs)

    if has_phase3_drugs:
        report["tier"]     = "Gold"
        report["score"]    = 1.0
        report["verdict"]  = f"Validated target. Phase 3+ drugs in Open Targets."
        report["strategy"] = tractable_modalities[0] if tractable_modalities else "Known Drug"
        if pli_score and pli_score >= 0.9:
            report["verdict"] += " ⚠️ WARNING: High genetic constraint (pLI≥0.9)."

    elif has_any_drugs:
        report["tier"]     = "Gold"
        report["score"]    = 0.9
        report["verdict"]  = (
            f"Validated target. "
            f"{len(dgidb_interactions)} DGIdb + {len(ot_drugs)} OT drugs found."
        )
        report["strategy"] = tractable_modalities[0] if tractable_modalities else "Drug Interaction"

    elif tractable_modalities:
        report["tier"]     = "Silver"
        report["score"]    = 0.8
        report["verdict"]  = f"Tractable target: {', '.join(tractable_modalities)}."
        report["strategy"] = tractable_modalities[0]

    elif any(kw in go_text for kw in ["kinase", "receptor", "ion channel"]):
        report["tier"]     = "Silver"
        report["score"]    = 0.7
        report["verdict"]  = "High potential: druggable protein class (GO)."
        report["strategy"] = "Small Molecule / Antibody"

    return report