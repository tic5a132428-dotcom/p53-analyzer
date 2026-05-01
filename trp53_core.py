import requests
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio.Align import PairwiseAligner

# --- BASE DE DADOS ---
VARIANT_DB = {
    "R175H": {
        "cancers": ["Mama", "Pulmão", "Ovário"],
        "phenotypes": ["Neoplasia", "Proliferação celular"],
        "expression": ["Mama", "Pulmão", "Ovário"],
        "gnomad_freq": 1e-5,
        "mouse_model": "Desenvolve tumores espontâneos",
        "clinical_trials": "Ensaios com terapias alvo p53",
        "article_desc": "o impacto funcional da mutação R175H na proteína p53 e a sua reativação farmacológica"
    }
}

TREATMENTS_DB = {
    "R175H": ["APR-246 (Eprenetapopt)", "Terapias alvo MDM2"],
    "R172H": ["Modelos pré-clínicos com APR-246"]
}

ARTICLES_DB = {
    "R175H": "https://pubmed.ncbi.nlm.nih.gov/12826609/"
}

CHROMOSOME_MAP = {
    "Humano": "17p13.1",
    "Rato": "11B2"
}

DOMINIOS = [
    {"nome": "TAD", "inicio": 1, "fim": 42, "cor": "#FF9999"},
    {"nome": "PRD", "inicio": 63, "fim": 94, "cor": "#99FF99"},
    {"nome": "DBD", "inicio": 95, "fim": 289, "cor": "#9999FF"},
    {"nome": "OD", "inicio": 323, "fim": 355, "cor": "#FFFF99"},
    {"nome": "CTD", "inicio": 356, "fim": 393, "cor": "#FFCC99"},
]

def get_seq(uid):
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=5)
        if r.status_code == 200:
            return "".join(r.text.split("\n")[1:])
    except:
        pass
    return None

mouse_seq = get_seq("P02340")
human_seq = get_seq("P04637")

aligner = PairwiseAligner(mode="global")
alignment = aligner.align(mouse_seq, human_seq)[0] if mouse_seq and human_seq else None

def map_position(pos, from_mouse=True):
    if alignment is None: return pos
    m_idx, h_idx = (0, 1) if from_mouse else (1, 0)
    blocks = alignment.aligned
    c1 = c2 = 0
    for b1, b2 in zip(blocks[m_idx], blocks[h_idx]):
        for _ in range(b1[1] - b1[0]):
            c1 += 1
            c2 += 1
            if c1 == pos: return c2
    return pos

def get_domain(pos):
    for d in DOMINIOS:
        if d["inicio"] <= pos <= d["fim"]: return d["nome"]
    return "Linker"

def get_clinvar(ref, pos, alt):
    try:
        q = f"TP53 {ref}{pos}{alt}"
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                         params={"db": "clinvar", "term": q, "retmode": "json"}, timeout=3).json()
        ids = r.get("esearchresult", {}).get("idlist", [])
        if ids:
            s = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                             params={"db": "clinvar", "id": ids[0], "retmode": "json"}, timeout=3).json()
            return s["result"][ids[0]]["clinical_significance"]["description"]
    except:
        pass
    return None

def compute_score(mutation, pos, clinvar, has_equivalent):
    score = 0
    details = []
    if 95 <= pos <= 289:
        score += 2
        details.append(("Domínio DBD", "+2", "Região crítica (máx: +2)"))
    if mutation in VARIANT_DB:
        score += 3
        details.append(("Mutação conhecida", "+3", "Evidência (máx: +3)"))
    if mutation in ["R175H", "R172H"]:
        score += 3
        details.append(("Hotspot", "+3", "Alta frequência (máx: +3)"))
    if has_equivalent:
        score += 2
        details.append(("Equivalente interespécies", "+2", "Impacto conservado (máx: +2)"))
    if clinvar and "pathogenic" in clinvar.lower():
        score += 3
        details.append(("ClinVar", "+3", "Clínico (máx: +3)"))
    return min(score, 11), details

def classify(score):
    if score >= 6: return "Maligno"
    if score >= 3: return "Potencialmente maligno"
    if score >= 1: return "Indeterminado"
    return "Benigno"

def plot_protein_map(pos, mutation_name):
    fig, ax = plt.subplots(figsize=(10, 2.5))
    ax.add_patch(patches.Rectangle((1, 0.4), 393, 0.2, color="lightgrey", zorder=1))
    for d in DOMINIOS:
        ax.add_patch(patches.Rectangle((d["inicio"], 0.4), d["fim"] - d["inicio"], 0.2, color=d["cor"], zorder=2))
        ax.text((d["inicio"] + d["fim"]) / 2, 0.5, d["nome"], ha='center', va='center', fontsize=8, fontweight='bold')
    ax.plot(pos, 0.5, marker='v', color='red', markersize=10, zorder=3)
    ax.annotate(mutation_name, xy=(pos, 0.6), xytext=(pos, 0.9),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1), ha='center', color='red')
    ax.set_xlim(0, 400)
    ax.set_ylim(0, 1.2)
    ax.set_yticks([])
    ax.set_xlabel("Aminoácidos")
    for s in ["top", "left", "right"]: ax.spines[s].set_visible(False)
    return fig

def generate_report(mutation, pos, clinvar, score, details, species, equiv, human_pos):
    db = VARIANT_DB.get(mutation, {})
    treatments = TREATMENTS_DB.get(mutation, ["Não disponível"])
    article = ARTICLES_DB.get(mutation, "Não disponível")
    art_desc = db.get("article_desc", "esta mutação específica")

    h_chrom = CHROMOSOME_MAP["Humano"]
    m_chrom = CHROMOSOME_MAP["Rato"]

    table = "\n".join([f"{d[0]} | {d[1]} | {d[2]}" for d in details])

    return f"""RELATÓRIO DE ANÁLISE - p53

Mutação Analisada: {mutation}
Posição Original: {pos}
Domínio: {get_domain(pos)}
Espécie de Origem: {species} ({CHROMOSOME_MAP[species]})

--- MAPEAMENTO INTERESPÉCIES ---
Equivalente em Rato: {equiv if species == "Humano" else mutation} | Posição: {pos if species == "Rato" else map_position(pos, False)} | Cromossoma: {m_chrom}
Equivalente em Humano: {mutation if species == "Humano" else equiv} | Posição: {human_pos} | Cromossoma: {h_chrom}

Cancros associados: {", ".join(db.get("cancers", ["N/A"]))}
Fenótipos: {", ".join(db.get("phenotypes", ["N/A"]))}
Frequência (gnomAD): {db.get("gnomad_freq", "N/A")}

Modelo de rato: {db.get("mouse_model", "N/A")}
Ensaios clínicos: {db.get("clinical_trials", "N/A")}
Tratamentos potenciais: {", ".join(treatments)}
Dados clínicos: {clinvar if clinvar else "Não disponível"}

--- SCORE ---
{table}

Score final: {score}/11
Classificação: {classify(score)}

--- FONTES ---
• ClinVar: Fornece a significância clínica e patogenicidade.
• gnomAD: Fornece a frequência da mutação em populações saudáveis.
• COSMIC: Base de dados de mutações somáticas em cancro humano.
• UniProt: Fornece a sequência proteica e anotação de domínios.
• MGI: Informação sobre fenótipos e modelos de rato.
• Alliance Genome: Integração de dados genómicos e ortologia interespécies.

--- ARTIGO ---
Se quiser mais informação, consulte: {article}
Este artigo aborda detalhadamente {art_desc}."""