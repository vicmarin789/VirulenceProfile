import streamlit as st
import pandas as pd
import requests

# ---------------------------
# Atualiza base automática da UniProt (com múltiplos nomes por gene)
# ---------------------------
def atualizar_base_virulencia():
    organismos = [
        ("Salmonella enterica", 99287, "neg"),
        ("Escherichia coli", 562, "neg"),
        ("Staphylococcus aureus", 1280, "pos"),
        ("Listeria monocytogenes", 1639, "pos"),
        ("Pseudomonas aeruginosa", 287, "neg"),
        ("Shigella flexneri", 623, "neg"),
        ("Klebsiella pneumoniae", 573, "neg"),
        ("Vibrio cholerae", 666, "neg"),
        ("Campylobacter jejuni", 197, "neg")
    ]

    genes = []
    for _, org_id, gram in organismos:
        query = f"(organism_id:{org_id}) AND (keyword:KW-0800)"
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=tsv&fields=accession,gene_names"
        r = requests.get(url)
        if r.status_code != 200:
            continue
        lines = r.text.strip().split("\n")[1:]
        for line in lines:
            parts = line.split("\t")
            if len(parts) < 2 or not parts[1]:
                continue
            gene_names = parts[1].strip().lower().split(" ")
            for name in gene_names:
                name = name.strip()
                if name:
                    genes.append({
                        "gene": name,
                        "categoria": "desconhecida",
                        "gram": gram,
                        "peso": 3,
                        "evidencia": True
                    })

    df = pd.DataFrame(genes).drop_duplicates(subset="gene")
    df.to_csv("base_virulencia.csv", index=False)

# ---------------------------
# Atualiza base completa (automática + manual) — manual tem prioridade
# ---------------------------
def atualizar_base_virulencia_completa():
    atualizar_base_virulencia()
    base_auto = pd.read_csv("base_virulencia.csv")

    try:
        base_manual = pd.read_csv("base_manual.csv")
    except FileNotFoundError:
        st.warning("⚠️ Arquivo 'base_manual.csv' não encontrado. Apenas a base automática será usada.")
        base_auto.to_csv("base_virulencia.csv", index=False)
        return len(base_auto)

    for df in [base_auto, base_manual]:
        df["gene"] = df["gene"].astype(str).str.strip().str.lower()
        df["gram"] = df["gram"].astype(str).str.strip().str.lower()
        df["categoria"] = df["categoria"].fillna("desconhecida")
        df["peso"] = pd.to_numeric(df["peso"], errors="coerce").fillna(0)
        df["evidencia"] = df["evidencia"].astype(str).str.strip().str.lower() == "true"

    # 🔹 Prioriza a base manual
    base_completa = pd.concat([base_manual, base_auto]).drop_duplicates(subset="gene", keep="first")
    base_completa.to_csv("base_virulencia.csv", index=False)
    return len(base_completa)

# ---------------------------
# Carrega base completa
# ---------------------------
@st.cache_data
def carregar_base_completa():
    base = pd.read_csv("base_virulencia.csv")
    base["gene"] = base["gene"].astype(str).str.strip().str.lower()
    base["gram"] = base["gram"].astype(str).str.strip().str.lower()
    base["categoria"] = base["categoria"].fillna("desconhecida")
    base["peso"] = pd.to_numeric(base["peso"], errors="coerce").fillna(0)
    base["evidencia"] = base["evidencia"].astype(str).str.strip().str.lower() == "true"
    return base

# ---------------------------
# Classificação refinada — preserva nome original
# ---------------------------
def classificar(df, gram="neg", base=None):
    pontuacao = 0
    detalhes = []

    df["gene_original"] = df["gene"]  # guarda para exibir
    df["gene"] = df["gene"].astype(str).str.strip().str.lower()
    df["identidade"] = pd.to_numeric(df["identidade"], errors="coerce").fillna(0)
    df["cobertura"] = pd.to_numeric(df["cobertura"], errors="coerce").fillna(0)

    for _, fator in df.iterrows():
        gene = fator["gene"]
        identidade = fator["identidade"]
        cobertura = fator["cobertura"]

        info = base[
            (base["gene"] == gene) &
            ((base["gram"] == gram) | (base["gram"] == "ambos"))
        ]

        if not info.empty:
            linha = info.iloc[0]
            if identidade >= 85 and cobertura >= 95:
                peso = linha["peso"]
                if linha["evidencia"]:
                    peso *= 1.5
                pontuacao += peso
                detalhes.append((fator["gene_original"], linha["categoria"], round(peso, 2)))
            else:
                detalhes.append((fator["gene_original"], linha["categoria"], 0))
        else:
            detalhes.append((fator["gene_original"], "não classificado", 0))

    if pontuacao >= 20:
        resultado = "Alta probabilidade de patogenicidade"
    elif pontuacao >= 10:
        resultado = "Potencial patogênico moderado"
    else:
        resultado = "Baixa probabilidade de ser patogênico"

    return pontuacao, resultado, detalhes

# ---------------------------
# Interface Streamlit
# ---------------------------
st.set_page_config(page_title="Analisador de Patogenicidade Bacteriana", page_icon="🧬")
st.title("🧬 Analisador de Patogenicidade Bacteriana")
st.write("Faça upload de um arquivo CSV com os genes identificados e selecione o tipo de bactéria.")

if st.button("🔄 Atualizar base de virulência"):
    total = atualizar_base_virulencia_completa()
    st.success(f"Base completa atualizada com {total} genes virulentos.")

modelo_csv = """gene,identidade,cobertura
inva,98.5,100
prgh,96.0,99
sipb,97.0,98
ssav,95.0,97
stn,94.0,96
ctx,97.0,99
hlya,96.0,97
sope2,95.0,96
phop,96.0,98
"""
st.download_button("📥 Baixar modelo de entrada (.csv)", modelo_csv, file_name="exemplo_entrada.csv")

base_virulencia = carregar_base_completa()

uploaded_file = st.file_uploader("Carregar arquivo CSV", type=["csv"])
gram = st.radio("Tipo de bactéria:", ("Gram-negativa", "Gram-positiva"))
gram = gram.strip().lower()
gram = "neg" if gram == "gram-negativa" else "pos"

if uploaded_file:
    df = pd.read_csv(uploaded_file)
    df.columns = df.columns.str.strip().str.lower()

    mapeamento = {}
    for col in {"gene", "identidade", "cobertura"}:
        for existente in df.columns:
            if col in existente:
                mapeamento[existente] = col
    df = df.rename(columns=mapeamento)

    st.subheader("Pré-visualização dos dados")
    st.dataframe(df)

    if st.button("Classificar"):
        if not {"gene", "identidade", "cobertura"}.issubset(df.columns):
            st.error("O arquivo CSV deve conter as colunas: 'gene', 'identidade' e 'cobertura'.")
        else:
            pontuacao, resultado, detalhes = classificar(df, gram, base_virulencia)

            st.subheader("Resultado da análise")
            st.write(f"**Pontuação total ajustada:** {round(pontuacao, 2)}")
            st.write(f"**Classificação:** {resultado}")

            st.subheader("Detalhes por gene")
            detalhes_df = pd.DataFrame(detalhes, columns=["Gene", "Categoria", "Pontuação"])
            st.dataframe(detalhes_df)

            st.subheader("Distribuição por categoria funcional")
            categoria_df = detalhes_df.groupby("Categoria").agg({
                "Pontuação": ["count", "sum"]
            }).reset_index()
            categoria_df.columns = ["Categoria", "Número de genes", "Pontuação total"]
            st.dataframe(categoria_df)
            st.bar_chart(categoria_df.set_index("Categoria")["Número de genes"])
