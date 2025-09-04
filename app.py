import streamlit as st
import pandas as pd
import requests

# ---------------------------
# Função para atualizar a base de virulência
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
    for org_name, org_id, gram in organismos:
        query = f"(organism_id:{org_id}) AND (keyword:KW-0800)"
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=tsv&fields=accession,gene_names,organism_name"
        r = requests.get(url)
        if r.status_code != 200:
            continue
        lines = r.text.strip().split("\n")[1:]
        for line in lines:
            parts = line.split("\t")
            gene_name = parts[1].split(" ")[0] if parts[1] else "desconhecido"
            genes.append({
                "gene": gene_name,
                "categoria": "desconhecida",
                "gram": gram,
                "peso": 3,
                "evidencia": True
            })

    df = pd.DataFrame(genes).drop_duplicates(subset="gene")
    df.to_csv("base_virulencia.csv", index=False)
    return len(df)

# ---------------------------
# Carrega base externa de genes virulentos
# ---------------------------
@st.cache_data
def carregar_base_virulencia():
    return pd.read_csv("base_virulencia.csv")

# ---------------------------
# Função de classificação refinada
# ---------------------------
def classificar(df, gram="neg", base=None):
    pontuacao = 0
    detalhes = []

    for _, fator in df.iterrows():
        gene = fator["gene"]
        identidade = fator["identidade"]
        cobertura = fator["cobertura"]

        info = base.query(f"(gene == @gene) and (gram == @gram or gram == 'ambos')")
        if not info.empty:
            linha = info.iloc[0]
            if identidade >= 85 and cobertura >= 95:
                peso = linha["peso"]
                if linha["evidencia"]:
                    peso *= 1.5
                pontuacao += peso
                detalhes.append((gene, linha["categoria"], round(peso, 2)))
            else:
                detalhes.append((gene, linha["categoria"], 0))
        else:
            detalhes.append((gene, "não classificado", 0))

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

# Botão para atualizar a base
if st.button("🔄 Atualizar base de virulência"):
    total = atualizar_base_virulencia()
    st.success(f"Base atualizada com {total} genes virulentos.")

# Botão para baixar modelo de entrada
modelo_csv = "gene,identidade,cobertura\ninvA,98.5,100\nfimH,92.3,97\nstn,88.0,95\n"
st.download_button("📥 Baixar modelo de entrada (.csv)", modelo_csv, file_name="exemplo_entrada.csv")

# Carrega a base atualizada
base_virulencia = carregar_base_virulencia()

uploaded_file = st.file_uploader("Carregar arquivo CSV", type=["csv"])
gram = st.radio("Tipo de bactéria:", ("Gram-negativa", "Gram-positiva"))
gram = "neg" if gram == "Gram-negativa" else "pos"

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

            # Visualização por categoria funcional
            st.subheader("Distribuição por categoria funcional")
            categoria_df = detalhes_df.groupby("Categoria").agg({
                "Pontuação": ["count", "sum"]
            }).reset_index()
            categoria_df.columns = ["Categoria", "Número de genes", "Pontuação total"]
            st.dataframe(categoria_df)
            st.bar_chart(categoria_df.set_index("Categoria")["Número de genes"])