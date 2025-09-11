import streamlit as st
import pandas as pd
from datetime import datetime

st.set_page_config(page_title="Analisador de Patogenicidade", layout="wide")

# -----------------------------
# Função para carregar a base
# -----------------------------
def carregar_base():
    try:
        df = pd.read_csv("base_virulencia.csv", sep=None, engine="python", encoding="utf-8")
    except Exception:
        df = pd.read_csv("base_virulencia.csv", sep=None, engine="python", encoding="latin-1")

    # 🔹 Linha de debug removida

    if "gene" not in df.columns:
        for col in df.columns:
            if col.strip().lower() == "gene":
                df.rename(columns={col: "gene"}, inplace=True)
                break

    if "gene" not in df.columns:
        st.error("❌ A base não possui coluna 'gene'. Verifique o cabeçalho do arquivo.")
        return pd.DataFrame(columns=["gene", "categoria", "gram", "peso", "evidencia"])

    df["gene"] = df["gene"].astype(str).str.strip().str.lower()
    return df

# -----------------------------
# Backup automático da base
# -----------------------------
def backup_base():
    data_str = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    try:
        df = pd.read_csv("base_virulencia.csv", sep=None, engine="python", encoding="utf-8")
    except Exception:
        df = pd.read_csv("base_virulencia.csv", sep=None, engine="python", encoding="latin-1")
    df.to_csv(f"base_virulencia_backup_{data_str}.csv", index=False)

# backup_base()  # Descomente se quiser backup automático

# -----------------------------
# Interface
# -----------------------------
st.title("🧬 Analisador de Patogenicidade Bacteriana")
st.write("Envie um arquivo CSV para iniciar a análise.")

base = carregar_base()
# 🔹 Linha de debug removida

# Modelo CSV para download
modelo_df = pd.DataFrame({
    "gene": ["geneA", "geneB", "geneC"],
    "identidade": [99.5, 87.3, 92.1],
    "cobertura": [100, 96, 98]
})
st.download_button(
    label="📥 Baixar modelo CSV de entrada",
    data=modelo_df.to_csv(index=False).encode("utf-8"),
    file_name="modelo_entrada.csv",
    mime="text/csv"
)

# Upload do arquivo de entrada
arquivo = st.file_uploader(
    "Envie seu arquivo CSV com colunas: gene, identidade, cobertura",
    type=["csv"]
)

if arquivo:
    df_input = pd.read_csv(arquivo)
    df_input["gene"] = df_input["gene"].astype(str).str.strip().str.lower()

    st.write(f"📄 Genes no arquivo enviado: {len(df_input)}")

    if st.button("Classificar"):
        resultados = []
        for _, row in df_input.iterrows():
            gene = row["gene"]
            identidade = row["identidade"]
            cobertura = row["cobertura"]

            match = base[base["gene"] == gene]
            if not match.empty and identidade >= 95 and cobertura >= 95:
                info = match.iloc[0]
                probabilidade = round((identidade / 100) * (cobertura / 100) * (info["peso"] / 3), 2)
                resultados.append({
                    "gene": gene,
                    "categoria": info["categoria"],
                    "gram": info["gram"],
                    "peso": info["peso"],
                    "evidencia": info["evidencia"],
                    "identidade": identidade,
                    "cobertura": cobertura,
                    "probabilidade": probabilidade
                })

        if resultados:
            df_resultados = pd.DataFrame(resultados)
            df_resultados["peso"] = pd.to_numeric(df_resultados["peso"], errors="coerce")
            pontuacao_total = df_resultados["peso"].sum()

            if pontuacao_total >= 30:
                classificacao = "Alta probabilidade de patogenicidade"
            elif pontuacao_total >= 20:
                classificacao = "Média probabilidade de patogenicidade"
            else:
                classificacao = "Baixa probabilidade de patogenicidade"

            st.markdown("## 📊 Resultado da análise")
            st.write(f"**Pontuação total ajustada:** {pontuacao_total}")
            st.write(f"**Classificação:** {classificacao}")

            # --- Filtros interativos ---
            st.markdown("## 🔍 Filtrar resultados")
            genes_unicos = sorted(df_resultados["gene"].unique())
            filtro_gene = st.multiselect("Filtrar por gene:", options=genes_unicos)

            categorias_unicas = sorted(df_resultados["categoria"].unique())
            filtro_categoria = st.multiselect("Filtrar por categoria:", options=categorias_unicas)

            gram_unicos = sorted(df_resultados["gram"].unique())
            filtro_gram = st.multiselect("Filtrar por Gram:", options=gram_unicos)

            df_filtrado = df_resultados.copy()
            if filtro_gene:
                df_filtrado = df_filtrado[df_filtrado["gene"].isin(filtro_gene)]
            if filtro_categoria:
                df_filtrado = df_filtrado[df_filtrado["categoria"].isin(filtro_categoria)]
            if filtro_gram:
                df_filtrado = df_filtrado[df_filtrado["gram"].isin(filtro_gram)]

            st.markdown("## 📋 Tabela filtrada")
            st.dataframe(df_filtrado, use_container_width=True)

        else:
            st.warning("Nenhum gene foi classificado com os critérios atuais.")
