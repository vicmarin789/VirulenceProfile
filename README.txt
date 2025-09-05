## 📩 Solicitação de uso oficial

Para acesso à **base completa de dados** e utilização da ferramenta em análises oficiais, entre em contato:

**Victor Augustus Marin**  
📧 **Email:** [victor.marin@unirio.br](mailto:victor.marin@unirio.br)

---

# 🧬 Analisador de Patogenicidade Bacteriana

Ferramenta para análise e classificação de genes de virulência bacteriana com base em identidade e cobertura, integrando dados automáticos da UniProt e base manual curada.

Desenvolvida para auxiliar na triagem e avaliação de risco microbiológico em diferentes contextos, com potencial de aplicação por órgãos de fiscalização, vigilância e pesquisa.

---

## 🚀 Principais Aplicações

- **MAPA** – Monitoramento de patógenos em alimentos e insumos agropecuários.
- **Anvisa** – Apoio em análises de risco e investigações sanitárias.
- **Ibama** – Avaliação de riscos microbiológicos em fauna, flora e ambientes sensíveis.

---

## 📂 Estrutura da Ferramenta

- **Base automática**: obtida da UniProt (genes com palavra-chave de virulência).
- **Base manual**: curada por especialista, com categorias funcionais e pesos ajustados.
- **Classificação**: baseada em identidade (%) e cobertura (%) dos genes encontrados.

---

## 🖥️ Como Testar

1. **Baixe o modelo de entrada** disponível no app.
2. Preencha com os genes identificados, identidade e cobertura.
3. Carregue o arquivo no aplicativo Streamlit.
4. Selecione o tipo de bactéria (Gram-negativa ou Gram-positiva).
5. Visualize a pontuação e a classificação final.

> ⚠️ **Atenção:** A versão pública não contém a base completa. Para uso oficial, solicite acesso conforme instruções no início deste documento.

---

## 📊 Exemplo de Saída

| Gene  | Categoria                                   | Pontuação |
|-------|---------------------------------------------|-----------|
| invA  | Estrutura de sistema de secreção tipo III   | 4.5       |
| prgH  | Estrutura de sistema de secreção tipo III   | 4.5       |
| sipB  | Efetor de sistema de secreção tipo III      | 4.5       |
| ssaV  | Estrutura de sistema de secreção tipo III   | 4.5       |

---

## 📜 Licença

Uso restrito para fins de pesquisa e desenvolvimento.  
Para aplicação em análises oficiais, é necessário contato prévio e autorização.

---
