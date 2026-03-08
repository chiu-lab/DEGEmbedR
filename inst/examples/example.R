#### Example workflow for DEGEmbedR ####

# Install DEGEmbedR from GitHub
remotes::install_github("chiu-lab/DEGEmbedR")

# Set working directory if needed
setwd("path/to/your/work/directory")

library(DEGEmbedR)

# Provide your OpenAI API key, required only when using GSAI
api_key <- "YOUR_OPENAI_API_KEY"

# 1. Run DEG-function analysis using the default built-in background (approximately 18k human protein-coding genes) and
#    curated functional annotations(e.g., GOBP, KEGG).
# Example DEG list: NeST: 79 (ATM-Dependent DNA Repair), sourced from https://github.com/ncbi-nlp/GeneAgent/tree/main/Datasets/NeST
NeST79 <- c('ATM', 'AURKA', 'BARD1', 'BLM', 'BRCA1',
            'BRCA2', 'BRIP1', 'BUB1B', 'CDC73', 'CHEK1',
            'FANCA', 'FANCD2', 'MDC1', 'MDM2', 'MRE11',
            'MSH6', 'NBN', 'PALB2', 'POLH', 'RAD50',
            'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'TOP2A',
            'TP53', 'WRN', 'XRCC2', 'XRCC3')

result_tb1 <- RunDEGEmbedR(
  degs = NeST79,
  category = "GOBP"
)


# 2. Generate and evaluate AI-derived functional hypotheses against DEGs (GSAI mode).
# Example DEG list: NeST: 105 (Ubiquitin Regulation of p53 Activity), sourced from https://github.com/ncbi-nlp/GeneAgent/tree/main/Datasets/NeST
NeST105 <- c('CUL3', 'ELOC', 'FBXW7', 'HSP90AA1', 'MDM2', 'SKP1', 'STK11', 'TNF', 'VHL')
result_tb2 <- RunDEGEmbedR(
  degs = NeST105,
  category = "GSAI",
  api_key = api_key
)

# View top results
head(result_tb1)
head(result_tb2)

