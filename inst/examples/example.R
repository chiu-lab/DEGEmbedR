#### Example workflow for DEGEmbedR ####

# Install DEGEmbedR from GitHub
remotes::install_github("chiu-lab/DEGEmbedR")

# Set working directory if needed
setwd("path/to/your/work/directory")

library(DEGEmbedR)

# Provide your OpenAI API key
api_key <- "YOUR_OPENAI_API_KEY"

#### A. Analyze built-in functional description databases ####

# Load example DEGs from visugromab treatment
# (Melero et al., Nature 2024; case study described in the manuscript)
load(system.file("examples", "example.rdata", package = "DEGEmbedR"))

# Inspect input genes
print(degs)
length(degs)

# 1. Run DEG-function analysis using the default built-in background
#    (approximately 18k human protein-coding genes)
result_tb1 <- RunDEGEmbedR(
  degs = degs,
  category = "GOBP",
  api_key = api_key
)

# 2. Run DEG-function analysis using a user-supplied background gene set
length(bkgs)
result_tb2 <- RunDEGEmbedR(
  degs = degs,
  bkgs = bkgs,
  category = "GOBP",
  api_key = api_key
)

# 3. Run analysis using an AI-derived functional hypothesis (GSAI mode)
result_tb3 <- RunDEGEmbedR(
  degs = degs,
  category = "GSAI",
  api_key = api_key
)

# View top results
head(result_tb1)
head(result_tb2)
head(result_tb3)
