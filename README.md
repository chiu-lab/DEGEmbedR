# DEGEmbedR
DEGEmbedR is an easy-to-use R package for hypothesis generation and statistical testing of gene functions using embeddings from large language models. It supports curated and custom embeddings, enabling enrichment analysis and discovery of novel LLM-inferred biological functions.

---

## Installation
```r
# install.packages("remotes")
remotes::install_github("chiu-lab/DEGEmbedR")

# For vignette building (optional)
# remotes::install_github("chiu-lab/DEGEmbedR", build_vignettes = TRUE)
```
**System requirements**
- R (>= 4.1)
- Internet access for OpenAI API calls (only for `Generate_*` functions)
- Suggested: macOS/Linux/Windows with 8GB+ RAM

**Key dependencies (installed automatically)**: `tibble`, `stringr`, `lsa`, `effsize`, `httr`, `jsonlite`.

---

## Data assets bundled with the package
The package ships curated cosine similarity matrices and embeddings covering:
- **GO BP** gene–term similarities
- **MSigDB C2** sub-collections: BIOCARTA, KEGG, PID, REACTOME, WP
- **MOA** gene–MOA/function similarities
- **Gene embeddings** and the built‑in **18K gene universe**

These are loaded on demand by `CompareGeneSetEmbeddings()`; users do **not** need to load them manually.

---
## Quick start
```r
library(DEGEmbedR)

# Example data included with the package
load(system.file("examples", "example.rdata", package = "DEGEmbedR"))
# Objects: `degs`, `embed_mat` (mock custom embeddings for illustration)

# 1) Compare DEGs vs background against curated pathways (e.g., GO BP)
res_gobp <- CompareGeneSetEmbeddings(
  degs     = degs,
  category = "GOBP"
)
head(res_gobp)

# 2) Use your own term embeddings ("Customized")

res_custom <- CompareGeneSetEmbeddings(
  degs            = degs,
  category        = "Customized",
  embedding_input = embed_mat
)
head(res_custom)
```

> A tab‑delimited results file named like `result_YYYY-MM-DD-HHMMSS.txt` is also written to your working directory.

---

## Typical workflows

### A) Hypothesis testing when no curated gene set exists
1. Draft a brief description for a pathway/MOA of interest (or let `Generate_PathwayDescription()` draft one).
2. Convert the description(s) to embeddings with `Generate_TextEmbedding()`.
3. Run `CompareGeneSetEmbeddings(..., category = "Customized", embedding_input = your_matrix)` to test enrichment (DEGs > background).

```r
# 1) Generate a description with OpenAI (requires API key)
api_key = "<your OPENAI_API_KEY >"
text <- Generate_PathwayDescription(
  pathway = "STING pathway",
  api_key = api_key
)

# 2) Turn text into an embedding matrix (rows = pathway, cols = 3072)
emb <- Generate_TextEmbedding(
  text   = text,
  api_key = api_key
)

# 3) Test whether DEGs show higher similarity than background
res <- CompareGeneSetEmbeddings(
  degs            = degs,
  category        = "Customized",
  embedding_input = emb
)

```

---

## Citation
If you use **DEGEmbedR** in your research, please cite:

> Tan, Y., Wang, L.-J., Liang, T., Lai, Y.-J., Shih, C.-H., Guo, Y., Yasaka, T. M., Tseng, G. C., & Chiu, Y.-C. (2025). _Large language model embeddings enable quantitative gene function representation and hypothesis testing._ **Under review**

You can obtain a BibTeX entry in R via:

```r
citation("DEGEmbedR")
```

---

