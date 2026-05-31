# DEGEmbedR
**DEGEmbedR** is an R package for **gene-set–free**, **embedding-based** functional analysis of differentially expressed genes (DEGs) using large language model (LLM) embeddings.

Instead of relying on predefined gene sets (e.g., GO or KEGG), DEGEmbedR embeds both genes and biological functions in a **continuous semantic space**, enabling quantitative statistical assessment of **DEG–function relationships** for curated or LLM-generated biological functions.

<p align="center">
  <img src="man/figures/DEGEmbedR_workflow.png" alt="DEGEmbedR workflow" width="800">
</p>

---

## Installation
**⚠️ We have recently updated to version `1.2.0`. Please re-install if you have previously installed it:**
```r
# install.packages("remotes")
remotes::install_github("chiu-lab/DEGEmbedR")

# For vignette building (optional)
# remotes::install_github("chiu-lab/DEGEmbedR", build_vignettes = TRUE)
```
**System requirements**
- R (>= 4.2.3)
- Internet access for OpenAI API calls (required only for the LLM generative mode `category = "GSAI"`)
- Suggested: macOS/Linux/Windows with 8GB+ RAM

**Key dependencies (installed automatically):**  `tibble`, `stringr`, `lsa`, `effsize`, `httr`, `jsonlite`.

---

## Quick Start with the Example Script (Strongly Recommended)

An example script named **`example.R`** is included with the package under the `examples/` directory.  
You can open it in your R editor with:

```r
file.edit(system.file("examples", "example.R", package = "DEGEmbedR"))
```

---

# Overview
The main function `RunDEGEmbedR()` provides a **gene-set–free statistical framework** for evaluating DEG-function relationships using gene–function cosine similarity distributions between LLM-derived gene and function embeddings. It supports analysis against curated functional annotations as well as AI-generated *de novo* functional hypotheses. **Prebuilt background genes** (approximately 18,000 human protein-coding genes) are included for contrast with the DEGs. Users can also upload customized background genes if the DEGs are derived from a targeted panel, or other specialized datasets.

The package supports two complementary workflows, both accessible via the main function `RunDEGEmbedR`:

1. **Analyze DEGs using built-in functional annotation databases**  
   GO-BP (Biological Processes) and CP (BioCarta, KEGG, PID, Reactome, WikiPathways).

2. **Generate *de novo* functional hypotheses of DEGs using an LLM and analyze their statistical significance (requires OpenAI API key)**  
   GSAI (Hu *et al*. *Nature Methods* 2025), re-implemented using GPT-4o.

---

# Key Features

### **✔ Gene-set–free functional modeling**
Biological functions are represented with **continuous LLM embeddings**, not fixed gene lists.

### **✔ Unified statistical DEG–function test**
DEGEmbedR reports:
- Median cosine similarity (DEGs vs background)
- Difference in medians
- One-tailed p-value from the Wilcoxon rank sum test
- Cliff’s delta and 95% CI
- Top DEGs driving the signal

### **✔ Works with curated functions and LLM-generated functional hypotheses**
You can test:
- Built-in biological processes and pathways for integration and comparison with existing gene set analysis methods
- LLM-generated biological hypotheses for *de novo* functions

### **✔ Reproducible & offline**
Only the LLM generation mode (`category = "GSAI"`) requires internet and OpenAI API key. All statistical tests run **offline**.

---

# Workflows

## A. Analyze DEGs Using Built-in Functional Databases

This workflow tests your DEG list against curated functional databases such as **GO Biological Processes** and **canonical pathways**.

```r
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
```

---

## B. Generate and evaluate AI-derived *de novo* functional hypotheses against DEGs (OpenAI API key required)

Use this workflow when you want to test a biological function that is **not** part of curated databases.

### Generate description
```r
# Example DEG list: NeST: 105 (Ubiquitin Regulation of p53 Activity), sourced from https://github.com/ncbi-nlp/GeneAgent/tree/main/Datasets/NeST
NeST105 <- c('CUL3', 'ELOC', 'FBXW7', 'HSP90AA1', 'MDM2', 'SKP1', 'STK11', 'TNF', 'VHL')

result_tb2 <- RunDEGEmbedR(
  degs = NeST105,
  category = "GSAI",
  api_key = api_key
)
```

---

# Function Reference

## RunDEGEmbedR()
Main function to perform statistical DEG–function testing using curated or *de novo* generated functional annotations.

---

## Citation
If you use **DEGEmbedR** in your research, please cite:

> Tan, Y., Wang, L.-J., Liang, T., Lai, Y.-J., Shih, C.-H., Guo, Y., Yasaka, T. M., Tseng, G. C., & Chiu, Y.-C. (2026). _An embedding-based framework enables statistical testing of gene-set function hypotheses inferred by large language models._ **Under review**

You can obtain a BibTeX entry in R via:

```r
citation("DEGEmbedR")
```

---

# Contact
Please open an Issue or Pull Request on GitHub for questions or contributions.
