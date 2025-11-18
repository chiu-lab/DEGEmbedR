# DEGEmbedR
**DEGEmbedR** is an R package for **gene-set–free**, **embedding-based** functional analysis of differentially expressed genes (DEGs) using large language model (LLM) embeddings.

Instead of relying on predefined gene sets (e.g., GO or KEGG), DEGEmbedR embeds both genes and biological functions in a **continuous semantic space**, enabling quantitative statistical assessment of **DEG–function relationships** for curated, user-defined, or LLM-generated biological functions.

---

## Installation
```r
# install.packages("remotes")
remotes::install_github("chiu-lab/DEGEmbedR")

# For vignette building (optional)
# remotes::install_github("chiu-lab/DEGEmbedR", build_vignettes = TRUE)
```
**System requirements**
- R (>= 4.2.3)
- Internet access for OpenAI API calls (only for `GeneratePathwayDescription()` and `GenerateTextEmbedding()`)
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
DEGEmbedR provides a **gene-set–free statistical framework** for evaluating DEG–function relationships using cosine similarity distributions between LLM-derived gene and function embeddings.

The package supports three complementary workflows:

1. **Analyze DEGs using built-in functional databases**  
   GO-BP, CP (BioCarta, KEGG, PID, Reactome, WikiPathways), and MOA (common mechanisms of action)
   → via `RunDEGEmbedR()`

2. **Generate functional descriptions from names using GPT-4o (optional)**  
   Example: “STING pathway in cancer immunotherapy”
   → via `GeneratePathwayDescription()`

3. **Convert custom or GPT-4o-generated functional descriptions to embeddings (optional)**
   → via `GenerateTextEmbedding()`

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

### **✔ Works with curated functions, custom names, or LLM-generated text**
You can test:
- Built-in GO, pathway, and MOA functions  
- Custom pathway names  
- LLM-generated biological hypotheses  
- Any text description you write

### **✔ Reproducible & offline**
Only text/embedding generation requires internet and OpenAI API key. All statistical tests run **offline**.

---

# Workflows

## A. Analyze DEGs Using Built-in Functional Databases

This workflow tests your DEG list against curated functional databases such as **GO Biological Processes**, **canonical pathways**, and **mechanisms of action (MOA)**.

### 1. Load example DEGs
```r
load(system.file("examples", "example.rdata", package = "DEGEmbedR"))
print(degs)
```

### 2. Run DEG–function analysis
```r
result_tb <- RunDEGEmbedR(
  degs     = degs,
  bkgs     = bkgs,     # background genes in the example dataset
  category = "GOBP"
)
```

---

## B. Generate Functional Descriptions From a Pathway Name (Optional)

Use this workflow when you want to test a biological function that is **not** part of curated databases.

### Generate description
```r
pathway_desc <- GeneratePathwayDescription(
  pathway = "STING pathway in cancer immunotherapy",
  api_key = api_key
)
```

This produces a concise, standardized biological description.

---

## C. Convert Custom or LLM-Generated Descriptions Into Embeddings (Optional)

Use this workflow to test **your own functional hypotheses**, including GPT-4o-generated summaries (from **Workflow B**) or manually written descriptions.

### 1. Convert text to embeddings
```r
embed_mat <- GenerateTextEmbedding(
  text    = pathway_desc,
  api_key = api_key
)
```

### 2. Run customized DEG–function analysis
```r
result_custom <- RunDEGEmbedR(
  degs            = degs,
  bkgs            = bkgs,
  category        = "Customized",
  embedding_input = embed_mat
)
```

---

# Function Reference

## RunDEGEmbedR()
Perform statistical DEG–function testing.

## GeneratePathwayDescription()
Generate a functional description using OpenAI **GPT-4o (gpt-4o-2024-08-06)**.

## GenerateTextEmbedding()
Convert text to embeddings using OpenAI **text-embedding-3-large**.

---

## Citation
If you use **DEGEmbedR** in your research, please cite:

> Tan, Y., Wang, L.-J., Liang, T., Lai, Y.-J., Shih, C.-H., Guo, Y., Yasaka, T. M., Tseng, G. C., & Chiu, Y.-C. (2025). _Large language model embeddings enable quantitative gene function representation and hypothesis testing._ **Under review**

You can obtain a BibTeX entry in R via:

```r
citation("DEGEmbedR")
```

---

# Contact
Please open an Issue or Pull Request on GitHub for questions or contributions.
