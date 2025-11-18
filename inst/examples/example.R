####Example code of DEGEmbedR
remotes::install_github("chiu-lab/DEGEmbedR")
setwd("path/to/your/work/directory")

library(DEGEmbedR)

####A. Analyze Built-in Functional Description Databases####
# Load DEGs from visugromab treatment from Melero et al. Nature 2024. (case study reported in the manuscript)
load(system.file("examples","example.rdata",package = "DEGEmbedR"))
print(degs)

# Run DEG-function analysis analysis all GO-BP terms, compared with the background genes probed in the IO panel
result_tb <- RunDEGEmbedR(
  degs = degs,
  bkgs = bkgs,
  category = "GOBP"
)

# Run DEG-function analysis for the STING pathway (pre-calculated and included in examples; see sections B and C below),
# compared with the background genes probed in the IO panel (Fig. 6B in the manuscript)
result_tb2 <- RunDEGEmbedR(
  degs = degs,
  bkgs = bkgs,
  category = "Customized",
  embedding_input = embed_mat
)

# Run DEG-function analysis for the STING pathway (pre-calculated and included in examples; see sections B and C below),
# compared with all built-in 18K background genes (all protein-coding genes in NCBI) (Fig. 6B in the manuscript)
result_tb3 <- RunDEGEmbedR(
  degs = degs,
  category = "Customized",
  embedding_input = embed_mat
)

# View top pathways
head(result_tb)
head(result_tb2)
head(result_tb3)



####B. Creating and Analyze Custom Functional Descriptions from Names####
# OpenAI key is required to generate descriptions and embeddings
api_key <-"YOUR_OPENAI_API_KEY"

# Single pathway - generate a single pathway description
desc_sting <- GeneratePathwayDescription(
  pathway = "STING Pathway in Cancer Immunotherapy",
  api_key
)


# Multiple pathways - list the pathway names you want to analyze
pathways <- c("STING Pathway in Cancer Immunotherapy", "Wnt Signaling Pathway", "Apoptosis", "MAPK Pathway")

# Generate descriptions for each pathway using lapply()
pathway_desc_list <- lapply(pathways, function(pw) {
  GeneratePathwayDescription(
    pathway = pw,
    api_key
  )
})

# Convert the list to a named character vector
pathway_desc <- setNames(unlist(pathway_desc_list), pathways)

# Check results
names(pathway_desc)
pathway_desc[1:2]



####C.Convert GPT-generated (from section B) or custom functional descriptions to embeddings####
embed_mat <- GenerateTextEmbedding(
  text = pathway_desc,
  api_key
)

# Check structure
dim(embed_mat)
rownames(embed_mat)

# Then run RunDEGEmbedR with category = "Customized" and embedding_input = embed_mat (refer to section A)

