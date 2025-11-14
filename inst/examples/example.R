####Example code of DEGEmbedR
remotes::install_github("chiu-lab/DEGEmbedR")

library(DEGEmbedR)
api_key <-"YOUR_OPENAI_API_KEY"
####1. Generate pathway descriptions with GPT-4o####
#--Single pathway--#
# Generate a single pathway description
desc_wnt <- Generate_PathwayDescription(
  pathway = "Wnt Signaling Pathway",
  api_key
)

#--Multiple pathways--#
# Define the pathway names you want to analyze
pathways <- c("Wnt Signaling Pathway", "Apoptosis", "MAPK Pathway")

# Generate descriptions for each pathway using lapply()
pathway_desc_list <- lapply(pathways, function(pw) {
  Generate_PathwayDescription(
    pathway = pathways,
    api_key
  )
})

# Convert the list to a named character vector
pathway_desc <- setNames(unlist(pathway_desc_list), pathways)

# Check results
names(pathway_desc)
pathway_desc[1:2]

####2. Convert text descriptions to embeddings####
# Convert the pathway descriptions into embedding vectors
embed_mat <- Generate_TextEmbedding(
  text = pathway_desc,
  api_key
)

# Check structure
dim(embed_mat)
rownames(embed_mat)

####3. Perform enrichment or similarity testing####
# Select random genes as an example
load(system.file("examples","example.rdata",package = "DEGEmbedR"))
print(degs)

# Run the comparison using the custom embeddings
result_tb <- CompareGeneSetEmbeddings(
  degs = degs,
  category = "Customized",
  embedding_input = embed_mat
)

# View top pathways
head(result_tb)

