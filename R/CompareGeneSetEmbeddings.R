#' Compare DEG vs. background cosine similarities across pathways/MOAs
#'
#' @description
#' For each pathway (or user-supplied term embedding), this function compares the
#' cosine similarity distribution of differentially expressed genes (DEGs) to that
#' of background genes and summarizes the statistical evidence for enrichment.
#'
#' The function returns the one-tailed Wilcoxon rank-sum test p-value
#' (testing whether DEGs > background), median cosine similarities for DEGs
#' and background genes, their median difference, Cliff’s delta with a 95% confidence
#' interval, and the top 10 DEGs with the highest cosine similarity.
#'
#' @details
#' This implements the analysis used in the paper's applications (e.g., testing
#' pathway/MOA hypotheses for STING pathway and drug mechanisms). It supports
#' curated similarity matrices distributed with the package (GO BP, MSigDB C2
#' subcollections BIOCARTA/KEGG/PID/REACTOME/WP, and MOA) as well as a
#' "Customized" mode that accepts user-supplied pathway/term embeddings.
#'
#'
#' @param degs Differentially expressed genes (DEGs). After intersecting with
#'   the built-in set of 18K genes, the number of matched DEGs must be in 15~500.
#' @param bkgs Background gene symbols (optional). If
#'   \code{NULL},  the built-in set of 18K genes is used.
#' @param category Select pathway/MOAs/customized categories to analyze
#'   one of \code{"GOBP"}, \code{"C2CP_all"}, \code{"BIOCARTA"}, \code{"KEGG"},
#'   \code{"PID"}, \code{"REACTOME"}, \code{"WP"}, \code{"MOA"}, or
#'   \code{"Customized"} (case-insensitive).
#' @param embedding_input Numeric matrix or data frame containing pathway or term
#'   embeddings (required when \code{category = "Customized"}). Rows represent terms,
#'   columns represent embedding dimensions (length = 3072). Row names are used as term labels.
#' @param output Logical. The result will be saved as a .txt file (Default: TRUE).
#'
#' @return A \link[tibble]{tibble} with one row per pathway/term:
#' \describe{
#'   \item{\code{name}}{Pathway/term name}
#'   \item{\code{p_value_MWN_one_tailed}}{One-tailed Wilcoxon p-value (DEGs > background)}
#'   \item{\code{median_cosine_similarity_degs}}{Median similarity in DEGs}
#'   \item{\code{median_cosine_similarity_bkgs}}{Median similarity in background}
#'   \item{\code{diff_cosine_similarity.}}{Median(DEGs) - Median(Background)}
#'   \item{\code{cliffs_delta}}{Cliff's delta effect size}
#'   \item{\code{cliffs_delta_ci_95}}{95% CI for Cliff's delta}
#'   \item{\code{cliffs_delta_magnitude}}{Magnitude category}
#'   \item{\code{top10_degs_with_highest_cosine_similarity}}{Top 10 DEGs with values}
#' }
#'
#' @section A tab-delimited results file named like:
#'  \code{"result_YYYY-MM-DD-HHMMSS.txt"} is written to the working directory.
#'
#' @examples
#' \dontrun{
#' # Customized example with mock embeddings (no on-disk .Rdata needed)
#' load(system.file("examples","example.rdata",package = "DEGEmbedR"))
#' head(embed_mat)
#' length(degs)
#' res <- CompareGeneSetEmbeddings(
#'   degs = degs,
#'   category = "Customized",
#'   embedding_input = embed_mat
#' )
#' head(res)
#' }
#'
#' @seealso
#' \code{\link[stats]{wilcox.test}}, \code{\link[effsize]{cliff.delta}},
#' \code{\link[tibble]{tibble}}, \code{\link[stringr]{str_extract}}
#' @importFrom stats wilcox.test
#' @importFrom stringr str_extract
#' @keywords enrichment embeddings similarity statistics effect-size
#' @export



CompareGeneSetEmbeddings <- function(degs,
                                   bkgs=NULL,
                                   category = c("GOBP","C2CP_all","BIOCARTA", "KEGG","PID","REACTOME", "WP", "MOA","Customized"),
                                   embedding_input=NULL,
                                   output = TRUE){


###Load data###
  bp_dt    <- readRDS(data_path("BP_15-500_similarity_text-embedding-3-large_2024110801.rds"))
  cp_dt    <- readRDS(data_path("CP_15-500_similarity_text-embedding-3-large_2024110801.rds"))
  moa_gene_dt <- readRDS(data_path("gene_geneset_moa_drug_function_similarity_text-embedding-3-large_2025031701.rds"))
  gene_dt  <- readRDS(data_path("gene_embedding_2024110801.rds"))
  gene_list <- readRDS(data_path("gene_list.rds"))

###Require R packages###
if (!requireNamespace("lsa", quietly = TRUE))     stop("Package 'lsa' is required.")
if (!requireNamespace("tibble", quietly = TRUE))  stop("Package 'tibble' is required.")
if (!requireNamespace("effsize", quietly = TRUE)) stop("Package 'effsize' is required.")
if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")

library(effsize)
library(tibble)
library(stringr)
library(lsa)

###Check gene symbol###
##bkgs##
if (is.null(bkgs)) {
  match_bkgs <- gene_list
} else {
  # setdiff+intersect are C-optimized; avoids which()/sum()>0 branching
  match_bkgs <- intersect(setdiff(bkgs, degs), gene_list)
}
message(sprintf("There are %d matched genes with build-in gene lists", length(match_bkgs)))



##degs##
match_degs <- intersect(degs, gene_list)
if(length(match_degs) < 15 | length(match_degs) > 500){
  stop("Insufficient number of genes to match")
}

###Check category#
if(is.null(category)){
  stop("Category is required")
}

###Customized###
##Generate cosine similarity table

if(toupper(category) == toupper("Customized")){
  if(is.null(embedding_input)){
    stop("Missing embedding input")
  }

  #Cosine similarity calculation
  # normalize rows
  gene_norm <- gene_dt / sqrt(rowSums(gene_dt * gene_dt))
  pth_norm  <- embedding_input / sqrt(rowSums(embedding_input * embedding_input))
  gene_norm <- as.matrix(gene_norm)
  mode(gene_norm) <- "numeric"

  pth_norm <- as.matrix(pth_norm)
  mode(pth_norm) <- "numeric"
  costom_matrix <- gene_norm %*% t(pth_norm)

}


###Output table###
results <- tibble::tibble(
  name                                      = character(),
  p_value_MWN_one_tailed                    = numeric(),
  median_cosine_similarity_degs             = numeric(),
  median_cosine_similarity_bkgs             = numeric(),
  diff_cosine_similarity.                   = numeric(),
  cliffs_delta                              = numeric(),
  cliffs_delta_ci_95                        = character(),
  cliffs_delta_magnitude                    = character(),
  top10_degs_with_highest_cosine_similarity = character()
)


###Select consine similarity matrix###
con_sim_mtrx <- switch(toupper(category),
                       GOBP       = bp_dt,
                       C2CP_all   = cp_dt,
                       BIOCARTA   = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "BIOCARTA")],
                       KEGG       = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "KEGG")],
                       PID        = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "PID")],
                       REACTOME   = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "REACTOME")],
                       WP         = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "WP")],
                       MOA        = moa_gene_dt,
                       CUSTOMIZED = costom_matrix,
                       stop("Unknown category")
)

tb <-  as.data.frame(con_sim_mtrx[rownames(con_sim_mtrx) %in% match_bkgs, , drop = FALSE])
for (i in 1:ncol(tb)) {


cos_sim_degs <- tb[match_degs, colnames(con_sim_mtrx)[i]]
names(cos_sim_degs) <- match_degs
cos_sim_bkgs <- tb[!row.names(tb) %in% match_degs, colnames(con_sim_mtrx)[i]]
###calculate pvalue###
   message("Calculating p-values using Wilcoxon rank-sum test")
   stat_wilcox <- wilcox.test(x = cos_sim_degs ,
                           y=cos_sim_bkgs,
                           alternative = "greater")

###Top10 genes with cosine similarity###
top10 <- head(sort(cos_sim_degs , decreasing = TRUE), 10)
top10_label <- paste0(names(top10), "(", round(top10, 4), ")", collapse = ", ")

###Median for degs and bkgs###
median_degs <- median(cos_sim_degs )
median_bkgs <- median(cos_sim_bkgs)

###difference between degs and bkgs###
diff_degs_minus_bkgs <- median_degs - median_bkgs

###Cliff delta effect size###
cliff_effect_size <- effsize::cliff.delta(cos_sim_degs ,
                                 cos_sim_bkgs)

results[i,] <- list(colnames(con_sim_mtrx)[i],
                    as.numeric(stat_wilcox$p.value),
                    median_degs,
                    median_bkgs,
                    diff_degs_minus_bkgs,
                   as.numeric( cliff_effect_size$estimate),
                  paste0("[",paste(round(cliff_effect_size$conf.int,2),collapse = ","),"]"),
                  cliff_effect_size$magnitude,
                  top10_label
                  )
}


results <- results[order(results$p_value_MWN_one_tailed,decreasing = F),]
if(output){
    write.table(results,file = paste("result", format(Sys.time(), "%Y-%m-%d-%H%M%S.txt"),sep = "_"),
                sep = "\t", col.names = T, row.names = F, quote = F)
}
  return(results)
}
