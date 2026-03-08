#' Compare DEG vs. background cosine similarities across curated or de novo biological functions or pathways
#'
#' @description
#' For each built-in or LLM-derived de novo biological function, this function compares the
#' gene-function cosine similarity distributions of differentially expressed genes (DEGs)
#' and background genes, summarizing the statistical evidence for DEG-function relationships.
#'
#' The function returns the one-tailed Wilcoxon rank-sum test p-value (testing whether DEGs >
#' background), the median cosine similarities for DEGs and background genes, their median
#' difference, Cliff’s delta with a 95% confidence interval, and the top 10 DEGs with the highest
#' similarity to the function.
#'
#' @details
#' \code{RunDEGEmbedR()} is the main function that implements the analysis workflow described
#' in the DEGEmbedR vignette and manuscript. It supports curated functional collections included
#' with the package, including MSigDB GO Biological Processes (BP) and MSigDB C2 pathway collections
#'  (BIOCARTA, KEGG, PID, REACTOME, and WikiPathways). In addition, the GSAI mode generates AI-derived
#'  de novo functional hypotheses from the input DEGs and then evaluates the similarity between the
#'  de novo functions to DEGs based on embeddings to statistically validate these functional hypotheses.
#'
#'
#' @param degs Character. Differential expressed genes (DEGs). After intersecting with
#'   the built-in gene universe (~18,000 genes), the number of matched DEGs must be between 15 and 500.
#' @param bkgs Character. Optional customized background genes such as targeted arrays. If `NULL`, a built-in gene
#'   universe (~18k human protein-coding genes from the NCBI Gene database) is used.
#' @param category Character. Functional database or mode to use. Must be one of:
#'   \code{"GOBP"}, \code{"C2CP_all"}, \code{"BIOCARTA"}, \code{"KEGG"}, \code{"PID"},
#'   \code{"REACTOME"}, \code{"WP"}, or \code{"GSAI"} (case-insensitive).
#' @param api_key Character (optional). OpenAI API key (https://openai.com/api/). Required only for the `GSAI` mode.
#' @param output Logical. Whether to save the result as a time-stamped tab-delimited \code{.txt} file.
#'   Default: \code{TRUE}.
#'
#' @return A \link[tibble]{tibble} with one row per pathway or term, including:
#' \describe{
#'   \item{\code{name}}{Function or pathway name}
#'   \item{\code{p_value_MWN_one_tailed}}{One-tailed Wilcoxon p-value (DEGs > background)}
#'   \item{\code{median_cosine_similarity_degs}}{Median gene-function cosine similarity among DEGs}
#'   \item{\code{median_cosine_similarity_bkgs}}{Median gene-function cosine similarity among background genes}
#'   \item{\code{diff_cosine_similarity}}{Difference: median(DEGs) – median(background)}
#'   \item{\code{cliffs_delta}}{Cliff's delta effect size}
#'   \item{\code{cliffs_delta_ci_95}}{95% confidence interval for Cliff's delta}
#'   \item{\code{cliffs_delta_magnitude}}{Magnitude category (e.g., negligible, small, medium, large)}
#'   \item{\code{top10_degs_with_highest_cosine_similarity}}{Top 10 DEGs ranked by cosine similarity}
#' }
#'
#' @section Output File:
#' A tab-delimited results file named like \code{"result_YYYY-MM-DD-HHMMSS.txt"} will be written
#' to the working directory if \code{output = TRUE}.
#'
#' @examples
#' \dontrun{
#' api_key <- "YOUR_OPENAI_API_KEY"
#'
#' # (A) Using a curated functional annotation collection
#' # NeST: 79 (ATM-Dependent DNA Repair)
#' NeST79 <- c('ATM', 'AURKA', 'BARD1', 'BLM', 'BRCA1',
#'             'BRCA2', 'BRIP1', 'BUB1B', 'CDC73', 'CHEK1',
#'              'FANCA', 'FANCD2', 'MDC1', 'MDM2', 'MRE11',
#'              'MSH6', 'NBN', 'PALB2', 'POLH', 'RAD50',
#'              'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'TOP2A',
#'              'TP53', 'WRN', 'XRCC2', 'XRCC3')
#'
#'result_tb1 <- RunDEGEmbedR(
#'degs = NeST79,
#'category = "GOBP"
#')
#'
#' head(res_tb1)
#'
#' # (B) Using AI-derived functional hypotheses
#' # NeST: 105 (Ubiquitin Regulation of p53 Activity)
#' NeST105 <- c('CUL3', 'ELOC', 'FBXW7', 'HSP90AA1', 'MDM2', 'SKP1', 'STK11', 'TNF', 'VHL')
#'
#' result_tb2 <- RunDEGEmbedR(
#'  degs = NeST105,
#'  category = "GSAI",
#'  api_key = api_key
#')
#'
#' head(res_tb2)
#'
#'
#' @seealso
#'   \code{\link[stats]{wilcox.test}},
#'   \code{\link[effsize]{cliff.delta}},
#'   \code{\link[tibble]{tibble}},
#'   \code{\link[stringr]{str_extract}}
#'
#' @importFrom stats wilcox.test
#' @importFrom stringr str_extract
#' @keywords enrichment embeddings similarity statistics effect-size
#' @export



RunDEGEmbedR <- function(degs,
                         bkgs = NULL,
                         category = c("GOBP","C2CP_all","BIOCARTA", "KEGG","PID","REACTOME", "WP", "GSAI"),
                         api_key = NULL,
                         output = TRUE){

  #Input checks
  if (missing(degs) || is.null(degs) || identical(degs, "")) stop("DEGs is required.")
  if (missing(category) || is.null(category) || identical(category, "")) stop("Category is required.")

  ###Load data###
  gsai_prompt <- readRDS(data_path("GSAIPrompt.rds"))
  bp_dt    <- readRDS(data_path("bp_function_retrieval_similarity.rds"))
  cp_dt    <- readRDS(data_path("cp_pathway_retrieval_similarity.rds"))
  rag_function <- readRDS(data_path("function_retrieval_2023_top40_embedding_2026020901.rds"))
  rag_pathway <- readRDS(data_path("pathway_retrieval_2023_top40_embedding_2026020901.rds"))

  ##category type##
  cat_upper <- toupper(category)

  if (cat_upper %in% c("GOBP", "GSAI")) {
    source_df <- rag_function
  } else if (cat_upper %in% c("C2CP_ALL", "BIOCARTA", "KEGG", "PID", "REACTOME", "WP")) {
    source_df <- rag_pathway
  } else {
    stop("Unknown category")
  }

  gene_dt   <- source_df[, -1]
  gene_list <- source_df$Name

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
    message(sprintf("No background gene list detected. Defaulting to %d built-in genes.", length(match_bkgs)))
  } else {
    # setdiff+intersect are C-optimized; avoids which()/sum()>0 branching
    match_bkgs <- intersect(setdiff(bkgs, degs), gene_list)
    message(sprintf("Background genes: There are %d matched genes with build-in gene list", length(match_bkgs)))
  }


  ##degs##
  match_degs <- intersect(degs, gene_list)
  if(length(match_degs) < 15 | length(match_degs) > 500){
    stop("Insufficient number of genes to match")
  }else{
    message(sprintf("Differentially expressed genes: There are %d matched genes with build-in gene list", length(match_degs)))
  }



  ###GSAI###

  if(toupper(category) == toupper("GSAI")){
  if (missing(api_key) || is.null(api_key) || identical(api_key, "")) stop("Missing API key.")

  ##Run GSAI##
    gsai <- RunGSAI(degs = match_degs,api_key = api_key,gsai_prompt = gsai_prompt,output = TRUE)

  ##Generate Text embedding##
    embed_mat <- GenerateTextEmbedding(
      text = gsai,
      api_key,
      output =FALSE
    )

    #Cosine similarity calculation
    # normalize rows
    gene_norm <- gene_dt / sqrt(rowSums(gene_dt * gene_dt))
    pth_norm  <- embed_mat / sqrt(rowSums(embed_mat * embed_mat))
    gene_norm <- as.matrix(gene_norm)
    mode(gene_norm) <- "numeric"

    pth_norm <- as.matrix(pth_norm)
    mode(pth_norm) <- "numeric"
    custom_matrix <- gene_norm %*% t(pth_norm)

  }


  ###Output table###
  results <- tibble::tibble(
    name                                      = character(),
    p_value_MWN_one_tailed                    = numeric(),
    median_cosine_similarity_degs             = numeric(),
    median_cosine_similarity_bkgs             = numeric(),
    diff_cosine_similarity                   = numeric(),
    cliffs_delta                              = numeric(),
    cliffs_delta_ci_95                        = character(),
    cliffs_delta_magnitude                    = character(),
    top10_degs_with_highest_cosine_similarity = character()
  )

  cp_path <- colnames(cp_dt)

  ###Select consine similarity matrix###
  con_sim_mtrx <- switch(toupper(category),
                         GOBP       = bp_dt,
                         C2CP_ALL   = cp_dt,
                         BIOCARTA   = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "BIOCARTA")],
                         KEGG       = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "KEGG")],
                         PID        = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "PID")],
                         REACTOME   = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "REACTOME")],
                         WP         = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "WP")],
                         GSAI       = custom_matrix,
                         stop("Unknown category")
  )

  tb <-  as.data.frame(con_sim_mtrx[rownames(con_sim_mtrx) %in% c(match_degs,match_bkgs), , drop = FALSE])
  for (i in 1:ncol(tb)) {


    cos_sim_degs <- tb[match_degs, colnames(con_sim_mtrx)[i]]
    names(cos_sim_degs) <- match_degs
    cos_sim_bkgs <- tb[!row.names(tb) %in% match_degs, colnames(con_sim_mtrx)[i]]
    ###calculate pvalue###
    if(i==1){
      message("Calculating p-values using Wilcoxon rank-sum test...")
    }
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
