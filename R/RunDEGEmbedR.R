#' Compare DEG vs. background cosine similarities across functions, pathways, or MOAs
#'
#' @description
#' For each built-in biological function or user-supplied term embedding, this function compares the
#' cosine similarity distributions of differentially expressed genes (DEGs) and background genes,
#' summarizing the statistical evidence for a functional relationship.
#'
#' The function returns the one-tailed Wilcoxon rank-sum test p-value (testing whether DEGs >
#' background), the median cosine similarities for DEGs and background genes, their median
#' difference, Cliff’s delta with a 95% confidence interval, and the top 10 DEGs with the highest
#' similarity to the function.
#'
#' @details
#' This function implements the analysis workflow described in the DEGEmbedR vignette and
#' manuscript (e.g., enrichment testing for the STING pathway or drug mechanism hypotheses). It
#' supports curated functional collections included with the package (MSigDB GO Biological
#' Processes; MSigDB C2 pathways including BIOCARTA, KEGG, PID, REACTOME, and WP; and mechanisms
#' of action [MOA]), as well as a \code{"Customized"} mode that accepts user-supplied pathway or
#' function embeddings.
#'
#' @param degs Character vector. Differentially expressed genes (DEGs). After intersecting with
#'   the built-in gene universe (~18,000 genes), the number of matched DEGs must be between 15 and 500.
#' @param bkgs Character vector. Background genes (optional). If \code{NULL}, the built-in gene
#'   universe—comprising ~18k human protein-coding genes from the NCBI Gene database—is used
#'   as the background.
#' @param category Character. Functional database or mode to use. Must be one of:
#'   \code{"GOBP"}, \code{"C2CP_all"}, \code{"BIOCARTA"}, \code{"KEGG"}, \code{"PID"},
#'   \code{"REACTOME"}, \code{"WP"}, \code{"MOA"}, or \code{"Customized"} (case-insensitive).
#' @param embedding_input Numeric matrix or data frame of term embeddings (only required if
#'   \code{category = "Customized"}). Each row is a term, each column an embedding dimension
#'   (length = 3072). Row names are used as term labels.
#' @param output Logical. Whether to save the result as a timestamped tab-delimited \code{.txt} file.
#'   Default: \code{TRUE}.
#'
#' @return A \link[tibble]{tibble} with one row per pathway or term, including:
#' \describe{
#'   \item{\code{name}}{Pathway or term name}
#'   \item{\code{p_value_MWN_one_tailed}}{One-tailed Wilcoxon p-value (DEGs > background)}
#'   \item{\code{median_cosine_similarity_degs}}{Median cosine similarity among DEGs}
#'   \item{\code{median_cosine_similarity_bkgs}}{Median cosine similarity among background genes}
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
#' # Example using user-supplied embeddings
#' load(system.file("examples", "example.rdata", package = "DEGEmbedR"))
#' head(embed_mat)
#' length(degs)
#' res <- RunDEGEmbedR(
#'   degs = degs,
#'   category = "Customized",
#'   embedding_input = embed_mat
#' )
#' head(res)
#' }
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
                         C2CP_all   = cp_dt,
                         BIOCARTA   = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "BIOCARTA")],
                         KEGG       = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "KEGG")],
                         PID        = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "PID")],
                         REACTOME   = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "REACTOME")],
                         WP         = cp_dt[,which(str_extract(cp_path, "^[^_]+") == "WP")],
                         MOA        = moa_gene_dt,
                         CUSTOMIZED = custom_matrix,
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
