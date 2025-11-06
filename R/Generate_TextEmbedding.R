#' Generate numerical embeddings for input text using the OpenAI API
#'
#' @description
#' Sends input text to the OpenAI Embeddings API endpoint and returns a numeric
#' matrix of text embeddings. Each row represents one text input, and each column
#' represents one embedding dimension (length = 3072).
#'
#' @param text Character vector of text strings to embed (pathway/MOA description from Generate_PathwayDescription() function).
#' @param api_key Character string. Your OpenAI API key.
#' @param api_url Character string. The API endpoint for embeddings.
#'   Defaults to \code{"https://api.openai.com/v1/embeddings"}.
#' @param request_timeout_sec Numeric. Timeout for the API request, in seconds.
#'   Defaults to \code{60}.
#' @param output Logical. The result will be saved as a .txt file (Default: TRUE).
#'
#' @details
#' This function interfaces with the OpenAI model \code{text-embedding-3-large},
#' which converts natural language text into dense vector representations.
#'
#' The output matrix has rows corresponding to the input texts and columns
#' corresponding to embedding dimensions (3,072 for this model).
#'
#' @return
#' A numeric matrix where each row corresponds to a text input and each column
#' corresponds to an embedding dimension. Row names are derived from
#' \code{names(text)} if provided, otherwise truncated text strings.
#'
#'
#' @examples
#' \dontrun{
#' # Generate embeddings
#' load(system.file("examples","example.rdata",package = "DEGEmbedR"))
#' txt <- c(desc1, desc2)
#' emb <- Generate_Embedding(
#'   text = txt,
#'   api_key = Sys.getenv("OPENAI_API_KEY")
#' )
#'
#' dim(emb)        # rows = 2, cols = 3072
#' head(emb[, 1:5]) # first few dimensions
#' }
#'
#' @seealso
#' \code{\link[httr]{POST}}, \code{\link[jsonlite]{fromJSON}}
#'
#' @references
#' OpenAI API documentation:
#' \url{https://platform.openai.com/docs/api-reference/embeddings}
#'
#' @keywords OpenAI embedding text-representation bioinformatics
#' @export

Generate_TextEmbedding <- function(
    text,
    api_key,
    api_url = "https://api.openai.com/v1/embeddings",
    request_timeout_sec = 60,
    output =TRUE
) {
  if (!requireNamespace("httr", quietly = TRUE)) stop("Package 'httr' is required.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required.")

  #Input checks
  if (missing(api_key) || is.null(api_key) || identical(api_key, "")) {
    stop("Missing API key.")
  }
  if (missing(text) || length(text) == 0) {
    stop("Missing text.")
  }
  if (!is.character(text)) {
    stop("'text' must be a character vector.")
  }

  # Ensure plain character vector without names
  input_vec <- unname(as.character(text))

  #Build body
  body <- list(
    model = "text-embedding-3-large",
    input = input_vec
  )

  #Request
  resp <- httr::POST(
    url = api_url,
    httr::add_headers(
      "Authorization" = paste("Bearer", api_key),
      "Content-Type"  = "application/json"
    ),
    body   = body,
    encode = "json",
    httr::timeout(request_timeout_sec)
  )

  #Parse / error handle
  parsed <- httr::content(resp, "parsed", encoding = "UTF-8")
  if (!is.null(parsed$error)) {
    stop(sprintf("API error: %s", parsed$error$message))
  }
  if (is.null(parsed$data) || length(parsed$data) == 0) {
    stop("No embeddings returned.")
  }

  #Build matrix: rows = inputs, cols = dims (length of embedding)
  emb_list <- lapply(parsed$data, function(x) as.numeric(x$embedding))
  # Validate equal lengths
  dims <- unique(vapply(emb_list, length, integer(1)))
  if (length(dims) != 1) stop("Inconsistent embedding sizes in response.")
  embed_mat <- do.call(rbind, emb_list)
  storage.mode(embed_mat) <- "double"

  #Row names: prefer names(text), else truncated text
  rn <- names(text)
  if (is.null(rn) || any(is.na(rn)) || any(rn == "")) {
    trunc_label <- function(s, n = 40) {
      s <- gsub("\\s+", " ", s)
      ifelse(nchar(s) <= n, s, paste0(substr(s, 1, n - 1), "…"))
    }
    rn <- vapply(as.character(text), trunc_label, character(1))
  }
  rownames(embed_mat) <- rn
  if(output){
  write.table(embed_mat,file = paste("embedding", format(Sys.time(), "%Y-%m-%d-%H%M%S.txt"),sep = "_"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  }
  return(embed_mat)
}
