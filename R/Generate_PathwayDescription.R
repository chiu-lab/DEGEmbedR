#' Generate a pathway/MOA description using OpenAI API
#'
#' @description
#' Queries the OpenAI Chat Completions API to generate a short descriptive paragraph summarizing the function of a given pathway.
#'
#' @param pathway Character. The name of the pathway/ MOA.
#' @param api_key Character. OpenAI API key.
#' @param api_url Character. OpenAI API URL. Default: \code{"https://api.openai.com/v1/chat/completions"}.
#' @param request_timeout_sec Numeric. Request timeout in seconds. Default: 60.
#' @param output Logical. The result will be saved as a .txt file (Default: TRUE).
#'
#' @return A character vector of length 1 containing the generated pathway description.
#' @details
#' The description is printed, returned, and saved as a timestamped text file.
#' @examples
#' \dontrun{
#' desc <- Generate_Description(
#'   pathway = "Gluconeogenesis",
#'   api_key = Sys.getenv("OPENAI_API_KEY")
#' )
#' print(desc)
#' }
#' @seealso \code{\link[httr]{POST}}, \code{\link[jsonlite]{fromJSON}}
#' @export


Generate_PathwayDescription <- function(
    pathway,
    api_key,
    api_url = "https://api.openai.com/v1/chat/completions",
    request_timeout_sec = 60,
    output = TRUE
) {
  if (!requireNamespace("httr", quietly = TRUE)) stop("Package 'httr' is required.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required.")

  #Input checks
  if (missing(api_key) || is.null(api_key) || identical(api_key, "")) stop("Missing API key.")
  if (missing(pathway)  || is.null(pathway)  || identical(pathway,  "")) stop("Please provide pathway name.")

  #Prompt
  prompt <- sprintf("Write a paragraph to describe what is known about the function of '%s'.", pathway)
  message("Prompt: ", prompt)
  message("Start Processing!")

  #Request body
  body <- list(
    model = "gpt-4o-2024-08-06",
    messages = list(
      list(role = "system", content = "You are a biomedical knowledge expert"),
      list(role = "user",   content = prompt)
    ),
    temperature = 1
  )

  #HTTP request
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

  #Parse & surface errors
  parsed <- httr::content(resp, "parsed", encoding = "UTF-8")

  # If HTTP error, try to show API message
  if (httr::http_error(resp)) {
    msg <- tryCatch(parsed$error$message, error = function(e) NULL)
    if (is.null(msg)) msg <- paste("HTTP", httr::status_code(resp), httr::http_status(resp)$message)
    stop(sprintf("API error: %s", msg))
  }

  # Extract text
  out <- parsed$choices[[1]]$message$content
  if (is.null(out)) stop("No content returned from API.")
  out <- trimws(out)

  # Name it with the pathway so downstream rownames work
  names(out) <- pathway
  cat("Description:\n", out, "\n")
if(output){
  cat( pathway, "\n\n", out, file = paste("description", format(Sys.time(), "%Y-%m-%d-%H%M%S.txt"),sep = "_"),fileEncoding = "UTF-8")
}
  return(out)
}
