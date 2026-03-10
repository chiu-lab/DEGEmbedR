# Generate a function to run GASI

RunGSAI <- function(
    degs,
    api_key,
    gsai_prompt,
    api_url = "https://api.openai.com/v1/chat/completions",
    request_timeout_sec = 60,
    model = "gpt-4o-2024-08-06",
    temperature = 1,
    output = FALSE,
    max_retries = 1,
    min_rows = 1
) {

  # Dependencies
  if (!requireNamespace("httr", quietly = TRUE)) stop("Package 'httr' is required.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package 'tibble' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")

  # Check inputs
  if (missing(api_key) || is.null(api_key) || identical(api_key, "")) stop("Missing API key.")
  if (missing(gsai_prompt) || is.null(gsai_prompt) || identical(gsai_prompt, "")) stop("Missing gsai_prompt.")
  if (missing(degs) || is.null(degs) || length(degs) == 0) stop("Please provide 'degs' (non-empty).")

  degs <- as.character(degs)

  # Prompt
  prompt <- paste(gsai_prompt, paste(degs, collapse = "\n"), sep = "\n")

  call_llm <- function(prompt_text) {
    body <- list(
      model = model,
      messages = list(
        list(role = "system", content = "You are an efficient and insightful assistant to a molecular biologist"),
        list(role = "user",   content = prompt_text)
      ),
      temperature = temperature
    )

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

    parsed <- httr::content(resp, "parsed", encoding = "UTF-8")

    if (httr::http_error(resp)) {
      msg <- tryCatch(parsed$error$message, error = function(e) NULL)
      if (is.null(msg)) msg <- paste("HTTP", httr::status_code(resp), httr::http_status(resp)$message)
      stop(sprintf("API error: %s", msg))
    }

    out <- parsed$choices[[1]]$message$content
    if (is.null(out) || identical(trimws(out), "")) stop("No content returned from API.")
    trimws(out)
  }

  # Helper: detect missing / null-like text
  is_missing_text <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    is.na(x) |
      x == "" |
      tolower(x) %in% c("na", "n/a", "null", "none")
  }

  parse_out <- function(out_text) {
    # split paragraphs
    parts <- stringr::str_split(out_text, "\n\\s*\n")[[1]]
    parts <- parts[trimws(parts) != ""]

    sep_class <- "[,:;\\-–—]"

    df <- tibble::tibble(text = parts) |>
      dplyr::mutate(
        confidence = suppressWarnings(
          as.numeric(stringr::str_extract(text, "(?<=confidence: )\\d+\\.?\\d*"))
        ),
        clean_text = stringr::str_remove(
          text,
          "^\\(LLM self-assessed confidence: [0-9.]+\\)\\s*"
        ),
        m = stringr::str_match(
          clean_text,
          paste0("^\\s*(.+?)\\s*", sep_class, "\\s*(.+)\\s*$")
        ),
        title = m[, 2],
        description = m[, 3]
      ) |>
      dplyr::select(-m) |>
      dplyr::mutate(
        title = stringr::str_squish(title),
        description = stringr::str_squish(description)
      ) |>
      dplyr::filter(
        !is_missing_text(title),
        !is_missing_text(description)
      ) |>
      dplyr::distinct(title, .keep_all = TRUE) |>
      dplyr::select(title, description, confidence)

    df
  }

  # Validate parsed output
  has_valid_result <- function(df, min_rows = 1) {
    if (is.null(df) || !is.data.frame(df)) return(FALSE)
    if (nrow(df) < min_rows) return(FALSE)
    if (!all(c("title", "description") %in% colnames(df))) return(FALSE)

    # reject if any title/description is NA, NULL-like, or empty
    if (any(is_missing_text(df$title))) return(FALSE)
    if (any(is_missing_text(df$description))) return(FALSE)

    TRUE
  }

  attempt <- 0
  last_out <- NULL
  last_df <- NULL

  while (attempt <= max_retries) {
    attempt <- attempt + 1

    prompt_try <- if (attempt == 1) {
      prompt
    } else {
      paste(
        prompt,
        "\n\nFORMAT STRICTLY as multiple paragraphs.",
        "\nEach paragraph MUST be exactly:",
        "\n(LLM self-assessed confidence: 0.xx) <TITLE> : <DESCRIPTION>",
        "\nRules:",
        "\n- Use exactly one ':' between title and description.",
        "\n- Do not use bullets or numbering.",
        "\n- TITLE cannot be NA, NULL, N/A, None, or empty.",
        "\n- DESCRIPTION cannot be NA, NULL, N/A, None, or empty.",
        "\n- If unsure, still provide the best valid pathway description.",
        sep = ""
      )
    }

    last_out <- tryCatch(
      call_llm(prompt_try),
      error = function(e) NULL
    )

    if (!is.null(last_out)) {
      last_df <- tryCatch(
        parse_out(last_out),
        error = function(e) NULL
      )
    } else {
      last_df <- NULL
    }

    if (has_valid_result(last_df, min_rows = min_rows)) {
      break
    }

    if (attempt <= max_retries) {
      message(sprintf(
        "Invalid or incomplete result (attempt %d/%d). Retrying LLM...",
        attempt, max_retries + 1
      ))
    }
  }

  # Final fallback
  if (!has_valid_result(last_df, min_rows = min_rows)) {
    final <- "No pathway avalibale"
    return(stats::setNames(final, final))
  }

  # Optional file output as a table
  if (isTRUE(output)) {
    fn <- paste0("gsai_table_", format(Sys.time(), "%Y-%m-%d-%H%M%S"), ".csv")
    utils::write.csv(last_df, fn, row.names = FALSE)
    message("Saved table to: ", fn)
  }

  # Return named vector
  final <- last_df$description
  names(final) <- last_df$title
  return(final)
}
