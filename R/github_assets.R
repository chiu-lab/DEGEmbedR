.pkg_name   <- "DEGEmbedR"
.repo_owner <- "chiu-lab"
.repo_name  <- "DEGEmbedR"

# RDS lists
.required_rds_files <- c(
  "BP_15-500_similarity_text-embedding-3-large_2024110801.rds",
  "CP_15-500_similarity_text-embedding-3-large_2024110801.rds",
  "gene_geneset_moa_drug_function_similarity_text-embedding-3-large_2025031701.rds",
  "gene_embedding_2024110801.rds",
  "gene_list.rds"
)

# tool function
.cache_dir <- function() {
  dir <- getOption("DEGEmbedR.cache_dir", tools::R_user_dir(.pkg_name, "cache"))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

.version_string <- function() as.character(utils::packageVersion(.pkg_name))
.stamp_path    <- function() file.path(.cache_dir(), sprintf("assets-%s.stamp", .version_string()))

# Construct the release tag and asset URL based on the current package version
.asset_url <- function(filename) {
  ver <- .version_string()
  sprintf("https://github.com/%s/%s/releases/download/v%s/%s",
          .repo_owner, .repo_name, ver, filename)
}

# get_data_root
get_data_root <- function() {
  .ensure_assets_downloaded(.required_rds_files)
  .cache_dir()
}

# read RDS data
data_path <- function(filename) {
  .ensure_assets_downloaded(filename)
  file.path(.cache_dir(), filename)
}

# clear_data_cache
clear_data_cache <- function() {
  files <- list.files(.cache_dir(), full.names = TRUE, all.files = TRUE, no.. = TRUE)
  unlink(files, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

# if no data or stamp is not match, download
.ensure_assets_downloaded <- function(filenames, force = FALSE) {
  if (is.null(filenames)) filenames <- character()
  if (!length(filenames)) return(invisible(TRUE))

  need_version_refresh <- force || !file.exists(.stamp_path())
  exist_mask <- file.exists(file.path(.cache_dir(), filenames))
  need_files <- filenames[!exist_mask | need_version_refresh]

  if (!length(need_files)) return(invisible(TRUE))

  message(sprintf("Fetching data assets for %s %s ...", .pkg_name, .version_string()))
  for (fn in unique(need_files)) {
    url <- .asset_url(fn)
    dst <- file.path(.cache_dir(), fn)
    tmp <- tempfile(fileext = ".bin")
    ok  <- try(utils::download.file(url, tmp, mode = "wb", quiet = TRUE), silent = TRUE)
    if (inherits(ok, "try-error")) {
      stop(sprintf("Download failed：%s\nURL: %s", fn, url), call. = FALSE)
    }
    # remove failed files
    file.rename(tmp, dst)
  }

  # update stamp
  old_stamps <- Sys.glob(file.path(.cache_dir(), "assets-*.stamp"))
  if (length(old_stamps)) file.remove(old_stamps)
  file.create(.stamp_path())
  invisible(TRUE)
}
