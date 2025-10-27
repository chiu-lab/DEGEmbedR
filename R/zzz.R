.onLoad <- function(libname, pkgname) {
  # ~/.cache/DEGEmbedR/（Linux）| ~/Library/Caches/DEGEmbedR/（macOS）| AppData\Local\…（Windows）
  options(DEGEmbedR.cache_dir = tools::R_user_dir(pkgname, "cache"))
}
