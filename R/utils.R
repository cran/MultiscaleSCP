#' Internal utilities for MultiscaleSCP
#'
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

safe_len <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) 0L else length(x)
}

normalize <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

force_epsg4326 <- function(r) {
  cr <- terra::crs(r, proj = TRUE)
  if (is.na(cr) || !grepl("EPSG:4326", cr)) terra::project(r, "EPSG:4326") else r
}
