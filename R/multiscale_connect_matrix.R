#' Build a multiscale H3 connectivity matrix
#'
#' @description
#' Construct a symmetric sparse connectivity matrix linking each finer H3
#' planning unit to its parent at the next coarser resolution. This is
#' typically used as the `data` input to [add_multiscale_connectivity_penalties()].
#'
#' @param maps A list as returned by [build_h3_maps()], containing at least
#'   `h3_vec`, `res_vec`, `res_levels`, and `row_idx_by_h3`.
#' @param symmetric Logical; if `TRUE` (default), the matrix is symmetrised
#'   so that each parent–child link appears in both directions.
#'
#' @return A `Matrix::dgCMatrix` connectivity matrix of size `n_pu × n_pu`.
#' 
#' @examples
#' # Minimal 2-resolution parent-child example
#' h3_child  <- "8a2a1072b59ffff"
#' h3_parent <- "872a1072bffffff"
#'
#' maps <- build_h3_maps(
#'   s_or_h3 = c(h3_parent, h3_child),
#'   res_vec = c(7L, 8L)
#' )
#'
#' conn <- build_multiscale_connectivity_matrix(maps)
#' conn
#' 
#' @export
build_multiscale_connectivity_matrix <- function(maps, symmetric = TRUE) {
  h3_vec   <- maps$h3_vec
  res_vec  <- maps$res_vec
  res_lvls <- sort(unique(as.integer(maps$res_levels)))
  row_idx  <- maps$row_idx_by_h3
  
  n <- length(h3_vec)
  if (length(res_lvls) < 2L) {
    stop("Need at least two H3 resolutions to build a multiscale connectivity matrix.")
  }
  
  ii <- integer(0)
  jj <- integer(0)
  
  # link each finer resolution to its parent at the next coarser resolution
  for (k in seq_len(length(res_lvls) - 1L)) {
    r_low  <- res_lvls[k]
    r_high <- res_lvls[k + 1L]
    
    idx_high <- which(res_vec == r_high)
    if (!length(idx_high)) next
    
    parents <- h3jsr::get_parent(h3_address = h3_vec[idx_high], res = r_low, simple = TRUE)#parents <- h3forr::h3_to_parent(h3_vec[idx_high], r_low)
    parent_idx <- as.integer(row_idx[parents])
    
    ok <- which(!is.na(parent_idx))
    if (!length(ok)) next
    
    child_idx  <- idx_high[ok]
    parent_idx <- parent_idx[ok]
    
    ii <- c(ii, child_idx)
    jj <- c(jj, parent_idx)
  }
  
  if (!length(ii)) {
    return(Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE))
  }
  
  conn <- Matrix::sparseMatrix(
    i    = ii,
    j    = jj,
    x    = 1,
    dims = c(n, n)
  )
  
  if (symmetric) {
    conn <- conn + Matrix::t(conn)
  }
  
  conn
}
