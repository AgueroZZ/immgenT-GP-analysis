filter_cells_by_total_membership <- function (L, max_val = 10, numiter = 10) {
  n <- nrow(L)
  rows <- 1:n
  for (iter in 1:numiter) {
    x <- rowSums(L)
    cat(sprintf("%d. Filtered out %d cells.\n",iter,sum(x > max_val)))
    i <- which(x <= max_val)
    L <- L[i,]
    rows <- rows[i]
    d <- apply(L,2,max)
    L <- scale_cols(L,1/d)
  }
  return(rows)
}

scale_cols <- function (A, b)
  t(t(A) * b)
