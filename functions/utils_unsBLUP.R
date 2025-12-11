#' unsBLUP
#'
#' @description Function to process BLUP from sommer unstructured model
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd

unsBLUP <- function(blups){
  l <- unlist(lapply(blups,function(x){length(x[[1]])}))
  lmin <- min(l); lmax <- max(l)
  indexCov1 <- 1:lmin
  indexCov2 <- (lmin+1):lmax
  ntraits <- length(blups[[1]])
  # blups follow the order of a lower triangula matrix
  # (n*(n-1))/2 = l
  # l*2 = n2 - n
  # n2 = l*2 - n
  n <- 1:100
  possibilities <- ((n*(n-1))/2) + n
  ntrue <- n[which(possibilities == length(l))]
  ## index to know how to add them up
  base <- matrix(NA,ntrue,ntrue)
  base[lower.tri(base, diag=TRUE)] <- 1:length(l)
  index <- which(!is.na(base), arr.ind = TRUE)
  index <- index[order(index[,1]), ]


  for(i in 1:ntrue){ # for each main blup
    main <- which(index[,1] == i & index[,2] == i, arr.ind = TRUE)
    cov1 <- which(index[,1] == i & index[,2] != i, arr.ind = TRUE)
    cov2 <- which(index[,1] != i & index[,2] == i, arr.ind = TRUE)
    for(itrait in 1:ntraits){
      start <- blups[[main]][[itrait]]
      for(icov1 in cov1){
        start <- start + blups[[icov1]][[itrait]][indexCov1]
      }
      for(icov2 in cov2){
        start <- start + blups[[icov2]][[itrait]][indexCov2]
      }
      # store adjusted blup adding covariance effects in the same structure
      blups[[main]][[itrait]] <- start
    }
  }
  return(blups)
}
