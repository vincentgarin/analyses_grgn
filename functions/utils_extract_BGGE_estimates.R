#' extract_BGGE_estimates
#'
#' @description Function to extract BGGE model estimates
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd

extract_BGGE_estimates <- function(m, type = "MM", Ne, Ng){

  # get the BLUP
  if (type == "MM"){

    BLUP <- matrix(m$K$G$u, nrow = Ng)

  } else if (type == "MDe"){

    BLUP_G <- matrix(m$K$G$u, nrow = Ng)

    BLUP_GE <- lapply(X = m$K[2:(Ne + 1)], `[[`, 1)
    BLUP_GE <- Reduce(`+`, BLUP_GE)
    BLUP_GE <- matrix(BLUP_GE, nrow = Ng)
    BLUP <- BLUP_G + BLUP_GE

  }

  BLUP <- data.frame(BLUP)
  colnames(BLUP) <- paste0("E", nm_0_lead(Ne))
  rownames(BLUP) <- paste0("G", nm_0_lead(Ng))

  # rest is the same for both models

  # variance components
  Se <- m$varE
  Sg <- m$K$G$varu
  # for now do not store those values specific to MDe ...
  # Sge <- unlist(lapply(X = m$K[2:(Ne + 1)], `[[`, 3))

  R_hat <- diag(rep(Se, Ne))
  SG_hat <- diag(rep(Sg, Ne))

  # get the fixed effects by difference yhat BLUP
  y_hat <- matrix(m$yHat, ncol = Ne)
  B <- y_hat - BLUP
  B <- colMeans(B)

  return(list(B = B, R_hat = R_hat, SG_hat = SG_hat, y_hat = y_hat, g_hat = BLUP))

}
