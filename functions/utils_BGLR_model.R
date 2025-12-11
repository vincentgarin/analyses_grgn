#' BGLR_model
#'
#' @description Calculate BGLR model
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd

BGLR_model <- function(y, Ng, Ne, K, type = "UN", nF = 1, nIter = 10000,
                       burnIn = 1000, thin = 10){

  if(type == "UN"){
    G_VCOV <- list(list(K = K, model="RKHS",
                        Cov = list(type = "UN", df0 = Ne)))
  } else if (type == "DIAG") {
    G_VCOV <- list(list(K = K, model="RKHS",
                        Cov = list(type = "DIAG", S0 = rep(1, Ne), df0 = rep(1, Ne))))

  } else if (type == "FA"){
    M <- matrix(nrow = Ne, ncol = nF, TRUE)
    G_VCOV <- list(list(K = K, model="RKHS",
                        Cov = list(type = "FA", M = M,
                                   S0 = rep(1, Ne),
                                   df0 = rep(1, Ne), var = 100)))

  }

  m <- tryCatch(suppressMessages(Multitrait(y = y,
                           ETA = G_VCOV,
                           nIter = nIter, burnIn = burnIn, thin = thin,
                           resCov = list(type = 'DIAG', S0 = rep(1, Ne),
                                         df0 = rep(1, Ne)),
                           verbose = FALSE, saveAt = file.path(tempdir(), "sim"))),
                error = function(x) NULL)

  if(!is.null(m)){

    # fixed effects
    B <- m$mu

    # VCOV components - error
    R_hat <- m$resCov$R

    # VCOV components - GxE
    SG_hat <- m$ETA[[1]]$Cov$Omega

    # predicted phenotype
    y_hat <- data.frame(m$ETAHat)
    colnames(y_hat) <- paste0('E', nm_0_lead(Ne))
    rownames(y_hat) <- paste0("G", nm_0_lead(Ng))

    # genetic effects
    BLUP <- data.frame(m$ETA[[1]]$u)
    colnames(BLUP) <- paste0('E', nm_0_lead(Ne))
    rownames(BLUP) <- paste0("G", nm_0_lead(Ng))

    Nobs <- c(y)
    Nobs <- length(Nobs[!is.na(Nobs)])

    return(list(B = B, R_hat = R_hat, SG_hat = SG_hat, y_hat = y_hat, g_hat = BLUP))

  } else {

    return(list(B = NA, R_hat = NA, SG_hat = NA, y_hat = NA, g_hat = NA))

  }


}
