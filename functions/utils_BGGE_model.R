#' BGGE_model
#'
#' @description Calculate BGGE models
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd

# model = models[m]
# y = Yij
# Ng = Ng
# Ne = Ne
# K = K
# type = "MDe"
# nIter = nIter
# burnIn = burnIn
# thin = thin

# ajout option pour random genotype
# option pour Kernel type (GK vs GB)
# Ajout plusieurs matrices

BGGE_model <- function(y, Ng, Ne, K, type = "MM",
                       nIter = 10000, burnIn = 1000, thin = 10){

  GID <- rep(paste0("G", nm_0_lead(Ng)), Ne)
  EID <- rep(paste0("E", nm_0_lead(Ne)), each = Ng)

  rownames(K) <- colnames(K) <- paste0("G", nm_0_lead(Ng))

  Y <- data.frame(env = EID, GID, y = c(y))
  Y$GID <- factor(x = Y$GID, levels = paste0("G", nm_0_lead(Ng)))
  Y$env <- factor(x = Y$env, levels = paste0("E", nm_0_lead(Ne)))

  Ker <- getK(Y = Y, setKernel = list(K), model = type)
  XE <- model.matrix(~ -1 + env, data = Y)

  m <- tryCatch(BGGE(y = Y[, 3], K = Ker, XF = XE, ne = rep(Ng, Ne),
                     ite = nIter, burn = burnIn, thin = thin),
                error = function(x) NULL)

  if(!is.null(m)){

    return(extract_BGGE_estimates(m = m, type = type, Ne = Ne, Ng = Ng))

  } else {

    return(list(B = NA, R_hat = NA, SG_hat = NA, y_hat = NA, g_hat = NA))

  }

}
