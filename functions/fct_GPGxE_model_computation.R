#' GPGxE_model_computation
#'
#' @description Calculate different GP GxE models
#'
#' @details
#' The function can calculate different GP GxE models ...
#'
#' @param model Character string indicating the type of model you want to
#' calculate. Must be one of : "BGGE_MM", "BGGE_MDe", "som_ID", "som_DIAG",
#' "som_UN", "som_FA", "BGLR_DIAG", "BGLR_UN", "BGLR_FA". Default = "BGGE_MM"
#'
#' @param nF number of parameters for the factor analytic model (FA(nF)). Default = 1
#'
#' @param pheno Numeric matrix of phenotypic data with N_ind x N_env. Validation set data
#' can be masked with NA.
#'
#' @param K Numeric relationship matrix. Later possibility to add a list of kernels
#'
#' @param nIter Numeric value specifying the number of iteration for Bayesian models.
#' Default = 10000
#'
#' @param burnIn Numeric value specifying the number of iteration that are discared
#' for Bayesian models. Default = 1000
#'
#' @param thin thin: Numeric value specifying thinning for Bayesian models
#' (posterior estimations are based on values sampled at 'thin' interval).
#' Default = 10
#'
#' @return list containing: the model estimates
#' (B: fixed effects intercepts; R_hat: error VCOV; SG_hat: genotypic VCOV;
#' g_hat: random genotypic effect estimates or BLUP) and the computational time
#'
#' @examples
#'
#' library(MTM)
#'
#' data("K", package = "GPGxE")
#' data("pheno", package = "GPGxE")
#'
#' m1 <- GPGxE_model_computation(model = "BGGE_MM", pheno = pheno, K = K)
#' # m2 <- GPGxE_model_computation(model = "BGGE_MDe", pheno = pheno, K = K)
#' # m3 <- GPGxE_model_computation(model = "som_ID", pheno = pheno, K = K)
#' # m4 <- GPGxE_model_computation(model = "som_DIAG", pheno = pheno, K = K)
#' # m5 <- GPGxE_model_computation(model = "som_UN", pheno = pheno, K = K)
#' # m6 <- GPGxE_model_computation(model = "som_FA", pheno = pheno, K = K)
#' # m7 <- GPGxE_model_computation(model = "BGLR_DIAG", pheno = pheno, K = K)
#' # m8 <- GPGxE_model_computation(model = "BGLR_UN", pheno = pheno, K = K)
#' # m9 <- GPGxE_model_computation(model = "BGLR_FA", pheno = pheno, K = K)
#'
#' @import BGGE
#' @import BGLR
#' @import dplyr
#' @import ggplot2
#' @import reshape2
#' @import sommer
#' @import stringr
#' @importFrom methods as
#' @importFrom stats as.formula cor model.matrix rnorm
#'
#' @export

GPGxE_model_computation <- function(model = "BGGE_MM", nF = 1, pheno, K,
                                    nIter = 10000, burnIn = 1000, thin = 10){

  # check arguments ----

  model_list <- c("BGGE_MM", "BGGE_MDe", "som_ID", "som_DIAG", "som_UN", "som_FA",
                  "BGLR_DIAG", "BGLR_UN", "BGLR_FA")
  if (!(model %in% model_list)){
    stop(paste0("model must be one of: ", paste(model_list, collapse = ", ")))
  }

  if(!is.matrix(pheno)){
    stop("pheno must be a matrix.")

  } else{

    if(!is.numeric(pheno)){
      stop("pheno must be numeric.")
    }

  }

  if(!is.matrix(K)){
    stop("K must be a matrix.")

  } else{

    if(!is.numeric(K)){
      stop("K must be numeric.")
    }

  }

  # collect arguments ----
  Ng <- nrow(pheno)
  Ne <- ncol(pheno)

  models_id <- c("BGGE_MM", "BGGE_MDe", "som_ID", "som_DIAG", "som_UN",
                 "som_FA", "BGLR_DIAG", "BGLR_UN", "BGLR_FA")

  pkg_type <- unlist(strsplit(x = model, split = "_"))
  pkg <- pkg_type[1]
  type <- pkg_type[2]

  # model calculation ----

  t1 <- Sys.time()

  if(pkg == "som"){

    m <- som_model(y = pheno, Ng = Ng, Ne = Ne, K = K, type = type, nF = nF)

  } else if (pkg == "BGLR"){

    m <- BGLR_model(y = pheno, Ng = Ng, Ne = Ne, K = K, type = type, nF = nF,
                   nIter = nIter, burnIn = burnIn, thin = thin)

  } else if (pkg == "BGGE"){

    m <- BGGE_model(y = pheno, Ng = Ng, Ne = Ne, K = K, type = type,
                   nIter = nIter, burnIn = burnIn, thin = thin)

  }

  t2 <- Sys.time()


  # results ----

  d_time <- as.numeric(difftime(t2, t1, units = "secs"))

  # return(list(m = m, time = d_time))
  return(m)

}
