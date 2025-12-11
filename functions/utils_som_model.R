#' som_model
#'
#' @description Calculate sommer models
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd

som_model <- function(y, Ng, Ne, K, type = "UN", nF = 1){

  GID <- rep(paste0("G", nm_0_lead(Ng)), Ne)
  EID <- rep(paste0("E", nm_0_lead(Ne)), each = Ng)

  rownames(K) <- colnames(K) <- paste0("G", nm_0_lead(Ng))

  d <- data.frame(GID, env = EID, y = c(y))
  d$GID <- factor(x = d$GID, levels = paste0("G", nm_0_lead(Ng)))
  d$env <- factor(x = d$env, levels = paste0("E", nm_0_lead(Ne)))

  if(type == "UN"){
    G_VCOV <- "~vsr(usr(env),GID,Gu=K)"
    R_VCOV <- "~vsr(dsr(env),units)"
  } else if (type == "DIAG"){
    G_VCOV <- "~vsr(dsr(env),GID,Gu=K)"
    R_VCOV <- "~vsr(dsr(env),units)"
  } else if(type == "ID"){
    G_VCOV = "~vsr(GID,Gu=K)"
    R_VCOV <- "~units"
  }


  if(type == "FA"){

    Ki <- as(solve(K + diag(1e-4,ncol(K),ncol(K))), Class="dgCMatrix")

    m <- mmec(y ~ env, random= ~vsc(usc(rrc(env, GID, y, nPC = 1)), isc(GID), Gu = Ki),
              rcov= ~vsc(dsc(env),isc(units)),
              data=d, verbose = FALSE, dateWarning = FALSE)


  } else {

    m <- tryCatch(mmer(y ~ env,
                       random= as.formula(G_VCOV),
                       rcov= as.formula(R_VCOV),
                       data=d, verbose = FALSE, dateWarning = FALSE),
                  error = function(x) NULL)

  }

  if(!is.null(m)){

    return(extract_som_estimates(m = m, type = type, nF = nF, Ne = Ne, Ng = Ng))

  } else {

    return(list(B = NA, R_hat = NA, SG_hat = NA, y_hat = y_hat, g_hat = NA))

  }

}
