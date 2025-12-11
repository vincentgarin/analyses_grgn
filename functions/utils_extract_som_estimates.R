#' extract_som_estimates
#'
#' @description function to extract the model estimates from sommer
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd

extract_som_estimates <- function(m, type = "UN", nF, Ne, Ng){

  if(type == "UN"){

    # fixed effects
    B <- m$Beta$Estimate
    B[2:Ne] <- B[2:Ne] + B[1]

    # VCOV components - error
    S_vec <- m$sigmaVector
    S_ge_term <- ((Ne*(Ne+1))/2)
    Se_hat <- S_vec[(S_ge_term + 1):length(S_vec)]
    R_hat <- diag(Se_hat)

    # VCOV components - GxE
    SG_hat <- matrix(NA, Ne, Ne)
    SG_hat[upper.tri(SG_hat, diag = TRUE)] <- S_vec[1:S_ge_term]
    SG_hat[lower.tri(SG_hat)] <- t(SG_hat)[lower.tri(t(SG_hat))]

    # genetic effects
    BLUP <- unsBLUP(blups = m$U)

    seq_v <- (1:Ne * (1:Ne + 1))/2

    BLUP <- matrix(unlist(BLUP[seq_v]), nrow = Ng)
    BLUP <- data.frame(BLUP)

  } else if (type == "DIAG"){

    # fixed effects
    B <- m$Beta$Estimate
    B[2:Ne] <- B[2:Ne] + B[1]

    # VCOV components - error
    S_vec <- m$sigmaVector
    Se_hat <- S_vec[(Ne+1):length(S_vec)]
    R_hat <- diag(Se_hat)

    # VCOV components - GxE
    SG_hat <- diag(S_vec[1:Ne])

    # genetic effects
    BLUP <-m$U
    BLUP <- matrix(unlist(BLUP), nrow = Ng)
    BLUP <- data.frame(BLUP)

  } else if (type == "ID"){

    # fixed effects
    B <- m$Beta$Estimate
    B[2:Ne] <- B[2:Ne] + B[1]

    # VCOV components - error
    S_vec <- m$sigmaVector
    R_hat <- diag(Ne) * S_vec[2]

    # VCOV components - GxE
    SG_hat <- diag(Ne) * S_vec[1]

    # genetic effects
    BLUP <-m$U
    BLUP <- matrix(unlist(BLUP), nrow = Ng)
    BLUP <- data.frame(matrix(1, nrow = 1, Ne) %x% BLUP)

  } else if (type == "FA"){

    # fixed effects
    B <- m$b[, 1]
    B[2:Ne] <- B[2:Ne] + B[1]

    # VCOV components - error
    S_vec <- m$sigma
    S_vec <- rev(rev(S_vec)[1:Ne])
    R_hat <- diag(S_vec)

    # VCOV components - GxE
    SG_hat <- matrix(NA, Ne, Ne)

    # genetic effects
    Gamma <- with(m$data, rrc(env, GID, y, returnGamma = TRUE, nPC = nF))$Gamma
    score.mat <- m$uList[[1]] # extract factor scores
    BLUP <- score.mat %*% t(Gamma) # BLUPs in all environments

  }

  XB <- matrix(rep(B, each = Ng), nrow = Ng)
  y_hat <- XB + BLUP

  colnames(BLUP) <- colnames(y_hat) <- paste0("E", nm_0_lead(Ne))
  rownames(BLUP) <- rownames(y_hat) <- paste0("G", nm_0_lead(Ng))

  return(list(B = B, R_hat = R_hat, SG_hat = SG_hat, y_hat = y_hat, g_hat = BLUP))

}
