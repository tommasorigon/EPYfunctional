# library(mvtnorm)

ldet <- function(x) {
  if (!is.matrix(x)) {
    return(log(x))
  }
  determinant(x, logarithm = TRUE)$modulus
}

EFDMP <- function(X, Hl, L, B, indexB, time = NULL, cc,
                  alpha, Sigma, a_sigma, b_sigma,
                  maxiter = 1000, tol = 1e-7, verbose = FALSE) {

  # Internal checks to verify that the prior is specified correctly.
  if (length(Hl) != L) stop("The number of functional classes must coincide with the length of H")
  if (length(alpha) != L) stop("The number of functional classes must coincide with the length of alpha")
  if (length(cc) != L) stop("The number of functional classes must coincide with the length of cc")
  if (ncol(indexB) != ncol(B)) stop("The column dimension of B must coincide with the column dimension of indexB")
  if (nrow(indexB) != L) stop("The row dimension of indexB must coincide with L")
  if (length(Sigma) != L) stop("The dimension of Sigma must coincide with L")

  # Fixed quantities
  n <- NROW(X)
  TT <- NCOL(X)
  H <- sum(Hl)
  M <- ncol(B)
  Ml <- rowSums(indexB)
  Findex <- rep(1:L, Hl)

  # Further consistency checks
  if (is.null(time)) {
    time <- 1:TT
    warning("The time dimension was not supplied and equally-spaced intervals are assumed.")
  }

  for (l in 1:L) {
    if (Ml[l] != ncol(Sigma[[l]]) | Ml[l] != nrow(Sigma[[l]])) {
      stop("The dimension of each Sigma must coincide with the corresponding dimension provided in indexB")
    }
  }
  if (length(time) != TT) {
    stop("The length of time and the dimension of X must coincide.")
  }
  if (length(time) != nrow(B)) {
    stop("The length of time and the dimension of B must coincide.")
  }

  colnames(X) <- time
  X <- as.matrix(X) # If not already done, convert it onto a matrix
  index_not_NA <- !is.na(X) # Identify the missing values
  nTT <- sum(index_not_NA) # Number of observed values
  n_time <- colSums(!is.na(X)) # How many not-missing values per column? (T_i in the paper notation)
  n_subject <- rowSums(!is.na(X)) # How many values per rows (each time-grid)?

  # VECTORIZATION
  Y <- c(t(X))[c(t(index_not_NA))] # Columns are ordered by subject
  # BB <- B; for(i in 2:n) {BB <- rbind(BB,B)} # It is pretty slow, although the operation is simple
  BB <- t(matrix(t(B), M, n * TT)) # This essentially avoid the above for-loop and it is faster.
  BB <- BB[c(t(index_not_NA)), ] # This select only the relevant values

  # Penalty terms
  traceBSigmaB <- matrix(0, TT, L)
  PSigma <- list()
  Pmu <- list()
  mu_tilde <- list()
  Sigma_tilde <- list()

  for (l in 1:L) {
    PSigma[[l]] <- solve(Sigma[[l]]) # Store the inverse of Sigmas
    traceBSigmaB[, l] <- rowSums(B[, indexB[l, ]] %*% Sigma[[l]] * B[, indexB[l, ]])

    mu_tilde[[l]] <- matrix(0, Hl[l], Ml[l])
    Sigma_tilde[[l]] <- array(0, c(Hl[l], Ml[l], Ml[l]))
  }

  # Verbose settings and output
  verbose_step <- 1

  rho <- matrix(0, n, H)
  rho_class <- matrix(0, n, L)
  predH <- matrix(0, H, TT) # Here everything is stacked
  E_residuals <- array(0, c(H, n, TT)) # Here as well
  E_logprob <- numeric(H)
  alpha_p <- list()
  alpha_p_tilde <- list()
  E_logp <- list()

  for (l in 1:L) {
    alpha_p[[l]] <- rep(cc[l] / Hl[l], Hl[l])
    alpha_p_tilde[[l]] <- alpha_p[[l]]
    E_logp[[l]] <- digamma(alpha_p_tilde[[l]]) - digamma(sum(alpha_p_tilde[[l]]))
  }

  alpha_tilde <- alpha
  E_logPi <- digamma(alpha_tilde) - digamma(sum(alpha_tilde))

  # Initialization of different pieces of the lowerbound
  lower1 <- lower4 <- lower7 <- lower8 <- lower9 <- lower10 <- 0
  lower3 <- lower5 <- numeric(H)
  lower2 <- lower6 <- numeric(L)
  lowerbound <- -Inf

  # -------------------------------------------------
  # Initialization settings
  # -------------------------------------------------

  if(verbose){ cat("Pre-allocating observations into groups...\n")}

  # Divide observations in cluster and perturbate the weights a bit.
  pre_clust  <- FB_clust(X=X,  Hl=Hl, L =L, B = B, indexB = indexB, time = time, prediction = TRUE)
  G          <- as.factor(pre_clust$cluster)
  G          <- factor(as.numeric(factor(G,levels(G)[order(table(G),decreasing = TRUE)])),levels=1:H)
  rho <- cbind(model.matrix(rep(1,n) ~ G - 1))

  # Mixing weights are essentially assigned at random
  # rho <- matrix(runif(n * H), n, H)
  rho <- rho / rowSums(rho)
  sums_rho <- colSums(rho)

  # This is a raw "estimate" for the parameters of the variance. The weights are not used at all
  a_sigma2_tilde <- a_sigma + nTT / 2
  b_sigma2_tilde <- b_sigma + sum((X - mean(X))^2, na.rm = TRUE) / 2

  E_tau <- a_sigma2_tilde / b_sigma2_tilde # Expected value of tau
  E_ltau <- digamma(a_sigma2_tilde) - log(b_sigma2_tilde)

  # -------------------------------------------------------
  # Starting the CAVI algorithm
  # -------------------------------------------------------

  if (verbose) {
    cat("Starting the Variational Bayes algorithm...\n")
  }
  for (r in 1:maxiter) {

    # Update the regression coefficients
    for (l in 1:L) {
      for (h in 1:Hl[l]) {
        # Set the h_id
        h_id <- sum(Hl[1:l]) - Hl[l] + h

        #        if (sum(rho[, h_id]) < 1e-12) { # Here it is h_id
        #          Sigma_tilde[[l]][h, , ] <- Sigma[[l]] # Here is just h
        #          mu_tilde[[l]][h, ] <- rep(0, Ml[l])
        #          predH[h_id, ] <- rep(0, TT)

        #  # Allocating the residuals
        #  for (i in 1:n) {
        #    E_residuals[h_id, i, ] <- X[i, ]^2 - 2 * X[i, ] * predH[h_id, ] + predH[h, ]^2 + traceBSigmaB[, l]
        #  }
        # } else {
        weights <- rep(rho[, h_id], n_subject) # This associate to each subject / observation a weight

        Sigma_tilde[[l]][h, , ] <- solve(E_tau * crossprod(BB[, indexB[l, ]] * sqrt(weights)) + PSigma[[l]])
        mu_tilde[[l]][h, ] <- Sigma_tilde[[l]][h, , ] %*% (crossprod(BB[, indexB[l, ]] * weights * E_tau, Y))

        # Predictions (the same for elements in the same cluster)
        predH[h_id, ] <- as.numeric(B[, indexB[l, ]] %*% mu_tilde[[l]][h, ])

        # Inefficient
        for (i in 1:n) {
          E_residuals[h_id, i, ] <- X[i, ]^2 - 2 * X[i, ] * predH[h_id, ] + predH[h_id, ]^2 + rowSums(B[, indexB[l, ]] %*% Sigma_tilde[[l]][h, , ] * B[, indexB[l, ]])
        }
        # }
      }
    }

    # Update Pi and p, for l=1,...,L
    for (l in 1:L) {
      alpha_tilde[l] <- alpha[l] + sum(rho[, which(Findex == l)])
      for (h in 1:Hl[l]) {
        h_id <- sum(Hl[1:l]) - Hl[l] + h
        alpha_p_tilde[[l]][h] <- alpha_p[[l]][h] + sum(rho[, h_id])
      }
      E_logp[[l]] <- digamma(alpha_p_tilde[[l]]) - digamma(sum(alpha_p_tilde[[l]]))
    }
    E_logPi <- digamma(alpha_tilde) - digamma(sum(alpha_tilde))

    # Store E_logprobs
    for (l in 1:L) {
      for (h in 1:Hl[l]) {
        h_id <- sum(Hl[1:l]) - Hl[l] + h
        E_logprob[h_id] <- E_logPi[l] + E_logp[[l]][h]
      }
    }

    # Variance
    a_sigma2_tilde <- a_sigma + 0.5 * nTT
    b_sigma2_tilde <- b_sigma + 0.5 * sum(rho * t(apply(E_residuals, c(1, 2), function(x) sum(x, na.rm = TRUE))))
    E_tau <- a_sigma2_tilde / b_sigma2_tilde # Expected value of tau
    E_ltau <- digamma(a_sigma2_tilde) - log(b_sigma2_tilde) # Expected value of its logarithm

    # Step  - Cluster allocation
    rho <- rho_update(X, E_logprob, E_residuals, E_tau)
    rho <- rho + 1e-16 # This step is needed for numerical reasons.
    rho <- rho / rowSums(rho)
    sums_rho <- colSums(rho)

    # Lower-bound for the likelihood
    lower1 <- sum(sums_rho * E_logprob) + 0.5 * nTT * E_ltau - 0.5 * E_tau * sum(rho * t(apply(E_residuals, c(1, 2), function(x) sum(x, na.rm = TRUE))))

    # Lowerbound for the prior
    lower4 <- -b_sigma * E_tau + (a_sigma - 1) * E_ltau

    # Lower-bound for the variational approximation
    lower7 <- a_sigma2_tilde * log(b_sigma2_tilde) - lgamma(a_sigma2_tilde) - b_sigma2_tilde * E_tau + (a_sigma2_tilde - 1) * E_ltau

    # Lower-bound for the discrete distribution
    lower8 <- sum(rho * log(rho))

    for (l in 1:L) {
      for (h in 1:Hl[l]) {
        h_id <- sum(Hl[1:l]) - Hl[l] + h
        # Lower bound for the priors
        lower3[h_id] <- -0.5 * (sum(diag(PSigma[[l]] %*% Sigma_tilde[[l]][h, , ])) + t(mu_tilde[[l]][h, ]) %*% PSigma[[l]] %*% (mu_tilde[[l]][h, ]))
        # lower3[h]    <-  - 0.5*lambda*(crossprod(mu_tilde[h,]) + sum(diag(Sigma_tilde[h,,])))
        lower5[h_id] <- -0.5 * ldet(Sigma_tilde[[l]][h, , ])
      }
    }

    for (l in 1:L) {
      # Lower bound for the priors
      lower2[l] <- sum((alpha_p[[l]] - 1) * E_logp[[l]])
      # Lower-bound for the variational distributions
      lower6[l] <- lgamma(sum(alpha_p_tilde[[l]])) - sum(lgamma(alpha_p_tilde[[l]])) + sum((alpha_p_tilde[[l]] - 1) * E_logp[[l]])
    }

    # Lower bound for the priors
    lower9 <- sum((alpha - 1) * E_logPi)
    # Lower-bound for the variational distributions
    lower10 <- lgamma(sum(alpha_tilde)) - sum(lgamma(alpha_tilde)) + sum((alpha_tilde - 1) * E_logPi)

    # Convergence checks
    lowerbound_new <- lower1 + sum(lower2) + sum(lower3) + lower4 + lower9 - sum(lower5) - sum(lower6) - lower7 - lower8 - lower10

    # Break the loop at convergence
    if (lowerbound_new - lowerbound < tol) {
      if (verbose) {
        cat(paste("Convergence reached after", r, "iterations. \n"))
      }
      break
    }

    # Otherwise continue
    lowerbound <- lowerbound_new

    # Display status
    if (verbose) {
      if (r %% verbose_step == 0) {
        cat(paste("Lower-bound: ", round(lowerbound, 7), ", iteration: ", r, "\n", sep = ""))
      }
    }
  }

  if (r == maxiter) {
    warning(paste("Convergence has not been reached after", r, "iterations. \n"))
  }

  # Output
  cluster <- factor(apply(rho, 1, which.max), levels = 1:H)

  for (l in 1:L) {
    rho_class[, l] <- rowSums(rho[, which(Findex == l)])
  }
  Fclass <- factor(apply(rho_class, 1, which.max), levels = 1:L)

  # Individual curves
  pred <- matrix(0, n, TT)
  for (i in 1:n) {
    pred[i, ] <- colSums(predH * rho[i, ])
  }

  # table(cluster)
  out <- list(
    mu_tilde = mu_tilde, Sigma_tilde = Sigma_tilde,
    a_tilde_sigma = a_sigma2_tilde, b_tilde_sigma = b_sigma2_tilde, prediction = pred,
    pred_cluster = predH, cluster = cluster, Fclass = Fclass, Findex = Findex, rho = rho, rho_class = rho_class, alpha_tilde = alpha_tilde, alpha_p_tilde = alpha_p_tilde, lowerbound = lowerbound
  )
  attr(out, "class") <- "EFDMP"
  return(out)
}


FB_clust <- function(X, Hl, L, B, indexB, time = NULL, prediction = TRUE, ...) {

  # Fixed quantities
  n <- NROW(X)
  TT <- NCOL(X)
  M <- ncol(B)
  Ml <- rowSums(indexB)

  if (is.null(time)) {
    time <- 1:TT
    warning("The time dimension was not supplied and equally-spaced intervals are assumed.")
  }
  if (length(time) != TT) {
    stop("The length of time and the dimension of X must coincide.")
  }
  if (length(time) != nrow(B)) {
    stop("The length of time and the dimension of B must coincide.")
  }

  colnames(X) <- time
  beta <- NULL

  cluster <- numeric(n)
  Fclass <- numeric(n)
  RMSE <- numeric(L)

  out <- list()

  if (prediction) {
    smoothed_values <- matrix(0, n, TT)
    pred <- matrix(0, sum(Hl), TT)
  }
  # Essentially I am calculating these things twice.
  for (i in 1:n) {
    id <- !is.na(as.matrix(X[i, ]))
    for (l in 1:L) {
      beta_temp <- solve(crossprod(B[id, indexB[l, ]]), crossprod(B[id, indexB[l, ]], X[i, id]))
      RMSE[l] <- sqrt(mean((X[i, id] - c(B[id, indexB[l, ]] %*% beta_temp))^2))
    }
    Fclass[i] <- which.min(RMSE)
  }

  for (l in 1:L) {
    beta <- matrix(0, sum(Fclass == l), Ml[l])
    idF <- which(Fclass == l)
    for (i in 1:sum(Fclass == l)) {
      beta[i, ] <- solve(crossprod(B[id, indexB[l, ]]), crossprod(B[id, indexB[l, ]], X[idF[i], id]))
      if (prediction) {
        smoothed_values[idF[i], ] <- c(B[id, indexB[l, ]] %*% beta[i, ])
      }
    }
    km <- kmeans(beta, min(Hl[l], nrow(beta) - 1))
    cluster[idF] <- sum(Hl[1:l]) - Hl[l] + km$cluster
    if (prediction) {
      for (h in 1:min(Hl[l], nrow(beta) - 1)) {
        pred[sum(Hl[1:l]) - Hl[l] + h, ] <- as.numeric(B[id, indexB[l, ]] %*% km$centers[h, ])
      }
    }
  }

  # Output
  out$Fclass <- Fclass
  out$cluster <- cluster

  if (prediction) {
    out$pred_cluster <- pred
    out$prediction <- smoothed_values
  }
  return(out)
}
