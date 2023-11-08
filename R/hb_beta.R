#' Small Area Estimation using Hierarchical Bayesian under Beta Distribution
#'
#' @param formula an object of class formula that contains a description of the model to be fitted. The variables included in the formula must be contained in the data.
#' @param iter.update Number of updates with default 3
#' @param iter.mcmc Number of total iterations per chain with default 10000
#' @param coef a vector contains prior initial value of Coefficient of Regression Model for fixed effect with default vector of 0 with the length of the number of regression coefficients
#' @param var.coef a vector contains prior initial value of variance of Coefficient of Regression Model with default vector of 1 with the length of the number of regression coefficients
#' @param thin Thinning rate, must be a positive integer with default 2
#' @param burn.in Number of iterations to discard at the beginning with default 2000
#' @param tau.u Prior initial value of inverse of Variance of area random effect with default 1
#' @param data data frame containing the variables named in \code{formula}
#' @param seed number used to initialize a pseudorandom number generator (default seed = 1). The random number generator method used is "base::Wichmann-Hill".
#' @param quiet if TRUE, then messages generated during compilation will be suppressed (default TRUE).
#' @param plot if TRUE, the autocorrelation, trace, and density plots will be generated (default TRUE).
#'
#' @return The function returns a list with the following objects : Estimation \code{Est}, random effect variance \code{refVar}, beta coefficient \code{Coefficient} and MCMC result \code{result_mcmc}
#' @export
#'
#' @examples
#' \dontrun{
#'  model1 <- hb_beta(
#'     y ~ x1 + x2 + x3,
#'     data = unemployment, seed = 1,
#'     burn.in = 2000, iter.mcmc = 10000, iter.update = 10,
#'  )
#' }
hb_beta <- function(formula, iter.update = 3, iter.mcmc = 10000, coef,
                    var.coef, thin = 3, burn.in = 2000, tau.u = 1, data, seed = 1, quiet = TRUE, plot = TRUE) {

  result <- list(Est = NA, refVar = NA, coefficient = NA, plot = NA)
  formuladata <- stats::model.frame(formula, data, na.action = NULL)

  if (any(is.na(formuladata[, -1]))) {
    stop("Auxiliary Variables contains NA values.")
  }
  auxVar <- as.matrix(formuladata[, -1])
  nvar <- ncol(auxVar) + 1

  if (!missing(var.coef)) {
    if (length(var.coef) != nvar) {
      stop("length of vector var.coef does not match the number of regression coefficients, the length must be ", nvar)
    }
    tau.b.value <- 1 / var.coef
  } else {
    tau.b.value <- 1 / rep(1, nvar)
  }

  if (!missing(coef)) {
    if (length(coef) != nvar) {
      stop("length of vector coef does not match the number of regression coefficients, the length must be ", nvar)
    }
    mu.b.value <- coef
  } else {
    mu.b.value <- rep(0, nvar)
  }

  if (iter.update < 3) {
    stop("the number of iteration updates at least 3 times")
  }


  if (quiet) {
    cli::cli_progress_bar(
      total = iter.update,
      format = "Update {iter}/{iter.update} | {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta} \n"
    )
    my_env <- environment()
    default_pb <- "none"
    update_progress <- function(iter, iter.update){
      Sys.sleep(2/10000)
      cli::cli_progress_update(set = iter, .envir = my_env)
    }
  }else{
    default_pb <- "text"
    update_progress <- function(iter, iter.update){
      cli::cli_h1('Update {iter}/{iter.update}')
    }
  }

  # rjags::load.module("lecuyer")
  if (!any(is.na(formuladata[, 1]))) {
    formuladata <- as.matrix(stats::na.omit(formuladata))
    if (any(formuladata[, 1] <= 0) || any(formuladata[, 1] >= 1)) {
      stop("response variable must be 0 < ", formula[2]," < 1")
    }

    x <- stats::model.matrix(formula, data = as.data.frame(formuladata))
    x <- as.matrix(x)
    n <- nrow(formuladata)
    mu.b <- mu.b.value
    tau.b <- tau.b.value
    tau.ub <- tau.ua <- phi.aa <- phi.ab <- phi.ba <- phi.bb <- a.var <- 1

    for (iter in 1:iter.update) {
      update_progress(iter, iter.update)

      dat <- list(
        n = n, nvar = nvar, y = formuladata[,1],
        x = as.matrix(formuladata[, 2:nvar]), mu.b = mu.b,
        tau.b = tau.b, tau.ua = tau.ua, tau.ub = tau.ub,
        phi.aa = phi.aa, phi.ab = phi.ab, phi.ba = phi.ba,
        phi.bb = phi.bb
      )
      inits <- list(
        u = rep(0, n), b = mu.b, tau.u = tau.u,
        .RNG.name = "base::Wichmann-Hill",
        .RNG.seed = seed
      )

      my_model <- "
        model{
  					for (i in 1:n) {
  							y[i] ~ dbeta(A[i],B[i])
  							A[i] <- mu[i] * phi[i]
                B[i] <- (1-mu[i]) * phi[i]
                logit(mu[i]) <- b[1] + sum(b[2:nvar]*x[i,]) + u[i]
  							u[i] ~ dnorm(0,tau.u)
  							phi[i] ~ dgamma(phi.a,phi.b)
  					}

  					for (k in 1:nvar){
  					    b[k] ~ dnorm(mu.b[k],tau.b[k])
  					}

  					phi.a ~ dgamma(phi.aa,phi.ab)
  					phi.b ~ dgamma(phi.ba,phi.bb)
  					tau.u ~ dgamma(tau.ua,tau.ub)
  					a.var <- 1 / tau.u

  			}"
      jags.m <- rjags::jags.model(
        textConnection(my_model),
        data = dat,
        inits = inits, n.chains = 1, n.adapt = 500,
        quiet = quiet
      )

      params <- c("mu", "a.var", "b", "phi.a", "phi.b", "tau.u")
      samps <- rjags::coda.samples(
        jags.m, params,
        n.iter = iter.mcmc,
        thin = thin, seed = seed,
        progress.bar = default_pb
      )

      samps1 <- stats::window(samps, start = burn.in + 1, end = iter.mcmc)
      result_samps <- summary(samps1)
      a.var <- result_samps$statistics[1]
      beta <- result_samps$statistics[2:(nvar + 1), 1:2]

      mu.b <- beta[, 1]
      tau.b <- 1 / (beta[, 2]^2)

      phi.aa <- result_samps$statistics[2 + nvar + n, 1]^2 / result_samps$statistics[2 + nvar + n, 2]^2
      phi.ab <- result_samps$statistics[2 + nvar + n, 1] / result_samps$statistics[2 + nvar + n, 2]^2
      phi.ba <- result_samps$statistics[3 + nvar + n, 1]^2 / result_samps$statistics[3 + nvar + n, 2]^2
      phi.bb <- result_samps$statistics[3 + nvar + n, 1] / result_samps$statistics[3 + nvar + n, 2]^2
      tau.ua <- result_samps$statistics[4 + nvar + n, 1]^2 / result_samps$statistics[4 +  nvar + n, 2]^2
      tau.ub <- result_samps$statistics[4 + nvar + n, 1] / result_samps$statistics[4 + nvar + n, 2]^2
    }

    result_samps <- summary(samps1)
    b.varnames <- c('intercept', all.vars(formula[[3]]))

    result_mcmc <- samps1[, c(2:(nvar + 1))]
    colnames(result_mcmc[[1]]) <- b.varnames
    a.var <- result_samps$statistics[1]
    beta <- result_samps$statistics[2:(nvar + 1), 1:2]
    rownames(beta) <- b.varnames
    mu <- result_samps$statistics[(nvar + 2):(1 + nvar + n), 1:2]

    Estimation <- data.frame(mu)
    Quantiles <- as.data.frame(result_samps$quantiles[1:(2 + nvar + n), ])
    q_mu <- Quantiles[(nvar + 2):(nvar + 1 + n), ]
    q_beta <- (Quantiles[2:(nvar + 1), ])
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta, q_beta)
    Estimation <- data.frame(Estimation, q_mu)
    colnames(Estimation) <- c(
      "MEAN", "SD", "2.5%", "25%",
      "50%", "75%", "97.5%"
    )
  } else {
    formuladata <- as.data.frame(formuladata)
    n <- nrow(formuladata)
    mu.b <- mu.b.value
    tau.b <- tau.b.value
    tau.ub <- tau.ua <- phi.aa <- phi.ab <- phi.ba <- phi.bb <- a.var <- 1
    formuladata$idx <- 1:n
    data_sampled <- stats::na.omit(formuladata)


    if (any(data_sampled[, 1] <= 0) || any(data_sampled[,1] >= 1)) {
      stop("response variable must be 0 < ", formula[2]," < 1")
    }
    data_nonsampled <- formuladata[-data_sampled$idx, ]
    r <- data_nonsampled$idx
    n1 <- nrow(data_sampled)
    n2 <- nrow(data_nonsampled)

    for (iter in 1:iter.update) {
      update_progress(iter, iter.update)

      dat <- list(
        n1 = n1, n2 = n2, nvar = nvar, y_sampled = data_sampled[,1],
        x_sampled = as.matrix(data_sampled[, 2:nvar]),
        x_nonsampled = as.matrix(data_nonsampled[, 2:nvar]),
        mu.b = mu.b, tau.b = tau.b, tau.ua = tau.ua,
        tau.ub = tau.ub, phi.aa = phi.aa, phi.ab = phi.ab,
        phi.ba = phi.ba, phi.bb = phi.bb
      )
      inits <- list(
        u = rep(0, n1), v = rep(0, n2), b = mu.b,
        tau.u = tau.u,
        .RNG.name = "base::Wichmann-Hill",
        .RNG.seed = seed
      )

      my_model <- "model {
					for (i in 1:n1) {
							y_sampled[i] ~ dbeta(A[i],B[i])
							A[i] <- mu[i] * phi[i]
              B[i] <- (1-mu[i]) * phi[i]
              logit(mu[i]) <- b[1] +  sum(b[2:nvar]*x_sampled[i,]) + u[i]
							u[i] ~ dnorm(0,tau.u)
							phi[i] ~ dgamma(phi.a,phi.b)

					}

					for (j in 1:n2) {
				    	y.nonsampled[j] ~ dbeta(A.nonsampled[j],B.nonsampled[j])
							A.nonsampled[j] <- mu.nonsampled[j] * phi.nonsampled[j]
              B.nonsampled[j] <- (1-mu.nonsampled[j]) * phi.nonsampled[j]
              logit(mu.nonsampled[j]) <- mu.b[1] +sum(mu.b[2:nvar]*x_nonsampled[j,])+v[j]
					    v[j]~dnorm(0,tau.u)
					    phi.nonsampled[j] ~ dgamma(phi.a,phi.b)

          }

					for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					}

					phi.a ~ dgamma(phi.aa,phi.ab)
					phi.b ~ dgamma(phi.ba,phi.bb)
					tau.u ~ dgamma(tau.ua,tau.ub)
					a.var <- 1 / tau.u

			}"
      jags.m <- rjags::jags.model(
        file = textConnection(my_model), data = dat,
        inits = inits, n.chains = 1, n.adapt = 500,
        quiet = quiet
      )

      params <- c(
        "mu", "mu.nonsampled", "a.var", "b",
        "phi.a", "phi.b", "tau.u"
      )
      samps <- rjags::coda.samples(
        jags.m, params,
        n.iter = iter.mcmc,
        thin = thin,
        progress.bar = default_pb
      )
      samps1 <- stats::window(samps, start = burn.in + 1, end = iter.mcmc)
      result_samps <- summary(samps1)
      a.var <- result_samps$statistics[1]
      beta <- result_samps$statistics[2:(nvar + 1), 1:2]

      mu.b <- beta[, 1]
      tau.b <- 1 / (beta[, 2]^2)

      phi.aa <- result_samps$statistics[2 + nvar + n, 1]^2 / result_samps$statistics[2 + nvar + n, 2]^2
      phi.ab <- result_samps$statistics[2 + nvar + n, 1] / result_samps$statistics[2 + nvar + n, 2]^2
      phi.ba <- result_samps$statistics[3 + nvar + n, 1]^2 / result_samps$statistics[3 + nvar + n, 2]^2
      phi.bb <- result_samps$statistics[3 + nvar + n, 1] / result_samps$statistics[3 + nvar + n, 2]^2
      tau.ua <- result_samps$statistics[4 + nvar + n, 1]^2 / result_samps$statistics[4 + nvar + n, 2]^2
      tau.ub <- result_samps$statistics[4 + nvar + n, 1] / result_samps$statistics[4 + nvar + n, 2]^2
    }

    result_samps <- summary(samps1)
    b.varnames <- c('intercept', all.vars(formula[[3]]))
    result_mcmc <- samps1[, c(2:(nvar + 1))]
    colnames(result_mcmc[[1]]) <- b.varnames
    a.var <- result_samps$statistics[1]
    beta <- result_samps$statistics[2:(nvar + 1), 1:2]
    rownames(beta) <- b.varnames
    mu <- result_samps$statistics[(nvar + 2):(1 + nvar + n1),1:2]
    mu.nonsampled <- result_samps$statistics[(2 + nvar + n1):(1 + nvar + n), 1:2]
    Estimation <- matrix(rep(0, n), n, 2)
    Estimation[r, ] <- mu.nonsampled
    Estimation[-r, ] <- mu
    Estimation <- as.data.frame(Estimation)
    Quantiles <- as.data.frame(result_samps$quantiles[1:(3 +
                                                           nvar + n), ])
    q_beta <- (Quantiles[2:(nvar + 1), ])
    q_mu <- Quantiles[(nvar + 2):(nvar + 1 + n1), ]
    q_mu.nonsampled <- (Quantiles[(2 + nvar + n1):(1 + nvar + n), ])
    q_Estimation <- matrix(0, n, 5)
    for (i in 1:5) {
      q_Estimation[r, i] <- q_mu.nonsampled[, i]
      q_Estimation[-r, i] <- q_mu[, i]
    }
    # q_Estimation[r, 1:5] <- q_mu.nonsampled[, 1:5]
    # q_Estimation[-r, 1:5] <- q_mu[, 1:5]

    rownames(q_beta) <- b.varnames
    beta <- cbind(beta, q_beta)
    Estimation <- data.frame(Estimation, q_Estimation)
    colnames(Estimation) <- c(
      "MEAN", "SD", "2.5%", "25%",
      "50%", "75%", "97.5%"
    )
  }
  result$Est <- Estimation
  result$refVar <- a.var
  result$coefficient <- beta
  result$result_mcmc <- result_mcmc
  result$dat <- dat
  result$inits <- inits
  result$formula <- formula
  result$formuladata <- formuladata
  class(result) <- 'saehb'

  if (plot) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))

    graphics::par(mar = c(2, 2, 2, 2))
    coda::autocorr.plot(result_mcmc, col = "brown2", lwd = 2)
    plot(result_mcmc, col = "brown2", lwd = 2)
  }

  cli::cli_h1('Coefficient')
  stats::printCoefmat(beta)
  return(result)
}



