#' @useDynLib wglmmModularised
#' @title Modularized Survey-weighted Generalized Linear Mixed Models
#' @aliases wglmm
#' @description Modularized function to estimate generalized linear mixed models with flexible random effects structure and survey weights incorporating options for E-Step and M-Step.
#' @usage wglmm(formula, data = NULL, family = gaussian(), weights = NULL,
#'   iter1  = 501, iter2 = 501,  MI = 500, tol1 = 2e-4,
#'   tol2 = 1e-5, trace = TRUE, nDecrease = 2, ...)
#' @param formula Formula for the Model to be estimated. Has the shape of lmer in the package lme4.
#' @param data Dataframe containing the variables
#' @param family Distribution from the exponential family. You can choose between gaussian(), binomial(),  poisson(),  Gamma() and  exponential(). However, only  gaussian() and binomial() are already extensively tested.
#' @param weights Survey weights
#' @param e_step_method Method for the E-Step ("e_step_importance_sampling" or "e_step_gibbs_sampling").
#' @param m_step_method Method for the M-Step ("m_step_gradient_descent" or "m_step_bfgs").
#' @param iter1 Maximum number of MCEM iterations
#' @param iter2 Maximum number of internal optimizations within a MCEM-step
#' @param MI Number of importance-sampled random numbers in E-step
#' @param tol1 Convergence criterion for MCEM algorithm
#' @param tol2 Convergence criterion for optimizations within a MCEM-step
#' @param ... Additional options for model set-up
#' @details For comparability between runs, set a seed due to stochastic optimization.
#' @export
#' @return A list with the following elements
#' \itemize{
#' \item coef Vector of fixed effects regression parameters
#' \item VarCov List of random effects variance-covariance matrices
#' \item scale Scale parameter for exponential family
#' \item RE_mat Matrix with simulated random effects from last E-step, including importance weights
#' \item RE_mode List with modes of the random effects
#' \item residuals List with different types of residuals
#' \item LLmod Joint maximum log-likelihood of observed data and the mode of random effects
#' \item LLexp Expected log-likelihood given the observed data
#' \item niter Number of iterations
#' \item convergence Has MCEM converged?
#' } or
#' \itemize{
#' \item DIC Deviance Information Criterion, a measure of model fit.
#' \item fixed.formula Formula for the fixed effects in the model.
#' \item random.formula Formula for the random effects in the model.
#' \item residual.formula Formula for the residual variance structure.
#' \item solutions Posterior means, confidence intervals, and effective sample sizes for the fixed effects.
#' \item Gcovariances Posterior means, confidence intervals, and effective sample sizes for the random effect covariances.
#' \item Gterms Joint Posterior modes and confidence intervals for the variance-covariance parameters of the random effects.
#' \item Rcovariances Posterior means, confidence intervals, and effective sample sizes for the residual variances.
#' \item Rterms Posterior modes and confidence intervals for the variance-covariance parameters of the residuals.
#' \item cstats Convergence diagnostics, including potential scale reduction factor (PSRF) and effective sample size.
#' \item cutpoints Cutpoints used for ordinal responses.
#' \item theta_scale Scale parameters for the random effects.
#' } based on method selected.
#' @references Burgard, Jan Pablo and Doerr, Patricia (2018). Survey-weighted Generalized Linear Mixed Models. Trier University. Working Paper.
#' @examples
#' \dontrun{
#' library(mvtnorm)
#' n <- 600              # Total population size per simulation
#' K <- 50               # Number of simulations
#' beta <- c(4, -2, -1)
#' covmat <- matrix(0.5, ncol = 2, nrow = 2)
#' diag(covmat) <- c(0.7, 1.3)
#' g1 <- 20              # Number of domains for d1
#' g2 <- 30              # Number of domains for d2
#'
#' # Preallocate storage for simulation results
#' population_data <- vector("list", K)
#' sampled_data_uninformative <- vector("list", K)
#' sampled_data_informative <- vector("list", K)
#'
#' # Perform simulations
#' for (k in 1:K) {
#'   X1           <- rnorm(n, mean = 2)
#'   X2           <- rexp(n, rate = 1)
#'
#'   group1 <- rep(1:g1, length.out = n)
#'   group2 <- rep(1:g2, length.out = n)
#'   re1    <- rnorm(g1, sd = 2)
#'   re2    <- rmvnorm(g2, sigma = covmat)
#'
#'   modX  <- model.matrix( ~ X1 + X2)
#'   modZ1 <- model.matrix( ~ -1 + as.factor(group1))
#'   modZ2A <- model.matrix( ~ -1 + as.factor(group2))
#'   modZ2B <- model.matrix( ~ -1 + X1:as.factor(group2))
#'
#'   eta    <- modX %*% beta + modZ1 %*% re1 +
#'     modZ2A %*% as.vector(re2[,1]) +
#'     modZ2B %*% as.vector(re2[,2])
#'
#'   lin <- eta + rnorm(n, sd = 2.3)
#'   prob <- 1 / (1 + exp(-eta))
#'   bin <- rbinom(n, size = 1, prob = prob)
#'
#' dfsamp <- data.frame(group1 = group1, group2 = group2, X1 = X1, X2 = X2, lin = lin, bin = bin, eta = eta)
#'   population_data[[k]] <- dfsamp
#' }
#'
#' # Step 2: Sampling
#' for (k in 1:K) {
#'   data_k <- population_data[[k]]
#'
#'   # Non-informative sampling
#'   sampled_data_uninformative[[k]] <- do.call("rbind", lapply(split(data_k, data_k$group1), function(df) {
#'     probs <- df$X2 / sum(df$X2)
#'
#'     sampled_indices <- sample(nrow(df), 5, replace = FALSE, prob = probs)
#'     sampled_df <- df[sampled_indices, ]
#'
#'     sampled_df$inclusion_prob <- probs[sampled_indices]
#'     sampled_df$weights <- 1 / sampled_df$inclusion_prob
#'
#'     return(sampled_df)
#'   }))
#'
#'   out <- wglmm( lin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#'               data = subset(sampled_data_uninformative[[k]], select = -c(inclusion_prob, weights)), trace = TRUE, iter1 = 250,
#'               iter2 = 1001, MI = 2000, family = gaussian(), e_step_method = "e_step_importance_sampling", weights = sampled_data_uninformative[[k]]$weights)
#'
#'
#'   # Informative sampling
#'   sampled_data_informative[[k]] <- do.call("rbind", lapply(split(data_k, data_k$group1), function(df) {
#'     residuals <- abs(df$lin - df$eta)
#'
#'     quantiles <- quantile(residuals, probs = c(0.25, 0.5, 0.75))
#'
#'     probs <- ifelse(residuals <= quantiles[1], 0.1,
#'                     ifelse(residuals <= quantiles[2], 0.2, 0.4))
#'
#'     sampled_indices <- sample(nrow(df), 5, replace = FALSE, prob = probs)
#'     sampled_df <- df[sampled_indices, ]
#'
#'     sampled_df$inclusion_prob <- probs[sampled_indices] * (5/30)
#'     sampled_df$weights <- 1 / sampled_df$inclusion_prob
#'
#'     return(sampled_df)
#'   }))
#'
#'   out <- wglmm( lin ~ X1 + X2 + (1|group1) + (1+X1|group2),
#'               data = subset(sampled_data_informative[[k]], select = -c(inclusion_prob, weights)), trace = TRUE, iter1 = 250,
#'               iter2 = 1001, MI = 2000, family = gaussian(), e_step_method = "e_step_importance_sampling", weights = sampled_data_informative[[k]]$weights)
#'
#' }
#' }
#' @import lme4
#' @rawNamespace import(Matrix, except = c(cov2cor, toeplitz, update) )
#' @import stats

wglmm <- function(formula, data = NULL, family = gaussian(), weights = NULL,
                  e_step_method = "e_step_importance_sampling", m_step_method = "m_step_bfgs",
                  iter1 = 501, iter2 = 501, MI = 500, tol1 = 2e-4, tol2 = 1e-5,
                  trace = TRUE, nDecrease = 2, ...) {

  # Validate E-Step and M-Step method inputs
  if (!exists(e_step_method)) {
    stop("Function '", e_step_method, "' is not implemented.")
  }

  if (!exists(m_step_method)) {
    stop("Function '", m_step_method, "' is not implemented.")
  }

  # Getting the parameters for the function call
  mc <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("data", "weights", "subset", "na.action", "offset", "contrasts"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::get_all_vars)
  mf$formula <- subbars(formula)

  # Data Preparation
  dat <- eval(mf, parent.frame())
  rownames(dat) <- 1:nrow(dat)

  nam <- c("data", "weights", "subset", "na.action", "contrasts", "offset")[c(TRUE, TRUE, m[-(1:2)] != 0)]

  if (m[2] == 0 || length(unique(dat$weights)) == 1) {
    w <- NULL
  } else {
    w <- dat$weights
  }

  # Sets up the model environment, preparing the variables and model settings.
  model_data <- setup_model(formula, family, tol1, MI, dat, w, nam, mf)

  # MCMC-algorithm for GLMMs with choice of type for e-step
  if (e_step_method %in% c("e_step_gibbs_sampling")) {
    mcmc_glmm_results <- mcmc_glmm_algorithm(formula, family, iter1, dat)
    return(mcmc_glmm_results)
  }

  # MCEM-algorithm for GLMMs with choice of type for e-step and m-step
  mcem_glmm_results <- mcem_glmm_algorithm(model_data, formula, family, tol1, tol2, iter1, iter2, MI, nDecrease, trace, dat, nam, mf, e_step_method, m_step_method)


  # Finalizing model parameters, optimization, gradient and Hessian calculations, and convergence checking
  X <- model_data$X
  y <- model_data$y
  Z <- model_data$Z
  w <- model_data$w
  qvec <- model_data$qvec
  ni <- model_data$ni
  n.r <- model_data$n.r
  est_results <- finalize_model_estimation(mcem_glmm_results, family, tol1, tol2, iter2, X, y, Z, w, qvec, ni, n.r, nDecrease)

  # Gathering all relevant results, including residuals and model diagnostics, to provide a comprehensive summary of the model's performance and outcomes
  ngrps <- model_data$ngrps
  nam <- model_data$nam
  conLin <- model_data$conLin
  i <- mcem_glmm_results$i
  final_results <- compile_model_results(est_results, formula, family, X, y, w, Z, ngrps, qvec, ni, nam, conLin, i)

  return(final_results)
}

setup_model <- function(formula, family, tol1, MI, dat, w, nam, mf) {
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame(2))
  if( is.function(family)) family <- family()

  nam <- c( "formula", nam)
  nam2 <- c( "formula", "dat", "w", nam[-(1:3)])
  env <- lapply(nam2, function(x) get(x))
  names(env) <- nam

  mf$formula <- formula
  scaleNeu <- NULL

  if (isTRUE(all.equal(family, gaussian()))) {
    # weights stored as heteroskedasticity correction, but correct ordering of elements
    conLin <- do.call(lFormula, env )
    env$weights <- w
    env$formula <- nobars( formula )
    mod <- do.call(lm, env )
    scale <- summary(mod)$sigma^2
  } else {

    env$family <- family
    conLin <- do.call(glFormula, env )
    env$weights <- w
    env$formula <- nobars( formula )
    mod <- do.call(glm, env )
    scale <- scaleNeu <- 1
    if(family$family == "Gamma") scale <- 1/summary(mod)$dispersion
  }
  # beta <- fixef(mod)
  # CovREalt <- CovRE <- lapply(VarCorr(mod), as.matrix)
  beta <- coef(mod)
  rm(mod)

  nam <- sapply( strsplit( names(conLin$reTrms$Ztlist), "| ", fixed = TRUE), function(x) return(x[2]) )
  ni <- sapply(conLin$reTrms$flist, function(x) length(unique(x)))

  ni <- ni[match(nam, names(conLin$reTrms$flist))]
  ni <- as.integer(ni)

  qvec <- mapply( function(z, ns) nrow(z)/ns, conLin$reTrms$Ztlist, ni)
  qvec <- as.integer(qvec)
  q <- sum(qvec)
  if(family$family == "Gamma")  CovREalt <- CovRE <- lapply(qvec, function(q) diag(q)*0.01) else CovREalt <- CovRE <- lapply(qvec, function(q) diag(q)*0.1)

  n.r <- sum(qvec*ni)

  ngrps <- length(qvec)

  if(is.null(w)) w <- as.numeric(1) else w <- as.numeric( conLin$fr$`(weights)` )
  w <- w/sum(w)*length(w)

  X <- as.matrix(conLin$X)
  y <- conLin$fr[,as.character(formula[2])]
  if(is.matrix(y) && ncol(y) > 1) y <- y[,1]/rowSums(y)
  if(is.matrix(y)) y <- y[,1]
  if(is.factor(y)) y <- as.numeric(y)-1
  Z <- t( conLin$reTrms$Zt )

  n <- length(y)
  p <- length(beta)

  i <- 0
  eps <- 1.1 * tol1
  MI2 <- MI0 <- MI
  mi <- round(MI*0.1)

  return(list(beta = beta, scale = scale, scaleNeu = scaleNeu, CovRE = CovRE, CovREalt = CovREalt, X = X, y = y, Z = Z, w = w, n = n, n.r = n.r, p = p, i = i, eps = eps, MI2 = MI2, mi = mi, qvec = qvec, ni = ni, ngrps = ngrps, nam = nam, conLin = conLin))
}

mcem_glmm_algorithm <- function(model_data, formula, family, tol1, tol2, iter1, iter2, MI, nDecrease, trace, dat, nam, mf, e_step_method, m_step_method) {
  # MCEM Algorithm logic

  # Initialization:
  beta <- model_data$beta
  scale <- model_data$scale
  scaleNeu <- model_data$scaleNeu
  CovRE <- model_data$CovRE
  CovREalt <- model_data$CovREalt
  X <- model_data$X
  y <- model_data$y
  Z <- model_data$Z
  w <- model_data$w
  n <- model_data$n
  n.r <- model_data$n.r
  p <- model_data$p
  i <- model_data$i
  eps <- model_data$eps
  MI2 <- model_data$MI2
  mi <- model_data$mi
  qvec <- model_data$qvec
  ni <- model_data$ni
  ngrps <- model_data$ngrps

  phi.alt <- phi.altalt <- beta
  modus <- rep(0,n.r)

  Chol <- 0.1 * diag(n.r)
  out <- .GenModusNeu(modus, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, MaxIt = iter2, tol = tol2, Chol = Chol, scale)
  Chol <- out$invChol
  LLalt <- .loglikTot(out$modus, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, scale)
  LLallTimeBest <- LL <- LLalt
  phiallTimeBest <- phi.alt
  CovREallTimeBest <- CovRE
  scaleallTimeBest <- scale
  counter <- 0

  while(i < iter1 && max(eps) > tol1){
    # em-loop
    k <- 0
    increaseLL <- FALSE

    # e-step
    e_step_arguments <- list(family = family, tol2 = tol2, iter2 = iter2, MI = MI, scale = scale, scaleNeu = scaleNeu, CovRE = CovRE, X = X, y = y, Z = Z, w = w, n = n, n.r = n.r, i = i, MI2 = MI2, mi = mi, qvec = qvec, ni = ni, phi.alt = phi.alt, Chol = Chol, out = out, LLalt = LLalt, LL = LL, k = k, increaseLL = increaseLL)
    e_step_results = do.call(e_step_method, e_step_arguments)

    # m-step
    LL = e_step_results$LL
    increaseLL = e_step_results$increaseLL
    LLalt = e_step_results$LLalt
    phi.neu = e_step_results$phi.neu
    outRE = e_step_results$outRE
    scaleNeu = e_step_results$scaleNeu

    m_step_arguments <- list(family = family, tol1 = tol1, tol2 = tol2, iter1 = iter1, iter2 = iter2, nDecrease = nDecrease, eps = eps, scale = scale, scaleNeu = scaleNeu, CovRE = CovRE, CovREalt = CovREalt, outRE = outRE, X = X, y = y, Z = Z, w = w, i = i, MI2 = MI2, mi = mi, qvec = qvec, ni = ni, phi.alt = phi.alt, phi.neu = phi.neu, Chol = Chol, out = out, LL = LL, increaseLL = increaseLL, LLallTimeBest = LLallTimeBest, phiallTimeBest = phiallTimeBest, CovREallTimeBest = CovREallTimeBest, scaleallTimeBest = scaleallTimeBest, counter = counter, trace = trace)
    m_step_results = do.call(m_step_method, m_step_arguments)

    LLalt <- m_step_results$LL
    scale <- m_step_results$scaleNeu
    phi.altalt <- m_step_results$phi.alt
    CovREalt <- m_step_results$CovRE
    CovRE <- m_step_results$outRE
    phi.alt <- m_step_results$phi.neu
    eps <- m_step_results$eps
    i <- m_step_results$i
    MI2 <- m_step_results$MI2
    Chol <- m_step_results$Chol
    out <- m_step_results$out
    LLallTimeBest <- m_step_results$LLallTimeBest
    phiallTimeBest <- m_step_results$phiallTimeBest
    CovREallTimeBest <- m_step_results$CovREallTimeBest
    scaleallTimeBest <- m_step_results$scaleallTimeBest
    counter <- m_step_results$counter
    br <- m_step_results$br

    if (br == TRUE) {
      break
    }
  }
  mcem_glmm_results <- m_step_results
  rm(e_step_results, m_step_results)
  return(mcem_glmm_results)
}

mcmc_glmm_algorithm <- function(formula, family, iter1, dat) {
  # MCMC Algorithm logic
  if (!require(MCMCglmm)) {
    install.packages("MCMCglmm", dependencies = TRUE)
  }
  library(MCMCglmm)

  if (isTRUE(all.equal(family, gaussian()))) {
    prior.lmm <- list(R = list(V = 1, nu = 0.002),  # For continous responses
                      G = list(G1 = list(V = diag(1), nu = 0.002),
                               G2 = list(V = diag(2), nu = 0.002)))

    # Edit formula by-hand
    fitLMM <- MCMCglmm(lin ~ X1 + X2 ,
                       random = ~ us(1):group1 + us(1+X1):group2,
                       rcov = ~idh(1):units,
                       prior = prior.lmm,
                       family = "gaussian",
                       data = dat,
                       nitt = iter1, burnin = 5, thin = 1, verbose = TRUE)

    return(summary(fitLMM))

  } else {
    prior.glmm <- list(R = list(V = 1, fix = 1),  # For binary responses
                       G = list(G1 = list(V = diag(1), nu = 0.002),
                                G2 = list(V = diag(2), nu = 0.002)))

    # Edit formula by-hand
    fitGLMM <- MCMCglmm(bin ~ X1 + X2 ,
                        random = ~ us(1):group1 + us(1+X1):group2,
                        rcov = ~idh(1):units,
                        family = "categorical",
                        prior = prior.glmm,
                        data = dat,
                        nitt = iter1, burnin = 5, thin = 1, verbose = TRUE)

    return(summary(fitGLMM))
  }
}

m_step_bfgs <- function(family, tol1, tol2, iter1, iter2, nDecrease, eps, scale, scaleNeu, CovRE, CovREalt, outRE, X, y, Z, w, i, MI2, mi, qvec, ni, phi.alt, phi.neu, Chol, out, LL, increaseLL, LLallTimeBest, phiallTimeBest, CovREallTimeBest, scaleallTimeBest, counter, trace) {
  tryCatch({
  # Maximization step with BFGS optimization logic

  br <- FALSE

  # Chol <- diag(n.r)/sqrt( sum((.GradientAll(out$modus, phi.neu, family$family, X, y, Z, outRE, w, qvec, ni, scaleNeu))^2) )
  out2 <- try( .GenModusNeu( out$modus, phi.neu, family$family, X, y, Z, outRE, w, qvec, ni, MaxIt = iter2, tol = tol2, Chol = Chol, scale = scale) )
  if(class(out2)=="try-error") break
  Chol <- out2$invChol

  MI2 <- MI2 + mi
  if(any(sapply(outRE, det)== 0)) break;

  if( increaseLL | i == 0){
    LLallTimeBest <- LL
    phiallTimeBest <- phi.neu
    CovREallTimeBest <- outRE
    scaleallTimeBest <- scaleNeu
  }else{
    counter <- counter + 1
  }
  if(counter == nDecrease){
    warning("Likelihood does not increase anymore!")
    break
  }

  out <- out2
  rm(out2)

  eps1 <- abs(phi.neu-phi.alt)/(abs(phi.alt) + 1)

  if( i > 0){
    eps2 <- unlist( mapply( function(mat1, mat2) as.vector(abs(mat1[lower.tri(mat1, diag = TRUE)]-mat2[lower.tri(mat1, diag = TRUE)])/(abs(mat2[lower.tri(mat2, diag = TRUE)]) + 1)), outRE, CovRE) )
    eps22 <- unlist( mapply( function(mat1, mat2) as.vector(abs(mat1[lower.tri(mat1, diag = TRUE)]-mat2[lower.tri(mat1, diag = TRUE)])/(abs(mat2[lower.tri(mat2, diag = TRUE)]) + 1)), outRE, CovREalt) )
  } else eps2 <- eps22 <- max(0.05, 1.1*tol1)
  eps <- c(eps1, eps2)

  if( MI2 %% 2 == 1 ) MI2 <- MI2 + 1

  i <- i+1

  if(trace){
    print(paste0("M-Step ", i, " of ", iter1, "  Maximum variable change in iteration: ", max(eps)*100, " %") )
    print(cbind(phi.alt, phi.neu) )
    print(mapply( function(mat1, mat2) cbind(mat1, mat2), CovRE, outRE, SIMPLIFY = FALSE))
    print(paste0("Increased LL: ", increaseLL))
  }
  return(list(LL = LL, scaleNeu = scaleNeu, phi.alt = phi.alt, phi.neu = phi.neu, CovRE = CovRE, outRE = outRE, eps = eps, i = i, MI2 = MI2, Chol = Chol, out = out, LLallTimeBest = LLallTimeBest, phiallTimeBest = phiallTimeBest, CovREallTimeBest = CovREallTimeBest, scaleallTimeBest = scaleallTimeBest, counter = counter, br = br))
  },error = function(e) {
    br <- TRUE
  }, finally = {
    return(list(LL = LL, scaleNeu = scaleNeu, phi.alt = phi.alt, phi.neu = phi.neu, CovRE = CovRE, outRE = outRE, eps = eps, i = i, MI2 = MI2, Chol = Chol, out = out, LLallTimeBest = LLallTimeBest, phiallTimeBest = phiallTimeBest, CovREallTimeBest = CovREallTimeBest, scaleallTimeBest = scaleallTimeBest, counter = counter, br = br))
  })
}

m_step_gradient_descent <- function(family, tol1, tol2, iter1, iter2, nDecrease, eps, scale, scaleNeu, CovRE, CovREalt, outRE, X, y, Z, w, i, MI2, mi, qvec, ni, phi.alt, phi.neu, Chol, out, LL, increaseLL, LLallTimeBest, phiallTimeBest, CovREallTimeBest, scaleallTimeBest, counter, trace) {
  # Gradient descent logic
  stop("M-Step method with 'Gradient descent logic' is not implemented.")

}

e_step_importance_sampling <- function(family, tol2, iter2, MI, scale, scaleNeu, CovRE, X, y, Z, w, n, n.r, i, MI2, mi, qvec, ni, phi.alt, Chol, out, LLalt, LL, k, increaseLL){
  # Expectation step with Importance sampling logic

  while(k < 10 & !increaseLL){
    Us <- .ImportanceSampling(MI = MI2 + k*mi, modus = out$modus, Chol = 1.05*Chol, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, scale)

    # Step 1 of optimization: Fixed effects
    outFix <- try( .NewCoefImp( Us[,-(n.r+1)], phi.alt, family$family, X, y, Z, weights = w, weightsImp = Us[,n.r+1], tol = tol2, MaxIt = iter2, scale) )
    if( class(outFix) == "try-error" || !outFix$convergence || any(is.infinite(outFix$phi)) || any(is.na(outFix$phi)) ) break else phi.neu <- outFix$phi

    # Step 2 of optimization: RE-covariance matrix
    outRE <- .NewCovREImp( Us[,-(n.r+1)], qvec, ni, weightsImp = Us[,n.r+1])

    if(family$family == "gaussian"){
      #out <- try( .GenModusNeu( out$modus, phi.neu, family$family, X, y, Z, outRE, w, qvec, ni, MaxIt = iter2, tol = tol2, Chol = out$invChol, scale = scale) )
      #if(class(out)=="try-error") break

      Xphi <- as.numeric(X %*% phi.neu)
      #linpred <- as.numeric(X %*% phi.neu + Z %*% out$modus)
      #res <- y - linpred
      scaleNeu <- sum( apply( Us[,1:n.r], 1, function(x) sum(w * (y - Xphi - as.numeric(Z %*% x))^2)/n ) * Us[,n.r+1] )
      #if(length(w) > 1) scaleNeu <- sum(w*res^2)/sum(w) else scaleNeu <- mean(res^2)
    }
    if(family$family == "Gamma"){
      # linpred <- as.numeric(X %*% phi.neu + Z %*% out2$modus)
      # res <- y - 1/linpred
      # if(length(w) > 1) scaleNeu <- sum(w)/sum(w*res^2*(linpred)^2) else scaleNeu <- 1/mean(res^2*(linpred)^2)
      Xphi <- as.numeric(X %*% phi.neu)
      scaleNeu <- sum( apply( Us[,1:n.r], 1, function(x){
        mu <- 1/(Xphi + as.numeric(Z %*% x))
        mu[mu <= 0] <- min(mu[mu > 0])*0.5
        # out <- sum(w*log(mu/y))*2/n # Deviance method
        out <- sum( w * (y - mu)^2/mu^2)/n # Pearson method
        # out <- sum( w * (y - mu)^2)/n
        return(out)} ) * Us[,n.r+1])
      # scaleNeu <- 1/scaleNeu
    }

    Us <- .ImportanceSampling(MI = MI2 + k*mi, modus = out$modus, Chol = 1.05*Chol, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, scale)

    LL <- .loglikTotImp(Us[,1:n.r], phi.neu, family$family, X, y, Z, outRE, w, qvec, ni, Us[,n.r+1], scaleNeu)
    LLalt <- .loglikTotImp(Us[,1:n.r], phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, Us[,n.r+1], scale)
    increaseLL <- LL < LLalt

    k <- k+1
    if(i == 0) break;
  }
  return(list(LL = LL, increaseLL = increaseLL, LLalt =  LLalt, phi.neu = phi.neu, outRE = outRE, scaleNeu = scaleNeu))
}

e_step_gibbs_sampling <- function(family, tol2, iter2, MI, scale, scaleNeu, CovRE, X, y, Z, w, n, n.r, i, MI2, mi, qvec, ni, phi.alt, Chol, out, LLalt, LL, k, increaseLL){
  # Gibbs sampling logic
  stop("E-Step method with 'Gibbs sampling logic' is implemented using 'mcmc_glmm_algorithm'.")

}

finalize_model_estimation <- function(mcem_glmm_results, family, tol1, tol2, iter2, X, y, Z, w, qvec, ni, n.r, nDecrease) {
  # Logic for convergence check

  # Initialization:
  phiallTimeBest <- mcem_glmm_results$phiallTimeBest
  CovREallTimeBest <- mcem_glmm_results$CovREallTimeBest
  scaleallTimeBest <- mcem_glmm_results$scaleallTimeBest
  out <- mcem_glmm_results$out
  CovRE <- mcem_glmm_results$CovRE
  MI2 <- mcem_glmm_results$MI2
  phi.alt <- mcem_glmm_results$phi.alt
  eps <- mcem_glmm_results$eps
  counter <- mcem_glmm_results$counter

  # Fixed effects editing
  coefs <- as.numeric(phiallTimeBest)
  CovRE <- CovREallTimeBest
  scale <- scaleallTimeBest

  out <- .GenModusNeu( out$modus, coefs, family$family, X, y, Z, CovRE, w, qvec, ni, MaxIt = iter2, tol = tol2, Chol = 1.05*out$invChol, scale)
  reHess1 <- tcrossprod(out$invChol)

  Us <- .ImportanceSampling(MI = MI2, modus = out$modus, Chol = 1.05*out$invChol, coefs, family$family, X, y, Z, CovRE, w, qvec, ni, scale = scale)

  names(coefs) <- names(beta)

  gr       <- as.numeric( .GradientBetaImp(Us[,-ncol(Us)], coefs, family$family, X, y, Z, w, weightsImp = Us[,ncol(Us)], scale = scale) )
  hess     <- .HessianBetaImp(Us[,-ncol(Us)], coefs, family$family, X, y, Z, w, weightsImp = Us[,ncol(Us)], scale = scale)
  gr2      <- apply(Us[,-ncol(Us)], 1, function(x) tcrossprod(.GradientBetaImp( matrix(x, nrow = 1), coefs, family$family, X, y, Z, w, 1, scale = scale)) )
  gr2      <- colSums( t(gr2) * Us[,ncol(Us)] )
  gr2      <- matrix(gr2, ncol = length(coefs))

  vcovBeta <- hess - gr2 + tcrossprod(gr)
  attr(coefs, "vcov") <- solve( vcovBeta )
  attr(coefs, "vcov.cstr") <- nearPD(solve(vcovBeta))$mat

  # Log Density of the Random Effects Vector
  LLmod <- -.loglikTot(out$modus, phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, scale)
  LLexp <- -.loglikTotImp(Us[,-(n.r+1)], phi.alt, family$family, X, y, Z, CovRE, w, qvec, ni, Us[,n.r+1], scale )

  convergence <- max(eps) <= tol1
  if( !convergence & counter != nDecrease) warning("No convergence!")

  return(list(coefs = coefs, out = out, reHess1 = reHess1, CovRE = CovRE, scale = scale, Us = Us, LLmod = LLmod, LLexp = LLexp, convergence = convergence))
}

compile_model_results <- function(est_results, formula, family, X, y, w, Z, ngrps, qvec, ni, nam, conLin, i) {
  # Output construction logic

  # Initialization:
  coefs <- est_results$coefs
  out <- est_results$out
  reHess1 <- est_results$reHess1
  CovRE <- est_results$CovRE
  scale <- est_results$scale
  Us <- est_results$Us
  LLmod <- est_results$LLmod
  LLexp <- est_results$LLexp
  convergence <- est_results$convergence

  lin.pred <- as.numeric(tcrossprod( X, t(coefs)))
  mu <- family$linkinv(lin.pred  + as.numeric(tcrossprod(Z, t(out$modus))))
  names(out$modus) <- colnames(Z)
  REs <- split(out$modus, rep(1:ngrps, ni*qvec))
  names(REs) <- nam

  resid <- sapply(0:length(REs), function(i) if(i == 0) return(numeric(nrow(X))) else return( as.numeric( crossprod(conLin$reTrms$Ztlist[[i]],  REs[[i]]) ) ) )
  resid <- t(apply(resid, 1, cumsum))

  REs  <- mapply( function(res, ns){
    m <- matrix(res, nrow = ns, byrow = TRUE)
    rownames(m) <- unique(names(res))
    return(m)}, REs, ni, SIMPLIFY = FALSE)
  nam2 <- names(conLin$reTrms$Ztlist)
  nam2 <- mapply(function(n1, n2) gsub(paste0(" [|] ", n1), "", n2), nam, nam2)
  nam2 <- gsub("1", "`(Intercept)'", nam2)
  nam2 <- strsplit(nam2, split = " + ", fixed = TRUE)
  REs  <- mapply(function(res, n){
    colnames(res) <- n
    return(res) }, REs, nam2, SIMPLIFY = FALSE)

  attr(REs, "invHess") <- reHess1

  pearson.resid <- (y - mu)/sqrt(family$variance(mu))
  if(!is.null(weights)) deviance.resid <- family$dev.resids(y, mu, w) else deviance.resid <- family$dev.resids(y, mu, rep(1, length(y)))

  resid <- y - family$linkinv(lin.pred + resid)
  if(!is.null(weights)) sigma <- sqrt(sum(w * resid[,ncol(resid)]^2)/sum(w)) else sigma <- sqrt(mean(resid[,ncol(resid)]^2))

  resid <- lapply(1:ncol(resid), function(i) resid[,i])
  attr(resid, "sigma") <- sigma

  V <- CovRE
  names(V) <- nam

  ret <- list(coef = coefs, VarCov = V, scale = scale, RE_mat = Us, RE_mode = REs, residuals = list(pearson = pearson.resid, deviance = deviance.resid, resid = resid), LLmod = LLmod, LLexp = LLexp, niter = i, convergence = convergence, formula = formula, family = family, data = conLin$fr)

  return(ret)
}
