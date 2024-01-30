softmax <- function(x) exp(c(0, x)) / sum(exp(c(0, x)))


LLinc <- function(beta, D, dimp, f0) {
  b     <- c(-log(sum(exp(D[, -1] %*% beta))), beta)
  p     <- as.vector(exp(D %*% b))
  denom <- Matrix::crossprod(p, dimp)
  -log(as.vector(denom)) %*% f0
}



fit <- function(formula, d0, lambda, stderr, cal,
                maxit, tol, seed, stepsize,
                lat, latlevels, lists, impute, lca, mse,
                numvars)
{

  set.seed(seed)

  # make d/f0 with observed rows/frequencies,
  # dobs excluding structural zeros and missings and,
  # (temporary) dcomp for completed data with latent variables (if any)
  if (mse)
  {
    s0   <- parse(text = paste(lists, "== 1", collapse = " | "))
    d    <- dplyr::filter(d0, eval(s0))
  } else {
    d    <- d0
  }

  f0    <- d$Freq
  n     <- sum(f0)
  dcomp <- dobs <- na.omit(d)
  dpop  <- na.omit(d0)

  # if lca, include latent variable(s) in dcomp
  if (lca)
  {
    nlat     <- length(lat)
    if (is.null(latlevels)) latlevels <- rep(2, nlat)
    nclasses <- prod(latlevels)
    prlat    <- softmax(-runif(1, .05, .5) - (2:nclasses) / nclasses)
    for (j in 2:nclasses)
    {
      dcomp <- rbind(dcomp, dobs)
      dpop  <- rbind(dpop, na.omit(d0))
    }

    X  <- matrix(0, nrow = nrow(dcomp), ncol = nlat, dimnames = list(NULL, sort(lat)))
    Y  <- matrix(0, nrow = nrow(dpop),  ncol = nlat, dimnames = list(NULL, sort(lat)))

    for (i in 1:nlat)
      {
      X[, i] <- rep(rep(1:latlevels[i], each = nrow(X) / prod(latlevels[1:i])), times = i)
      Y[, i] <- rep(rep(1:latlevels[i], each = nrow(Y) / prod(latlevels[1:i])), times = i)
    }
    dcomp <- cbind(mutate(as.data.frame(X), across(everything(), .fns = as.factor)), dcomp)
    dpop  <- cbind(mutate(as.data.frame(Y), across(everything(), .fns = as.factor)), dpop)
  }

  # make imputation matrix dimp
  if (impute)
  {
    mobs <- dobs %>% select(-Freq) %>% data.matrix()
    mmis <- d[!complete.cases(d), ] %>% select(-Freq) %>% data.matrix()

    imp  <- matrix(1, nrow(mobs), nrow(mmis))
    for (k in 1:nrow(mmis))
      for (j in 1:ncol(mobs))
        if (!is.na(mmis[k, j]))
          imp[, k] <- imp[, k] * (mobs[, j] == mmis[k, j])

    imp  <- cbind(Diagonal(nrow(mobs)), imp)

    if (lca)
    {
      dimp <- imp
      for (j in 2:nclasses)
        dimp <- rbind(dimp, imp)
    } else {
      dimp <- imp
    }
  } else {
    if (!lca) {
      dimp <- Diagonal(nrow(dobs))
    } else {
      imp <- Diagonal(nrow(dobs))
      for (j in 2:nclasses)
        dimp <- rbind(imp, Diagonal(nrow(dobs)))
    }
  }

  # design matrix
  D      <- Matrix::Matrix(model.matrix(formula, dcomp), sparse = T)

  # obtain random starting probabilities
  p <- if (lca) rep(prlat, each = nrow(dobs)) / nrow(dobs) else softmax(rnorm(nrow(dcomp) - 1, 0, .1))

  # compute starting values for dcomp and mu
  num    <- Matrix::t(p * dimp) * f0
  denom  <- Matrix::crossprod(p, dimp)
  mcomp  <- as.vector(denom^-1 %*% num)
  dcomp  <- mutate(dcomp, Freq = mcomp)
  b0     <- suppressWarnings(glm(formula = formula, poisson, mutate(dcomp, Freq = mcomp))$coef) / 2
  b      <- c(-log(sum(exp(D[, -1] %*% b0[-1]))), b0[-1])
  mu     <- n * as.vector(exp(D %*% b))

  iter   <- 0
  conv   <- FALSE
  LLold  <- log(as.vector(denom)) %*% f0
  hist   <- data.frame(iter = 0, LL = LLold, conv = NA)

  while (!conv)
  {
    db    <- Matrix::crossprod(D, mcomp) - Matrix::crossprod(D, mu) - 2 * lambda * c(0, b[-1])
    dmu   <- Diagonal(x = mu)
    I     <- Matrix::crossprod(D, dmu) %*% D + Diagonal(x = lambda * c(0, rep(2, length(b) - 1)))
    cvc   <- tryCatch(chol2inv(chol(I)), error = function(e) {
      stop("Hessian complete data loglikelihood not invertible, try increasing lambda")})
    bt    <- b + stepsize * cvc %*% db
    b     <- c(-log(sum(exp(D[, -1] %*% bt[-1]))), bt[-1])
    p     <- as.vector(exp(D %*% b))
    mu    <- n * p
    num   <- Matrix::t(p * dimp) * f0
    denom <- as.vector(Matrix::crossprod(p, dimp))
    mcomp <- as.vector(denom^-1 %*% num)

    LL    <- log(denom) %*% f0
    iter  <- iter + 1
    dLL   <- LL - LLold
    hist  <- rbind(hist, c(iter, LL, dLL))
    conv  <- dLL < tol
    LLold <- LL

    if (iter == maxit) break
  }

  # compute SE's based on complete data log-likelihood
  vcv <- chol2inv(chol(Matrix::crossprod(D, dmu) %*% D))
  se  <- sqrt(diag(vcv))

  # use incomplete data log-likelihood if requested and if impute and/or lca
  if (stderr[1] == "incomplete" & (impute | lca)){
    H  <- tryCatch(pracma::hessian(LLinc, b[-1], D = D, dimp = dimp, f0 = f0),
                   error = function(e) stderr <- "complete")
    if (is.matrix(H))
      vcv    <- chol2inv(chol(H))
      se[-1] <- sqrt(diag(vcv))
  }
  b[1] <- b[1] + log(n)

  coefs <- data.frame(beta = b, se,
                      zval = b / se,
                      pval = 2 * pnorm(-abs(b / se)),
                      row.names = colnames(D))

  fitted <- exp(model.matrix(formula, dpop) %*% b)
  Nhat   <- sum(fitted)
  dfct   <- mutate(dpop,
                   across(all_of(numvars), as.factor),
                   phat = fitted / Nhat) %>%
    select(-Freq)

  if (lca)
  {
    probs <- list()

    for(j in 1:nlat)
    {
      probs[[paste0("P(", colnames(dfct[j]), ")")]] <- xtabs(paste("phat ~", colnames(dfct[j])), data = dfct)
      q     <- NULL
      names <- NULL

      for (k in (length(lat) + 1):(ncol(dfct) - 1))
      {
        form  <- paste("phat ~ ", colnames(dfct)[k], "+", colnames(dfct)[j])
        q     <- rbind(q, xtabs(form, dfct)/matrix(colSums(xtabs(form, dfct)),
                                                   nlevels(dfct[, k]),
                                                   nlevels(dfct[, j]),
                                                   byrow = T))
        names <- c(names, rep(colnames(dfct)[k], nlevels(dfct[, k])))
      }

      rownames(q) <- paste(names, "=",  rownames(q))
      colnames(q) <- paste(colnames(dfct)[j], "=", 1:nlevels(dfct[, j]))
      probs[[paste0("P(.|", colnames(dfct)[j],")")]] <- q
    }
    hilo <- order(probs[[1]], decreasing = T)
    probs[[1]] <- probs[[1]][hilo]
    probs[[2]] <- probs[[2]][, hilo]

  }

  #############################################################
  cat("\n\n")
  print(cal)
  cat("\n\n")
  cat("MODEL \n")
  cat("=======================================================")
  cat("\n")
  cat("Imputation missing values   =", impute)
  cat("\n")
  cat("Latent variable(s)          =", lca)
  cat("\n")
  cat("Multiple Systems Estimation =", mse)
  cat("\n\n")
  cat("STATISTICS \n")
  cat("=======================================================")
  cat("\n")
  cat("Seed                  =", seed)
  cat("\n")
  cat("Lambda                =", lambda)
  cat("\n")
  cat("Sample size           =", n)
  cat("\n")
  cat("Iterations            =", iter)
  cat("\n")
  cat("Convergence criterion =", dLL)
  cat("\n")
  cat("Convergence           =", ifelse(dLL < 0, FALSE, conv))
  cat("\n")
  cat("Log-likelihood        =", LL)
  cat("\n")
  cat("Number of parameters  =", length(b))
  cat("\n")
  cat("Degrees of freedom    =", nrow(dcomp) - length(b))
  cat("\n")
  cat("Computation SE's      =", paste(stderr[1], "data log-likelihood"))
  cat("\n\n")
  cat("ESTIMATES \n")
  cat("=======================================================")

  if (mse)
  {
    cat("\n\n")
    cat("Population size estimate")
    cat("\n\n")
    print(data.frame(nobs = n, Nhat = Nhat, unobserved = Nhat - n))
  }

  if (lca){
    cat("\n\n")
    cat("Latent probabilities")
    cat("\n\n")
    print(lapply(probs, function(x) round(x, 4)))
  }

  cat("\n\n")
  cat("Parameter estimates")
  cat("\n\n")
  print(round(data.frame(coefs), 3))

  # to return by main function
  list(formula  = formula,
       coefs    = coefs,
       dtable   = d0,
       dfitted  = mutate(dpop, fitted = fitted) %>% select(-Freq),
       D        = as.matrix(D),
       vcv      = vcv,
       dimp     = dimp,
       hist     = hist)
}
