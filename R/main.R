#' A multi-purpose log linear model
#'
#' @description \code{loglinear} fits a log linear model with the EM-algoritm. The model
#' that allows for imputation of missing values, the inclusion of latent variables,
#' multiple systems estimation and/or regularization.
#'
#' @param formula formula of the model to be fitted. If the data contain individual records
#' the lhs of the formula should be left empty.
#' @param data a data frame with the data, either with individual records or in the form
#' of a contingency table. See "Details" for the specification of the variables.
#' @param latlevels a vector with the number of levels of the latent variable(s) in
#' the model. If present, the default \code{NULL} sets the levels of each latent variable to 2.
#' @param lambda a scalar denoting the size of the ridge penalty \eqn{\lambda\sum\beta^2}. Defaults to 0.
#' @param stderr determines whether the standard error are computed using the incomplete or the
#' complete data log-likelihood. The latter option is much faster when the model is complex, but
#' underestimates the standard errors if there are missings and/or latent
#' variables. If not, both are identical and the complete data log likelihood is used.
#' @param control a list with parameters controlling the fitting process:
#' \itemize{
#' \item \code{maxit} the maximum number of iterations of the EM algorithm. Defaults to 500.
#' \item \code{tol} The algorithm converges when \eqn{LL - LL_{old} < tol}. Defaults to 1e-5.
#' \item \code{seed} seed for obtaining the random starting values. Defaults to a random seed.
#' \item \code{stepsize} step size for the Newton-Raphson updates of the parameter estimates.
#' }
#' @return a list with the following components:
#' \item{formula}{formula used to fit the model.}
#' \item{coefs}{log linear parameter estimates, se's, and t- and p-values.}
#' \item{dtable}{data frame in the form of a contingency table used to fit the model.}
#' \item{dfitted}{data frame the fitted frequencies, including those of the structural zeros.}
#' \item{D}{the model matrix of the fitted model}
#' \item{vcv}{variance-covariance matrix of the parameter estimates.}
#' \item{dimp}{sparse matrix that specifies how the frequencies of the observed (incomplete) data are
#' to be distributed over the completed data.}
#' \item{hist}{iteration history.}
#' @details
#' The function \code{loglinear} interprets the role of the variables in the formula
#' by their names.
#'
#' The response variable in the rhs of the formula contains the observed frequencies if the is
#' in the form of a contingency table, e.g. \code{Freq ~ A + B}. If the data contains individual
#' records the lhs of the formula  should be left empty, e.g. \code{~ A + B}.
#'
#' Variables starting with upper case letters A-W are interpreted as population
#' registers, and these variables should be coded 0 (not observed in the register) and 1 (observed
#' in the register). If two or more are population registers are present in the model multiple systems
#' estimation will performed, and a population size estimate will be provided.
#'
#' Variables named X, Y and/or Z are interpreted as latent variables. If present in the model, the
#' \code{latlevels} specifies their respective number of levels.
#'
#' Variables that are neither latent nor population registers should start with a lower case letter.
#' Variables of class \code{numeric} will be treated as such, but they should have a limited number of
#' unique values.
#'
#' @importFrom stats complete.cases model.matrix na.omit optim pnorm xtabs rnorm rmultinom
#' glm poisson runif update terms coef
#' @importFrom dplyr mutate left_join %>% across all_of select filter everything filter rename
#' @importFrom Matrix Matrix Diagonal crossprod tcrossprod
#' @importFrom pracma hessian
#' @export

loglinear <- function(formula, data,
                      latlevels = NULL,
                      lambda  = 0,
                      stderr  = c("incomplete", "complete"),
                      control = list())
{

  cal    <- match.call()

  maxit     <- ifelse(exists("maxit", control), control$maxit, 500)
  tol       <- ifelse(exists("tol", control), control$tol, 1e-5)
  seed      <- ifelse(exists("seed", control), control$seed, sample(1e+5, 1))
  stepsize  <- ifelse(exists("stepsize", control), control$stepsize, 1)


  # get the terms object of the formula
  fterms <- terms(x = formula, data = data)

  # if present, expand rhs dot as in Freq ~ .
  formula <- formula(fterms)

  # if present, rename response variable to "Freq" in formula and data
  if (attr(fterms, "response") == 1)
  {
    response <- all.vars(formula)[1]
    pos      <- which(names(data) == response)
    names(data)[pos] <- "Freq"
    formula  <- update(formula, Freq  ~ .)
  }


  # get the names of the variables in the model
  vars   <- all.vars(formula)

  # select manifest variables and check for numeric variables
  data     <- data[, vars[!vars %in% LETTERS[24:26]]]
  varclass <- sapply(data, class)
  numvars  <- names(varclass)[varclass %in% c("integer", "numeric") &
                                !names(varclass) %in% c(LETTERS[1:26], "Freq")]

  # construct contingency table with factors,
  # convert numeric variables back to numeric,
  # and delete rows with missings and zero frequency
  if ("Freq" %in% vars){
    d0 <- as.data.frame(xtabs(Freq ~ ., data = data, addNA = TRUE)) %>%
      mutate(across(all_of(numvars), as.numeric)) %>%
      filter(!(!complete.cases(.) & Freq == 0))
  } else {
    d0 <- as.data.frame(xtabs(     ~ ., data = data, addNA = TRUE)) %>%
      mutate(across(all_of(numvars), as.numeric)) %>%
      filter(!(!complete.cases(.) & Freq == 0))
    # include response "Freq" in original formula
    formula  <- update(formula, Freq  ~ .)
  }

  # identify names and number of latent variables and lists in formula
  lat    <- vars[vars %in% LETTERS[24:26]]
  nlat   <- length(lat)
  lca    <- nlat > 0

  # identify population registers in formula
  lists  <- vars[substr(vars, 1, 1) %in% LETTERS[1:23]][-1]
  mse    <- length(lists) > 1

  # check if there are missing values in the data
  # sort on complete and incomplete cases
  impute <- sum(complete.cases(d0)) < nrow(d0)
  if (impute)
    d0 <- rbind(na.omit(d0), d0[!complete.cases(d0), ])


  #####################################
  # glm for ordinary log linear model #
  #####################################

  if (!mse & !lca & !impute)
    return(glm(formula = formula, poisson, data = d0))

  #####################################
  # imputation model                  #
  #####################################

  ret <- fit(formula = formula, d0 = d0, cal = cal,
             lambda = lambda, stderr = stderr, maxit = maxit,
             tol = tol, seed = seed, stepsize = stepsize,
             lists = lists, lat = lat, latlevels = latlevels, numvars = numvars,
             impute = impute, lca = lca, mse = mse)



  invisible(ret)

}

