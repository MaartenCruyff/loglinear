#' New Zealand population data of 2013 to 2020.
#'
#' Contingency table containing population registers and
#' ethnic classifications over the years 2013 to 2020. The frequencies are
#' rounded to give effect to the confidentiality provisions of the
#' Statistics Act of 1975.
#'
#' @format Data frame with 7455 rows and 10 variables.
#'
#' \describe{
#'   \item{A}{Register A (0 = case not observed in A, 1  = case observed in A)}
#'   \item{B}{Register B (0 = case not observed in B, 1  = case observed in B)}
#'   \item{C}{Register D (0 = case not observed in C, 1  = case observed in C)}
#'   \item{B}{Register D (0 = case not observed in D, 1  = case observed in D)}
#'   \item{a}{ethnic classifications in register A (levels 1, 2, 3, 4)}
#'   \item{b}{ethnic classifications in register B (levels 1, 2, 3, 4)}
#'   \item{c}{ethnic classifications in register C (levels 1, 2, 3, 4)}
#'   \item{d}{ethnic classifications in register D (levels 1, 2, 3, 4)}
#'   \item{y}{the year the observation was made in}
#'   \item{Freq}{The observed frequencies}
#' }
#'
#' @details The variables \code{a, b, c, d} are partially misclassified and
#' contain missing values.
#'
#' @source Stats New Zealand.
"nzy"

#' Updated version of \code{\link{nzy}} with variables A-D as factors.
#'
#' @format Data frame with 7028 rows and 10 variables.
#'
#' @source Stats New Zealand.
"nzy2"


#' New Zealand population data of 2013.
#'
#' An extraction of the year 2013 from the data \code{\link{nzy}}. The
#' variable \code{Y} is not included in the extraction.
#'
#' @format Data frame with 923 rows and 9 (categorical) variables.
#'
#' @details The \code{class} of variables with population registers and ethnic
#' classifications has been changed form \code{numeric} to \code{factor}. See
#' \code{\link{nzy}} for the variable specifications.
#'
#' @source Stats New Zealand.
"nz13"
