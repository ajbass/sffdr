#' Subset of p-values from the UK Biobank
#'
#' A dataset containing a subset of p-values from the UK Biobank.
#'
#' @format A data frame with 10,000 rows and 4 columns:
#' \describe{
#'   \item{bmi}{Body mass index}
#'   \item{bfp}{Body fat percentage}
#'   \item{cholesterol}{Cholesterol}
#'   \item{triglycerides}{Triglycerides}
#' }
#'
#' @examples
#' # Import data
#' data(bmi)
#'
#' # Separate main p-values and informative p-values
#' p <- sumstats$bmi
#' z <- as.matrix(sumstats[, -1])
#'
#' @keywords datasets
#' @name bmi
#' @aliases sumstats
NULL
