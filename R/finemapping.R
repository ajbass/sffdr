#' @title
#' Functional fine mapping
#'
#' @description
#' Perform functional fine mapping with a set of functional local FDRs in a region of interest (assuming a single causal variant).
#'
#' @param flfdr A vector of functional local FDRs of a region of interest
#' @param fpi0 An estimate of the function proportion of null tests
#'
#' @return
#' A list of object type "sffdr" containing:
#' \item{BF}{The functional local Bayes' factors.}
#' \item{PP}{Posterior probability of a SNP being causal.}
#'
#' @export
ffinemap <- function(flfdr, fpi0) {
  lfdr <- pmin(flfdr,1)
  BF <- (mean(fpi0) / (1 - mean(fpi0))) * ((1-lfdr ) / lfdr)
  PP <- BF / sum(BF)
#  o <- order(-PIP)
 # CS <- cumsum(PIP[o])[order(o)]
  # subset <- CS <= credible_set
  # if (sum(credible_set) == 0 & (max(CS) >= credible_set)) {
  #   subset[which.max(PIP)] <- TRUE
  # }
  return(data.frame(BF = BF,
                    PP = PP))
}
