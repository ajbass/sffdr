ffinemap <- function(flfdr, fpi0, CI = 0.95) {
  lfdr <- pmin(flfdr,1)
  BF <- (mean(fpi0) / (1-mean(fpi0))) * ((1-lfdr ) / lfdr)
  PIP <- BF / sum(BF)
  o <- order(-PIP)
  CS <- cumsum(PIP[o])[order(o)]
  credible_set <- CS <= CI
  if (sum(credible_set) == 0 & (max(CS) >= CI)) {
    credible_set[which.max(PIP)] <- TRUE
  }
  return(data.frame(bayes_factor = BF,
                    PIP = BF / sum(BF),
                    credible_set = credible_set))
}
