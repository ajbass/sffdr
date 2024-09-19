#' @title
#' Estimate the functional proportion of null tests
#'
#' @description
#' he function \code{\link{fpi0est}} estimates the functional proportion of null tests given
#' a set of informative variables.
#'
#' @details
#' This code extends the function from the fFDR package to handle
#' multiple informative variables and linkage disequilibrium.
#'
#' @param p A vector of p-values.
#' @param z A vector of informative variables
#' @param pi0_model Model formula corresponding to \code{z} for the functional proportion of truly null tests.
#' @param indep_snps A boolean vector (same size as p) specifying the set of independent tests. Default is NULL and all tests are treated independently.
#' @param lambda A vector of values between [0,1] to estimate the functional proportion of truly null tests.
#' @param method Either the "gam" (generalized additive model) or "glm" (generalized linear models) approach. Default is "gam".
#' @param maxit The maximum number of iterations for "glm" approach. Default is 1000.
#' @param pi0.method.control A user specified set of parameters for convergence for either "gam" or "glm". Default is NULL. See \code{\link{gam.control}} or \code{\link{glm.control}}.
#' @param \ldots Additional arguments passed to \code{\link{gam}} or \code{\link{glm}}.
#'
#' @return
#' A list of object type "fpi0" containing:
#' \item{fpi0}{A table containing the functional proportion of truly null tests.}
#' \item{tableLambda}{Functional proportion of null tests at the lambda values}
#' \item{MISE}{MISE values.}
#' \item{lambda.hat}{The chosen lambda value.}
#'
#' @examples
#' # import data
#' data(bmi)
#'
#' # separate main p-values and conditioning p-values
#' p <- sumstats$bmi
#' z <- as.matrix(sumstats[, -1])
#'
#' # apply pi0_model to create model
#' knots <- c(0.005, 0.01, 0.025, 0.05, 0.1)
#' fmod <- pi0_model(z, knots = knots)
#'
#' # Estimate functional pi0
#' fpi0_out <- fpi0est(p, z = fmod$zt, pi0_model = fmod$fmod)
#' fpi0 <- fpi0_out$fpi0
#'
#' # See relationship of BFP/cholesterol/triglycerides and fpi0
#' plot(fmod$zt$bfp, fpi0)
#' plot(fmod$zt$cholesterol, fpi0)
#' plot(fmod$zt$triglycerides, fpi0)
#'
#' @author Andrew J. Bass, David G. Robinson (author of original function)
#' @seealso \code{\link{sffdr}}, \code{\link{plot.sffdr}}
#' @keywords fpi0est
#' @aliases fpi0est
#' @importFrom stats density binomial dnorm ecdf fitted.values formula glm glm.control optimize pnorm predict qnorm quantile smooth.spline
#' @export
fpi0est <- function(p,
                    z,
                    pi0_model,
                    indep_snps = NULL,
                    lambda = seq(0.05, 0.9, 0.05),
                    method = "gam",
                    maxit = 1000,
                    pi0.method.control = NULL,
                    ...) {
  k <- phi.hat <- omega <- delta.sq <- chosen <- . <- NULL

  if (min(p) < 0 || max(p) > 1) {
    stop("P-values not in valid range")
  }
  if (is.null(pi0.method.control)) {
    if (method == "gam") {
      control <- gam::gam.control(epsilon = 1e-9, bf.epsilon = 1e-9, maxit= 200, bf.maxit = 200)
    } else {
      control <- glm.control(maxit = maxit)
    }
  }

  pi0_out <- p_in <- p
  fpi0_out <- rep(1, length(p))
  z.full <- z
  na.p <- is.na(p)
  na.z <- apply(z, 1, anyNA)
  rm_na <- !na.z & !na.p
  p <- p[rm_na]
  z <- z[rm_na,]
  indep_snps <- indep_snps[rm_na]
  if (!is.null(indep_snps)) {
    p.fit <- p[indep_snps]
    z.fit <- z[indep_snps,, drop = F]
  } else {
    p.fit = p
    z.fit = z
  }
  # Model matrix
  fm <- formula(paste("phi", paste(pi0_model, collapse = " ")))
  pi0hat_func <- function(lambda) {
    z.fit$phi <- as.numeric(p.fit >= lambda)
    fit <- NULL
    if (method == "glm") {
      fit <- suppressWarnings(glm(fm, data = z.fit,
                                  family = constrained.binomial(1 - lambda),
                                  control = control,
                                  maxit = maxit,
                                  ...))

    } else if (method == "gam") {
      fit <- gam::gam(fm,
                      family = constrained.binomial(1 - lambda),
                      data = z.fit,
                      control = control,
                      ...)

    }
    return(fit)
  }

  # Store as model and then predict in dplyr!
  df <- data.frame(p.value = p.fit)

  fpi0s <- df %>%
    dplyr::mutate(i = 1:dplyr::n()) %>%
    tidyr::crossing(lambda = lambda) %>%
    dplyr::select(-i) %>%
    dplyr::group_by(lambda) %>%
    do(fpi0 = pi0hat_func(.$lambda[1]))

  tt = sapply(1:length(lambda), FUN = function(x)  pmin(fitted.values(fpi0s$fpi0[[x]]) / (1 - lambda[x]), 1))
  get_fit <- data.frame(lambda = rep(lambda, each = length(p.fit)),
                        fpi0 = as.numeric(tt))

  ref <- get_fit %>%
    dplyr::ungroup() %>%
    filter(lambda == min(lambda)) %>%
    .$fpi0

  fpi0s2 <- get_fit %>%
    group_by(lambda) %>%
    mutate(k = optimize(function(k) mean((ref - k * (1 - ref) - fpi0) ^ 2),
                        interval = c(-1, 1))$minimum,
           phi.hat = ref - k * (1 - ref),
           fpi0 = pmin(fpi0, 1),
           phi.hat = pmin(phi.hat, 1))

  pi0.S <- qvalue::pi0est(p.fit, pi0.method = "bootstrap")$pi0
  stats <- fpi0s2 %>%
    group_by(lambda) %>%
    dplyr::summarize(omega = mean((fpi0 - phi.hat)^2),
                     delta.sq = (max(mean(fpi0) - pi0.S, 0))^2) %>%
    mutate(MISE = omega + delta.sq)

  lambda.hat = stats$lambda[which.min(stats$MISE)]
  fpi0s <- fpi0s %>% mutate(chosen = (lambda == lambda.hat))
  fpi0 <- fpi0s %>% filter(chosen == TRUE)

  # Large LD need to do chunk-wise
  msize <- 100000
  if (nrow(z) > msize) {
    fpi0.return <- rep(1, nrow(z))
    chunks <- ceiling(dim(z)[1] / msize)
      for (i in 1:chunks) {
        ind1 <- msize * (i-1) + 1
        ind2 <- pmin(msize * (i), nrow(z))
        fpi0.return[ind1:ind2] <- predict(fpi0$fpi0[[1]], newdata = z[ind1:ind2,],  type = "response") / (1 - fpi0$lambda)
      }
  } else {
    fpi0.return <- pmin(predict(fpi0$fpi0[[1]], newdata = z,  type = "response") / (1 - fpi0$lambda),1)
  }
  fpi0_out[rm_na] <- fpi0.return

  ret <- list(fpi0 = fpi0_out,
              tableLambda = fpi0s,
              MISE = stats,
              lambda = lambda.hat)
  class(ret) <- "fpi0"
  ret
}

# From fFDR package
constrained.binomial = function(maximum) {
  link <- structure(list(name = paste0("constrained.logit (0, ", maximum, ")"),
                         linkfun = function(mu) log((mu / maximum) / (1 - mu / maximum)),
                         linkinv = function(eta) maximum / (1 + exp(-eta)),
                         mu.eta = function(eta) maximum / (1 + exp(-eta)),
                         valideta = function(eta) TRUE), class = "link-glm")

  fam <- binomial(link)
  fam$validmu <- function(mu) all(mu > 0) && all(mu < maximum)
  fam$family <- paste0("constrained.binomial (0, ", maximum, ")")

  # for mgcv
  fam$d2link <- function(mu) { 1 / (1 - (mu / maximum)) ^ 2 - 1 / (mu / maximum) ^ 2 }
  fam$d3link <- function(mu) { 2 / (1 - (mu / maximum)) ^ 3 + 2 / (mu / maximum) ^ 3 }
  fam$d4link <- function(mu) { 6 / (1 - (mu / maximum)) ^ 4 - 6 / (mu / maximum) ^ 4 }
  fam$dvar <- function(mu) rep.int(1, length(mu))
  fam$d3var <- fam$d2var <- function(mu) rep.int(0, length(mu))

  # new addition to initialization: mu cannot be greater than maximum
  new.line <- substitute(mustart <- mustart * maximum, list(maximum = maximum))
  # convert call to expressiion so it be initialized by eval
  fam$initialize <- as.expression(c(fam$initialize, new.line))

  fam
}
