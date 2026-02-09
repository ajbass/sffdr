#' @title Plotting function for sffdr object
#' @description
#'  Graphical display of the sffdr object
#'
#' @param x A sffdr object.
#' @param rng Significance region to show. Optional.
#' @param \ldots Additional arguments. Currently unused.
#'
#' @return
#' Nothing of interest.
#'
#' @examples
#' # import data
#' data(bmi)
#'
#' # Separate main p-values and conditioning p-values
#' p <- sumstats$bmi
#' z <- as.matrix(sumstats[, -1])
#'
#' # Apply pi0_model to create model
#' fmod <- pi0_model(z)
#'
#' # Estimate functional pi0
#' fpi0_out <- fpi0est(p, z = fmod$zt, pi0_model = fmod$fmod)
#' fpi0 <- fpi0_out$fpi0
#'
#' # Apply sffdr
#' sffdr_out <- sffdr(p, fpi0)
#'
#' # Plot significance results
#' plot(sffdr_out, rng = c(0, 5e-4))
#'
#' @author Andrew J. Bass
#' @seealso \code{\link{sffdr}}
#' @keywords plot
#' @aliases plot, plot.sffdr
#' @import patchwork ggplot2
#' @export
plot.sffdr <- function(x, rng = c(0, 5e-8), ...) {
  z <- pvalue <- fpvalue <- threshold <- sig <- fqvalue <- NULL
  rm_na <- !is.na(x$pvalues)
  pvalues <- x$pvalues[rm_na]
  qvalues <- x$fqvalues[rm_na]
  fpvalues <- x$fpvalues[rm_na]

  fp.ord <- fpvalues[order(fpvalues)]
  if (min(fp.ord) > rng[2]) {
    rng <- c(min(fp.ord), quantile(fp.ord, 1e-2))
  }

  p.ord <- pvalues[fp.ord]

  pi0 <- x$fpi0

  pi0.df <- data.frame(
    z = rank(pi0, ties.method = "random")[which(
      fpvalues >= rng[1] & fpvalues <= rng[2]
    )] /
      length(pi0),
    pi0 = pi0[which(fpvalues >= rng[1] & fpvalues <= rng[2])]
  )

  # surrogate versus pi0
  p1 <- ggplot(pi0.df, aes(x = z, y = pi0)) +
    #geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_log10() +
    xlab("surrogate variable") +
    ylab("prior probability of null") +
    scale_color_brewer(palette = "Set1")

  # functional p versus raw p
  p2 <- ggplot(
    data.frame(
      pvalue = pvalues[which(fpvalues >= rng[1] & fpvalues <= rng[2])],
      fpvalue = fpvalues[which(fpvalues >= rng[1] & fpvalues <= rng[2])]
    ),
    aes(x = pvalue, y = fpvalue)
  ) +
    ylab("functional p-value") +
    scale_x_log10() +
    scale_y_log10() +
    xlab("p-value") +
    geom_point() +
    theme_bw()

  # Functional p versus discoveries
  df1 <- data.frame(
    threshold = fp.ord[fp.ord >= rng[1] & fp.ord <= rng[2]],
    sig = (1 + sum(fp.ord < rng[1])):sum(fp.ord <= rng[2])
  )
  df2 <- data.frame(
    threshold = sort(pvalues[pvalues >= rng[1] & pvalues <= rng[2]]),
    sig = (1 + sum(pvalues < rng[1])):sum(pvalues <= rng[2])
  )
  if (max(df2$threshold) > max(df1$threshold)) {
    df1 <- rbind(
      df1,
      data.frame(threshold = max(df2$threshold), sig = max(df1$sig))
    )
  } else {
    df2 <- rbind(
      df2,
      data.frame(threshold = max(df1$threshold), sig = max(df2$sig))
    )
  }
  if (min(df2$threshold) < min(df1$threshold)) {
    df1 <- rbind(
      df1,
      data.frame(threshold = min(df2$threshold), sig = min(df1$sig))
    )
  } else {
    df2 <- rbind(
      df2,
      data.frame(threshold = min(df1$threshold), sig = min(df2$sig))
    )
  }

  p3 <- ggplot(df1, aes(x = threshold, y = sig, linetype = "Functional P")) +
    xlab("threshold") +
    ylab("significant tests") +
    geom_line() +
    geom_line(data = df2, aes(x = threshold, y = sig, linetype = "Raw P")) +
    theme_bw() +
    scale_linetype_manual("", values = c("solid", "dashed")) +
    theme(
      strip.text.x = element_blank(),
      strip.background = element_rect(fill = NA),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = c(.25, .9),
      legend.text = element_text(size = 7)
    )

  # functional p versus q
  p4 <- ggplot(
    data.frame(
      fpvalue = fpvalues[which(fpvalues >= rng[1] & fpvalues <= rng[2])],
      fqvalue = qvalues[which(fpvalues >= rng[1] & fpvalues <= rng[2])]
    ),
    aes(x = fpvalue, y = fqvalue)
  ) +
    xlab("functional p-value") +
    scale_x_log10() +
    scale_y_log10() +
    ylab("functional q-value") + #geom_point() +
    theme_bw() +
    geom_line() +
    theme_bw() +
    scale_x_log10()

  (p1 + p2) /
    (p3 + p4) +
    plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")")
}
