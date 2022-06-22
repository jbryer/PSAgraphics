#' Compare balance graphically of a continuous covariate as part of a PSA
#'
#' Given predefined strata and two level treatment for a continuous covariate
#' from a propensity score analysis, \code{box.psa} draws pairs of side by side
#' boxplots corresponding to control and treatment for each stratum.
#'
#' Draws a pair of side by side boxplots for each stratum of a propensity score
#' analysis.  This allows visual comparisons within strata of the distribution
#' of the given continuous covariate, and comparisons between strata as well.
#' The number of observations in each boxplot are given below each box, and the
#' means of paired treatment and control groups are connected.
#'
#' @param continuous Vector or N X 3 dataframe or matrix. If a vector, then
#' represents the quantitative covariate that is being balanced within strata
#' in a PSA. If \code{continuous} has three columns, then the second and third
#' are assumed to be the \code{treatment} and \code{strata} respectively.
#' Missing values are not allowed.
#' @param treatment Binary vector of same length as \code{continuous}
#' representing the two treatments; can be a character vector or factor.
#' @param strata A vector or factor of same length as \code{continuous}
#' indicating the derived strata from estimated propensity scores. May be
#' numeric or character vector, or factor. Strata are ordered lexicographically
#' in plot.
#' @param boxwex Numeric; controls width of boxes. Default = 0.17
#' @param offset Numeric; controls distance between the two boxes in each
#' stratum. Default = 0.17
#' @param col Default = \code{c("yellow", "orange", "black", "red",
#' "darkorange3")}. Color vector for the control boxes, treatment boxes, and
#' line connecting their means.
#' @param xlab Label for the x-axis; default = \code{"Stratum"}.  Other
#' standard labels may be used as well.
#' @param legend.xy Binary vector giving coordinates of the legend. By default
#' the legend is placed to the top left.
#' @param legend.labels Vector of labels for the legend; default is essentially
#' \code{c("Treatment (first)", "Treatment (second)", "Treatment Means
#' Compared", "KS p-values", "Strata-Treatment Size")} where treatment names
#' are taken from \code{treatment}.  Vector has four elements if \code{balance
#' = FALSE}, ommitting "KS p-values".
#' @param pts Logical; if \code{TRUE} then (jittered) points are added on top
#' of the boxplots.
#' @param balance Logical; if \code{TRUE} then \code{bal.ms.psa} provides a
#' histogram of a permutation distribution and reference statstic to assess
#' balance across strata; \code{bal.ks.psa} adds p-values to the graph derived
#' from 2-sample Kologmorov-Smirnov tests of equivalence of control/treatment
#' distributions within each stratum.
#' @param trim If \code{balance=TRUE}, defines fraction (0 to 0.5) of
#' observations to be trimmed from each end of stratum-treatment level before
#' the mean is computed. See \code{\link{mean}}, \code{\link{bal.ms.psa}}.
#' @param B Passed to \code{bal.ms.psa} if necessary, determines number of
#' randomly generated comparison statistics.  Default =1000.
#' @param \dots Other graphical parameters passed to \code{boxplot}.
#' @author James E. Helmreich \email{James.Helmreich@@Marist.edu}
#'
#' Robert M. Pruzek \email{RMPruzek@@yahoo.com}
#' @seealso \code{bal.ks.psa}, \code{bal.ms.psa}, \code{cat.psa}
#' @keywords hplot htest
#' @examples
#'
#' continuous<-rnorm(1000)
#' treatment<-sample(c(0,1),1000,replace=TRUE)
#' strata<-sample(5,1000,replace=TRUE)
#' box.psa(continuous, treatment, strata)
#'
#' data(lindner)
#' attach(lindner)
#' lindner.ps <- glm(abcix ~ stent + height + female +
#'       diabetic + acutemi + ejecfrac + ves1proc,
#'       data = lindner, family = binomial)
#' ps<-lindner.ps$fitted
#' lindner.s5 <- as.numeric(cut(ps, quantile(ps, seq(0, 1, 1/5)),
#'       include.lowest = TRUE, labels = FALSE))
#' box.psa(ejecfrac, abcix, lindner.s5, xlab = "ejecfrac",
#'       legend.xy = c(3.5,110))
#'
#' lindner.s10 <- as.numeric(cut(ps, quantile(ps, seq(0, 1, 1/5)),
#'       include.lowest = TRUE, labels = FALSE))
#' box.psa(height, abcix, lindner.s10, xlab="height",
#'       boxwex = .15, offset = .15, legend.xy = c(2,130), balance = TRUE)
#'
#' @export box.psa
box.psa <- function(continuous,
                    treatment = NULL,
                    strata = NULL,
                    boxwex = 0.17,
                    offset = 0.17,
                    col = c("yellow", "orange", "black", "red", "darkorange3"),
                    xlab = "Stratum",
                    legend.xy = NULL,
                    legend.labels = NULL,
                    pts = TRUE,
                    balance = FALSE,
                    trim = 0,
                    B = 1000,
                    ...) {
  #Plots side by side boxplot of C/T for each strata.  Means may be connected if desired.
  #Point clouds for each boxplot may be plotted.  The nonparametric Kolmogorov-Smirnov test
  #may be used to find p-values for the test of equivalence of distribution between C/T in
  #each stratum.  In legend.lables, note that treatment levels A and B are actually taken from the treatment factor.

  #If "continuous" has three columns, treat as m, t, s.
  if (dim(as.data.frame(continuous))[2] == 3) {
    treatment   <- continuous[, 2]
    strata     <-
      continuous[, 3]
    continuous <-
      continuous[, 1]
  }
  cts <- data.frame(continuous, treatment, strata)

  #Sort the unique treatment levels.  To be used in legend as well.
  sut <- sort(unique(treatment))

  #Getting the legend labels sorted out.
  leg.test <- is.null(legend.labels)
  if (balance) {
    if (leg.test) {
      legend.labels <-
        c(
          paste("Treatment", sut[1]),
          paste("Treatment", sut[2]),
          "Treatment Means Compared",
          "KS p-values",
          "Strata-Treatment Size"
        )
    }
  } else{
    legend.labels <- legend.labels
  }
  if (!balance) {
    if (leg.test) {
      legend.labels <-
        c(
          paste("Treatment", sut[1]),
          paste("Treatment", sut[2]),
          "Treatment Means Compared",
          "Strata-Treatment Size"
        )
    }
  } else{
    legend.labels <- legend.labels
  }

  cs.0 <- subset(cts, treatment == sut[1], select = c(continuous, strata))
  cs.1 <- subset(cts, treatment == sut[2], select = c(continuous, strata))
  size <- table(treatment, strata)

  if (balance) {
    bal.ms <- bal.ms.psa(continuous, treatment, strata, B, trim = trim)
    cat("Press <enter> for bar chart...")
    readline()
  }

  s.d <- dim(table(strata))
  rgy <- range(continuous)
  ht <- rgy[2] - rgy[1]
  if (is.null(legend.xy)) {
    legend.xy <- c(1, rgy[2] + .22 * ht)
  }

  boxplot(
    continuous ~ strata,
    boxwex = boxwex,
    at = 1:s.d - offset,
    axes = FALSE,
    ylim = c(rgy[1] - .13 * ht, rgy[2] + .17 * ht),
    subset = treatment == sut[1],
    col = col[1],
    xlab = xlab,
    ...
  )

  boxplot(
    continuous ~ strata,
    add = TRUE,
    axes = FALSE,
    boxwex = boxwex,
    at = 1:s.d + offset,
    subset = treatment == sut[2],
    col = col[2],
    ...
  )

  #Add point clouds to boxplots
  if (pts) {
    for (i in 1:s.d) {
      cs.0.strat <- subset(cs.0, cs.0[, 2] == sort(unique(strata))[i])
      points(
        jitter(rep(i, length(cs.0.strat[, 1])) - offset, amount = .25 * boxwex),
        jitter(cs.0.strat[, 1], amount = .10 * boxwex),
        col = "yellow3",
        cex = .8
      )
    }
    for (i in 1:s.d) {
      cs.1.strat <- subset(cs.1, cs.1[, 2] == sort(unique(strata))[i])
      points(
        jitter(rep(i, length(cs.1.strat[, 1])), amount = .25 * boxwex) + offset,
        jitter(cs.1.strat[, 1], amount = .10 * boxwex),
        col = "orange3",
        cex = .8
      )
    }
  }

  #Add Kolgomorov-Smirnov p-values
  if (balance) {
    p.ks <- round(bal.ks.psa(continuous, treatment, strata), 2)
    for (i in 1:s.d) {
      text  (i,
             rgy[1] - (rgy[2] - rgy[1]) / 10,
             p.ks[i],
             col = col[4],
             cex = .7)
    }
  }

  axis(1, 1:s.d, labels = sort(unique(strata)))
  axis(2, at = NULL)

  #Size of strata
  for (i in 1:s.d) {
    text (i - offset,
          rgy[1] - (rgy[2] - rgy[1]) / 7.4,
          size[1, i],
          col = col[5],
          cex = .7)
  }
  for (i in 1:s.d) {
    text (i + offset,
          rgy[1] - (rgy[2] - rgy[1]) / 7.4,
          size[2, i],
          col = col[5],
          cex = .7)
  }

  strata.factor <- as.factor(strata)
  levels(strata.factor) <- 1:s.d
  strata.vector <- as.vector(strata.factor)

  #Create lines matching means
  x.t <-
    as.vector(rbind(
      1:s.d - offset + boxwex / 2,
      1:s.d + offset - boxwex / 2,
      rep("NA", s.d)
    ))
  last <- length(x.t)
  x.coord <- x.t[-last]
  y.t <-
    aggregate(
      cts[, 1],
      by = list(treatment = treatment, strata = strata),
      FUN = mean,
      trim = trim
    )
  y.tt <- NULL
  for (i in 1:s.d) {
    y.tt <- c(y.tt, y.t[(2 * i - 1), 3], y.t[2 * i, 3], "NA")
  }
  y.coord <- y.tt[-last]
  lines(
    x.coord,
    y.coord,
    col = col[3],
    lty = 1,
    type = 'o',
    lwd = 2
  )

  #Legend
  if (balance) {
    legend(
      x = legend.xy[1],
      y = legend.xy[2],
      legend = legend.labels,
      col = col,
      pch = c(15, 15, -1, 15, 15),
      lty = c(-1, -1, 1, -1, -1),
      lwd = c(2, 2, 2, 2),
      cex = .8
    )
  }
  else{
    legend(
      x = legend.xy[1],
      y = legend.xy[2],
      legend = legend.labels,
      col = col[c(1, 2, 3, 5)],
      pch = c(15, 15, -1, 15),
      lty = c(-1, -1, 1, -1),
      lwd = c(2, 2, 2, 2),
      cex = .8
    )
  }
}
