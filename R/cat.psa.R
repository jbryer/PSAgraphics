#' Compare balance graphically of a categorical covariate as part of a PSA
#'
#' Given predefined strata and two level treatment for a categorical covariate
#' from a propensity score analysis, \code{cat.psa} draws pairs of side by side
#' barplots corresponding to control and treatment for each stratum.
#'
#' Pairs of bars are graphed side by side so that comparisons may be made
#' within each stratum and across strata.  If \code{balance} is \code{TRUE},
#' then the histogram represents an ad hoc balance measure of the given strata
#' as compared to randomly generated strata.  The p-values provided on the
#' bargraph are bootstrapped in a standard fashion via randomly generated
#' treatment divisions within given strata.  For continuous covariates use
#' \code{box.psa}.
#'
#' @param categorical Vector or N X 3 dataframe or matrix. If a vector, then
#' represents a categorical covariate that is being balanced within strata in a
#' PSA. If \code{categorical} has three columns, then the second and third are
#' assumed to be the \code{treatment} and \code{strata} respectively.  Missing
#' values are not allowed. May be factor or numeric.
#' @param treatment Binary vector or factor of same length as \code{continuous}
#' representing the two treatments.
#' @param strata A vector or factor of same length as \code{continuous}
#' indicating the derived strata from estimated propensity scores. Strata are
#' ordered lexicographically in plot.
#' @param catnames List of names in order of the categories; used in the plot
#' legend. Default is \code{1:n}.
#' @param catcol List of colors used for the categories, default is
#' \code{terrain.colors}.
#' @param width Controls width of bars, default = 0.25.
#' @param barlab Binary list of single \code{treatment} character labels for
#' the bars, default is \code{c("A", "B")}.  These are defined in a legend by
#' \code{barnames}.
#' @param barnames Binary list of treatment names used in the legend; by
#' default names are taken from \code{treatment}.
#' @param rtmar Numeric.  Governs size of right margin allocated for legend.
#' Default = 1.5
#' @param balance Logical.  If \code{TRUE} a call is made to functions
#' \code{bal.cs.psa} and \code{bal.cws.psa}. The former provides a reference
#' histogram and ad hoc balance statistic, the second provides bootstrapped
#' p-values for the two-way table formed in each statum.  Default is
#' \code{FALSE}.
#' @param B Numeric; passed to \code{bal.cs.psa} governing size of reference
#' histogram generated.  Default is 100.
#' @param tbl Logical; if \code{TRUE}, then a matrix of the proportions used in
#' the creation of the bargraph is returned.
#' @param cex.leg Numeric; value of \code{cex} (governing font size) passed to
#' legend.  Default = 1.
#' @param \dots Other graphical parameters passed to plot.
#' @return If \code{tbl} is \code{TRUE}, then a matrix is returned containing
#' the proportions of each category, and in each treatment level and stratum
#' that were used to draw the bargraph.
#' @author James E. Helmreich \email{ James.Helmreich@@Marist.edu}
#'
#' Robert M. Pruzek \email{RMPruzek@@yahoo.com}
#' @seealso \code{bal.cs.psa}, \code{bal.cws.psa}, \code{box.psa}
#' @keywords hplot htest
#' @examples
#'
#' categorical<-sample(1:7,1000,replace=TRUE)
#' treatment<-sample(c(0,1),1000,replace=TRUE)
#' strata<-sample(5,1000,replace=TRUE)
#' cat.psa(categorical,treatment,strata)
#'
#' data(lindner)
#' attach(lindner)
#' lindner.ps <- glm(abcix ~ stent + height + female +
#'       diabetic + acutemi + ejecfrac + ves1proc,
#'       data = lindner, family = binomial)
#' ps<-lindner.ps$fitted
#' lindner.s5 <- as.numeric(cut(ps, quantile(ps, seq(0, 1, 1/5)),
#'       include.lowest = TRUE, labels = FALSE))
#' cat.psa(stent, abcix, lindner.s5, xlab = "stent")
#'
#' lindner.s10 <- as.numeric(cut(ps, quantile(ps, seq(0, 1, 1/10)),
#'       include.lowest = TRUE, labels = FALSE))
#' cat.psa(ves1proc,abcix, lindner.s10, balance = TRUE, xlab = "ves1proc")
#'
#' #Using a rpart tree for strata
#' library(rpart)
#' lindner.rpart<-rpart(abcix ~ stent + height + female + diabetic +
#'       acutemi + ejecfrac + ves1proc, data=lindner, method="class")
#' lindner.tree<-factor(lindner.rpart$where, labels = 1:6)
#' cat.psa(stent, abcix, lindner.tree, xlab = "stent")
#' cat.psa(ves1proc, abcix, lindner.tree, xlab = "ves1proc")
#'
#' @export cat.psa
#' @importFrom rpart rpart
cat.psa <- function(categorical,
                    treatment = NULL,
                    strata = NULL,
                    catnames = NULL,
                    catcol = "terrain.colors",
                    width = .25,
                    barlab = c("A", "B"),
                    barnames = NULL,
                    rtmar = 1.5,
                    balance = FALSE,
                    B = 1000,
                    tbl = TRUE,
                    cex.leg = 1,
                    ...)
{
  #categorical should be a categorical variable with levels 1:n
  #treatment should be the binary treatment variable from a PSA
  #strata should be the numerical strata that subjects belong to (from a logistic PSA, usually 1:5)
  #width will be width of bars

  #If "balance" is TRUE, then two different heuristic bootstrap distributions
  #and balance statistics are found. First a histogram of balance measures is
  #given from randomly generated strata, and compared with the measure for the
  #actual strata. Second, within strata, the covariate balance between treatments
  #is bootstrapped and measured (in a standard way[I think - need to talk to Bob]),
  #and then a "p-value" is calculated and given in the bar chart at the bottom of
  #each stratum in red.

  #Ok, more on the balance statistics calculated:  For the first, randomly
  #generated strata are generated. In a given stratum, proportions in each
  #treatment are calculated, and then the absolute difference is between
  #these proportions in each category is found, and summed across categories.
  #This value is summed across strata for the final balance statistic.  The
  #same statistic from the original is calculated and compared to the
  #randomly generated values via a histogram and ranking.  In the second
  #measure, the strata are considered fixed, and the covariate values in a
  #given statum are randomly sampled without replacement and assigned to
  #the two treatments.  Absolute differences in size of categories are
  #found, and summed across categories, to be compared with the true
  #covariate distribution.

  #If "categorical" has three columns, treat as c, t, s.
  if (dim(as.data.frame(categorical))[2] == 3) {
    treatment   <- categorical[, 2]
    strata      <-
      categorical[, 3]
    categorical <-
      categorical[, 1]
  }
  cat.f <- as.factor(categorical)
  cat.levels <- levels(cat.f)

  table.cts <- table(categorical, treatment, strata)
  cat.dim <- dim(table.cts)[1]
  strata.dim <- dim(table.cts)[3]
  if (is.null(catnames) &
      !is.factor(categorical)) {
    catnames <- sort(unique(categorical))
  }
  if (is.null(catnames)) {
    if (is.factor(categorical)) {
      catnames <- levels(categorical)
    } else{
      catnames <- c(1:cat.dim)
    }
  }
  if (catcol[1] == "terrain.colors") {
    catcol <- terrain.colors(cat.dim)
  }

  #Creating the x coordinates of the rectangles
  minlim <- 1:strata.dim - width
  midlim <- 1:strata.dim
  maxlim <- 1:strata.dim + width

  x.left <-
    as.vector(matrix(rep(minlim, cat.dim), nrow = cat.dim, byrow = TRUE))
  x.mid <-
    as.vector(matrix(rep(midlim, cat.dim), nrow = cat.dim, byrow = TRUE))
  x.right <-
    as.vector(matrix(rep(maxlim, cat.dim), nrow = cat.dim, byrow = TRUE))

  #x.l/x.r have left x coordinates of all 0/1 rectangles for each stratum
  x.l <- c(x.left, x.mid)
  x.r <- c(x.mid, x.right)

  #Want to create the y coordinates of the rectangles.
  #Finding the numbers of 0/1 by stratum in sum.strata.01

  sum.strata.01 <- NULL
  for (j in 1:strata.dim)
    sum.strata.01 <- rbind(sum.strata.01, apply(table.cts[, , j], 2, sum))

  #Calculating the vectors of bottom and top y values for the rectangles in each treatment
  y.b.0 <- NULL
  y.t.0 <- NULL
  y.b.1 <- NULL
  y.t.1 <- NULL
  prop.zo <- NULL
  for (j in 1:strata.dim) {
    # prop.i is a vector of the proportion of 0/1 by categorical level
    # level for stratum j
    prop.0 <- table.cts[, 1, j] / sum.strata.01[j, 1]
    prop.1 <- table.cts[, 2, j] / sum.strata.01[j, 2]
    C <- prop.0
    T <- prop.1
    prop.zo <- cbind(prop.zo, C)
    prop.zo <- cbind(prop.zo, T)
    # prop.i.sum is a vector of the cumulative proportions
    prop.0.sum <- NULL
    prop.1.sum <- NULL
    for (i in 1:cat.dim)
      prop.0.sum <- cbind(prop.0.sum, sum(prop.0[1:i]))
    for (i in 1:cat.dim)
      prop.1.sum <- cbind(prop.1.sum, sum(prop.1[1:i]))
    # The bottom y's should start at 0 and end at the penultimate entry of prop.i.sum.
    # The top y's should start with first entry of prop.i.sum and end at 1, which is last entry of prop.i.sum
    y.b.0 <- c(y.b.0, 0, prop.0.sum[-cat.dim])
    y.t.0 <- c(y.t.0, prop.0.sum)
    y.b.1 <- c(y.b.1, 0, prop.1.sum[-cat.dim])
    y.t.1 <- c(y.t.1, prop.1.sum)
  }
  y.b <- c(y.b.0, y.b.1)
  y.t <- c(y.t.0, y.t.1)

  if (balance) {
    bal.cs <- bal.cs.psa(categorical, treatment, strata, B = B)
    cat("Histogram of Random Strata Balance. Press <enter> for next chart...")
    readline()
    bal.fe <- bal.fe.psa(categorical, treatment, strata)
  }

  #Creating Graph
  xlimits <- c(.5, strata.dim + rtmar)

  plot(
    c(.5, 2),
    type = "n",
    log = "",
    axes = FALSE,
    xlim = xlimits,
    ylim = c(-.05, 1.05),
    ...
  )
  rect(x.l, y.b, x.r, y.t, col = catcol, lty = "solid")
  axis(1, at = 1:strata.dim, labels = sort(unique(strata)))
  axis(2, at = c(0, .5, 1))

  if (is.null(barnames)) {
    barnames <- unlist(dimnames(table.cts)[2])
  }

  for (i in 1:strata.dim) {
    text  (i - .125, 0, barlab[1], cex = .7, pos = 1)
  }
  for (i in 1:strata.dim) {
    text  (i + .125, 0, barlab[2], cex = .7, pos = 1)
  }
  table.ts <- table(treatment, strata)
  for (i in 1:strata.dim) {
    text (
      i - width / 2,
      1,
      table.ts[1, i],
      cex = .7,
      pos = 3,
      col = 4
    )
  }
  for (i in 1:strata.dim) {
    text (
      i + width / 2,
      1,
      table.ts[2, i],
      cex = .7,
      pos = 3,
      col = 4
    )
  }
  if (balance) {
    lnd <- c("No. Obs.", "F.E. p-val", barnames)
    ppch <- c("#", "#", barlab)
    ccol <- c(4, 2, 1, 1)
    for (i in 1:strata.dim) {
      text  (i,
             -.067,
             round(bal.fe[i], 2),
             col = 2,
             cex = .7)
    }
  } else{
    lnd <- c("No. Obs.", barnames)
    ppch <- c("#", barlab)
    ccol <- c(4, 1, 1)
  }
  legend(
    strata.dim + .4,
    .3,
    legend = lnd,
    pch = ppch,
    col = ccol,
    bty = "n",
    cex = cex.leg
  )
  legend(
    strata.dim + .4,
    .97,
    c("Levels", catnames),
    pch = c(0, rep(15, cat.dim)),
    col = c(0, catcol),
    bty = "n",
    cex = cex.leg
  )

  ##Output a table of percents for catagorical levels within each treatment, stratum.
  if (tbl) {
    Levels <- NULL
    numprop <- NULL
    tns <- NULL
    t.names <- unlist(dimnames(table.cts)[2])
    s.names <- unlist(dimnames(table.cts)[3])
    for (i in 1:strata.dim) {
      for (j in 1:2) {
        tns <- c(tns, paste(t.names[j], ":", s.names[i], sep = ""))
      }
    }
    strata.names <- unlist(dimnames(table.cts)[3])
    treatment.stratum.proportions <- round(prop.zo, 3)
    colnames(treatment.stratum.proportions) <- tns
    out <- list(treatment.stratum.proportions)
    names(out) <-
      c(paste("treatment", ":", "stratum", ".proportions", sep = ""))
    return(out)
  }
}
