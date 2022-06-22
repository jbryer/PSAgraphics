#' Supports Multiple Methods for Defining and Visualizing (PS) Strata
#'
#' Given propensity scores, allows strata to be directly user defined, possibly
#' to: equalize sizes of strata, equalize the ranges of propensity scores, or
#' to specify cut points on the unit interval. Once strata are created, a
#' simple graphic is generated to visualize or judge strata for overlap and
#' appropriateness. If a regression tree has been used, propensity scores are
#' defined for each leaf of the tree.
#'
#'
#' @param treatment Binary vector or factor defining the two treatments
#' @param propensity Vector of same length as \code{treatment} containing
#' estimated propensity scores.
#' @param strata Either a vector of same length as \code{treatment} of
#' predefined stratum number, or one integer \code{n} used to assign rows to
#' \code{n} strata \code{propensity} scores, each of approximately the same
#' number of cases.  If relatively few unique propensity scores have been
#' defined (as from a classification tree) then the logical \code{tree} should
#' be set equal to \code{TRUE}.
#' @param int Either a number \code{m} used to divide \code{[0,1]} into
#' \code{m} equal length subintervals, or a vector containing cut points
#' between 0 and 1 that define subintervals (perhaps as suggested by
#' loess.psa).  In either case the subintervals define strata, for which sizes
#' can differ.
#' @param tree Logical, default \code{FALSE}.  If there are few unique
#' propensity scores, say from a recursively partitioned tree, then \code{TRUE}
#' forces strata to be defined by the unique propensity scores.
#' @param minsize Smallest allowable stratum-treatment size.  If violated, rows
#' in the stratum are removed.  User may wish to redefine strata.
#' @param graphic Logical, default \code{TRUE}.  If set to \code{FALSE} the
#' graphic is not provided.
#' @param colors 2-ary color vector.  Sets the colors of the points in the
#' graphic.  Default = \code{c("blue", "orange")}
#' @param xlab Label for the x axis; default = \code{"Estimated Propensity
#' Scores with Random Heights"}.
#' @param pch 2-ary vector; determines the shape of points in the graphic.
#' Default = \code{c(16, 16)}.
#' @return \item{Original.Strata}{Table of strata-treatment sizes before
#' \code{minsize} evaluation.} \item{Used.Strata}{Table of strata-treatment
#' sizes after \code{minsize} evaluation.} \item{strata}{Vector of the same
#' length as \code{treatment}, indicating either the strata input by user or
#' those created by the function.}
#' @author James E. Helmreich \email{ James.Helmreich@@Marist.edu}
#'
#' Robert M. Pruzek \email{RMPruzek@@yahoo.com}
#'
#' KuangNan Xiong \email{harryxkn@@yahoo.com}
#' @seealso \code{\link{cv.bal.psa}}, \code{\link{loess.psa}}
#' @keywords hplot
#' @examples
#'
#' data(lindner)
#' attach(lindner)
#' lindner.ps <- glm(abcix ~ stent + height + female +
#'       diabetic + acutemi + ejecfrac + ves1proc,
#'       data = lindner, family = binomial)
#' ps <- lindner.ps$fitted
#' cstrata.psa(abcix, ps, strata = 5)
#' cstrata.psa(abcix, ps, strata = 10)
#' cstrata.psa(abcix, ps, int = c(.37, .56, .87, 1))
#'
#' @export cstrata.psa
cstrata.psa <- function(treatment,
						propensity,
						strata = NULL,
						int = NULL,
						tree = FALSE,
						minsize = 2,
						graphic = TRUE,
						colors = c("dark blue", "dark green"),
						xlab = "Estimated Propensity Scores with Random Heights",
						pch = c(16, 16))

{
	#treatment == binary vector of 0s and 1s (necessarily? what if character, or 1, 2?)
	#propensity == PS scores from some method or other.
	#strata == either a vector of strata number for each row of covariate, or one number n in which case it is attempted to group rows by ps scores into n strata of size approximately 1/n.
	#This does not seem to work well in the case of few specific propensity values, as from a tree.
	#int == a vector of cut points between 0 and 1 defining the subintervals (perhaps as suggested by loess.psa).  In either case these subintervals define strata, so strata can be of any size.
	#tree == logical, if unique ps scores are few, as from a recursively partitioned tree, then TRUE will force each ps value to define a stratum.
	#minsize == smallest allowable stratum-treatment size.  If violated, strata is removed.  Is this the best idea?
	#colors = 2-ary color vector for the points
	#xlab == x axis label
	#ylab == y axis label
	#pch == size of points

	######################################################## BEGIN A

	if (tree &
		(!is.null(strata) |
		 !is.null(int))) {
		stop("If 'tree == TRUE' then 'strata' and 'int' must be NULL")
	}
	if (tree) {
		int <- sort(unique(propensity))
	}
	if (is.null(strata) &
		is.null(int)) {
		stop("One of 'strata' or 'int' must be defined.")
	}
	if (!is.null(strata) &
		!is.null(int)) {
		stop("Only one of 'strata' or 'int' may be defined.")
	}
	if (!is.null(int))
		int <- sort(unique(int))

	# If strata is one number n, then rows are grouped into n similarly sized (length(treatment)/n) strata
	if (length(strata) == 1)
	{
		strata <-
			findInterval(propensity,
						 vec = quantile(propensity, seq(0, 1, 1 / strata)),
						 rightmost.closed = TRUE)
	}

	#J: Defines sho (counts table) if strata are given explicitly
	if (!is.null(strata))
	{
		sho <- table(treatment, strata)
		psct <- strata
		if (length(unique(strata)) > min(26, nrow(treatment) / 5))
		{
			stop("The number of strata (before minsize adjustment) cannot exceed 26 or N/5")
		}
	}

	#J: If strata is null, rely upon int
	if (is.null(strata)) {
		if (is.null(int))
		{
			int = sort(unique(propensity))
		}
		ed <- length(int)
		if (ed > min(26, nrow(treatment) / 5))
		{
			stop("The number of strata (before minsize adjustment) cannot exceed 26 or N/5")
		}
		#   if(ed == 1)
		#     {strat <- findInterval(propensity, vec = seq(0,1,1/int), rightmost.closed = TRUE)
		#      sho = table(treatment,strat)
		#     }
		if (ed == 1) {
			stop("'int' should be a vector of propensity scores defining the strata.")
		}
		if (ed > 1)
		{
			if (int[1] > 0 & int[ed] == 1) {
				int <- c(0, int)
			}
			if (int[1] == 0 & int[ed] < 1) {
				int <- c(int, 1)
			}
			if (int[1] > 0 & int[ed] < 1) {
				int <- c(0, int, 1)
			}
			strat <-
				findInterval(propensity,
							 vec = int,
							 rightmost.closed = FALSE)
			sho = table(treatment, strat)
		}
		psct <- strat
	}

	#This is to find the cutting points of propensity score between different strata
	#It is named as cup finally
	sn = sort(unique(psct))
	alv = NULL
	for (i in 1:length(sn))
	{
		alv[i * 2 - 1] = min(propensity[psct == sn[i]])
		alv[i * 2] = max(propensity[psct == sn[i]])
	}
	alv <- sort(alv)
	cup = NULL
	for (j in 1:(length(sn) - 1))
	{
		cup[j] = (alv[j * 2] + alv[j * 2 + 1]) / 2
	}

	######################################################## END A

	######################################################## BEGIN B

	shom = matrix(sho, nrow = 2)
	dimnames(shom) = list(dimnames(sho)$treatment, dimnames(sho)$strat)

	#J: Remove strata with number of observations < 'minsize'.
	cutstr <- NULL
	for (i in 1:dim(shom)[2]) {
		if (min(shom[, i]) < minsize) {
			cutstr <- c(cutstr, i)
		}
	}
	if (length(cutstr) > 0) {
		shomn <- shom[, -cutstr]
	} else{
		shomn <- shom
	}

	#Define 'som' as matrix of strata counts for treatment/control and strata sums after minsize restriction is enforced.
	som = rbind(shomn, colSums(shomn))
	rownames(som)[3] <- "Strata.Size"
	colnames(som) <- colnames(shomn)

	######################################################## END B

	######################################################## BEGIN C

	if (graphic)
	{
		#the following codes are adapted from loess.psa,
		#minor changes are made to show sample size difference between treatment and control
		#the y axis value in this graph is randomly generated.
		height = runif(length(treatment), 1, 100)
		htp <- data.frame(height, treatment, propensity)
		dimnames(htp)[2] <-
			list(c("height", "treatment", "propensity"))
		sut <- sort(unique(treatment))
		if (!length(strata) == 0) {
			xlab <- "Strata as Input (random heights)"
		}
		if (length(strata) == 0) {
			if (ed == 1) {
				xlab <-
					paste(int,
						  "Equally Spaced Strata (random heights)",
						  sep = " ")
			} else{
				xlab <- "Strata Cutpoints as Input (random heights)"
			}
		}
		par(mfrow = c(2, 1))
		plot(
			htp$p,
			htp$h,
			type = "n",
			xlim = range(htp$p),
			xlab = "",
			ylab = "",
			main = paste("Treatment =", sut[1]),
			yaxt = "n"
		)
		points(htp$p[treatment == sut[1]], htp$h[treatment == sut[1]],
			   col = colors[1], pch = pch[1])
		abline(v = cup)

		plot(
			htp$p,
			htp$h,
			type = "n",
			xlim = range(htp$p),
			xlab = xlab,
			ylab = "",
			main = paste("Treatment =", sut[2]),
			yaxt = "n"
		)

		points(htp$p[treatment == sut[2]], htp$h[treatment == sut[2]],
			   col = colors[2], pch = pch[2])
		abline(v = cup)
		par(mfrow = c(1, 1))
	}
	######################################################## END C
	out <- list(shom, som, psct)
	names(out) <- c("Original.Strata", "Used.Strata", "strata")
	return(out)
}
