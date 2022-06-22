#' Kolgomorov-Smirnov 2 sample tests for multiple strata
#'
#' Automates the Kolgomorov-Smirnov 2-sample nonparametric test of equivalence
#' of two distrbutions across multiple pairs of sample distributions.
#'
#' Makes multiple calls to \code{ks.test}, returning a vector of p-values
#' associated with strata from a Propensity Score Analysis.
#'
#' @param continuous Quantitative covariate that is being balanced within
#' strata in a PSA. If \code{continuous} has three columns, then the second and
#' third are assumed to be the treatment and strata respectively.  Missing
#' values are not allowed.
#' @param treatment Binary variable of same length as \code{continuous};
#' generally 0 for 'control,' 1 for 'treatment.'
#' @param strata Integer variable (usually 1 - 5); A vector of same length as
#' continuous indicating the derived strata from estimated propensity scores.
#' Generally 5 or 6 strata are used, but graph works reasonably well at least
#' up to 10 strata.
#' @return Returns a vector of same length as the number of strata containing
#' the p-values from the KS-test of equivalence of distributions for each
#' stratum-treatment pair.
#' @author James E. Helmreich \email{James.Helmreich@@Marist.edu}
#'
#' Robert M. Pruzek \email{RMPruzek@@yahoo.com}
#' @seealso \code{bal.ms.psa}, \code{bal.cs.psa}, \code{bal.cws.psa}
#' @keywords htest
#' @examples
#'
#' continuous<-rnorm(1000)
#' treatment<-sample(c(0,1),1000,replace=TRUE)
#' strata<-sample(5,1000,replace=TRUE)
#' bal.ks.psa(continuous,treatment,strata)
#'
#' @export bal.ks.psa
bal.ks.psa <- function(continuous,
					   treatment = NULL,
					   strata = NULL) {
	#Computes two sample Kolmogorov-Smirnov statistics and p-values for
	#the comparison of the two distrubition of a continuous covariate
	#within each defined stratum in a PSA.

	#If "continuous" has three columns, treat as m, t, s.
	if (dim(as.data.frame(continuous))[2] == 3) {
		treatment   <- continuous[, 2]
		strata      <-
			continuous[, 3]
		continuous <-
			continuous[, 1]
	}
	nstrat <- dim(table(strata))
	strat.f <- as.factor(strata)
	levels(strat.f) <- 1:nstrat

	ks <- NULL
	for (j in 1:nstrat) {
		kol.sm <-
			ks.test(continuous[treatment == unique(treatment)[1] &
							   	strat.f == j], continuous[treatment == unique(treatment)[2] &
							   							  	strat.f == j])
		ks[j] <- kol.sm$p
	}
	return(ks)
}
