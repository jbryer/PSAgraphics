#' Fisher's Exact Test for Independence of Treatments within Strata
#'
#' Simple function that calls fisher.test repeatedly for each strata, testing
#' the independence of treatements for the given covariate within strata.
#'
#' This function makes repeated calls to fisher.test, Fisher's Exact test, to
#' test whether the distribution of the covariate categorical is independent of
#' treatment within each stratum; a list of p-values for the test for each
#' stratum are returned.
#'
#' @param categorical Categorical covariate that is being balanced within
#' strata in a PSA. If \code{categorical} has three columns, then the second
#' and third are assumed to be the treatment and strata respectively.  Missing
#' values are not allowed. May be factor or numeric.
#' @param treatment Binary variable of same length as \code{categorical};
#' generally 0 for 'control,' 1 for 'treatment.'
#' @param strata Integer variable; a vector of same length as
#' \code{categorical} indicating the derived strata from estimated propensity
#' scores.
#' @param FB Numeric; number of replications sent to fisher.test.
#' @return Returns list of the same lenght as the number of strata containing
#' p-values for the indpendence of treatment within each stratum derived from
#' Fisher's Exact test.
#' @author James E. Helmreich \email{ James.Helmreich@@Marist.edu}
#'
#' Robert M. Pruzek \email{RMPruzek@@yahoo.com}
#' @seealso \code{bal.cs.psa}, \code{bal.ms.psa}, \code{bal.ks.psa}
#' @keywords htest
#' @examples
#' #Everything random
#' categorical<-sample(4, 1000, replace = TRUE)
#' treatment<-sample(c(0,1), 1000, replace = TRUE)
#' strata<-sample(5, 1000, replace = TRUE)
#' bal.fe.psa(categorical, treatment, strata)
#'
#' #Perfect balance on 80%, random on last 20%
#' categorical<-rep(sample(5,1000, replace=TRUE), 2)
#' treatment<-c(rep(0,1000), rep(1,1000))
#' strata<-sample(6, 1200, replace=TRUE)
#' strata<-c(strata[1:1000], strata[1:800], strata[1001:1200])
#' bal.fe.psa(categorical, treatment, strata)
#'
#' @export bal.fe.psa
bal.fe.psa <-
	function(categorical,
			 treatment = NULL,
			 strata = NULL,
			 FB = 2000) {
		#This function makes repeated calls to fisher.test, Fisher's Exact test,
		#a test of whether the distribution of categorical is
		#independent of treatment within each stratum.
		#p-values for the test for each stratum are returned.

		#If "categorical" has three columns, treat as c, t, s.
		if (dim(as.data.frame(categorical))[2] == 3) {
			treatment   <- categorical[, 2]
			strata      <-
				categorical[, 3]
			categorical <-
				categorical[, 1]
		}

		if (is.factor(treatment)) {
			treatment <- factor(treatment,
								levels = levels(treatment),
								labels = c(0, 1))
		}
		categorical <- as.factor(categorical)

		n <- length(categorical)
		nstrat <- dim(table(strata))



		contingency <- table(categorical, treatment, strata)
		fe <- NULL
		for (i in 1:nstrat) {
			fe <- c(fe, fisher.test(contingency[, , i], B = FB)$p)
		}

		return(fe)
	}
