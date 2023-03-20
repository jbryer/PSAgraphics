#' Data on 996 initial Percutaneous Coronary Interventions (PCIs) performed in
#' 1997 at the Lindner Center, Christ Hospital, Cincinnati.
#'
#' Data from an observational study of 996 patients receiving a PCI at Ohio
#' Heart Health in 1997 and followed for at least 6 months by the staff of the
#' Lindner Center. This is a landmark dataset in the literature on propensity
#' score adjustment for treatment selection bias due to practice of evidence
#' based medicine; patients receiving abciximab tended to be more severely
#' diseased than those who did not receive a IIb/IIIa cascade blocker.
#'
#'
#' @name lindner
#' @docType data
#' @format A data frame with 996 observations on the following 10 variables, no NAs.
#' \describe{
#' \item{lifepres}{Mean life years preserved due to survival for at least 6 months following PCI; numeric value of either 11.4 or 0.}
#' \item{cardbill}{Cardiac related costs incurred within 6
#' months of patient's initial PCI; numeric value in 1998 dollars; costs were
#' truncated by death for the 26 patients with lifepres == 0.}
#' \item{abcix}{Numeric treatment selection indicator; 0 implies usual
#' PCI care alone; 1 implies usual PCI care deliberately augmented by either
#' planned or rescue treatment with abciximab.}
#' \item{stent}{Coronary
#' stent deployment; numeric, with 1 meaning YES and 0 meaning NO.}
#' \item{height}{Height in centimeters; numeric integer from 108 to
#' 196.}
#' \item{female}{Female gender; numeric, with 1 meaning YES and 0
#' meaning NO.}
#'  \item{diabetic}{Diabetes mellitus diagnosis; numeric,
#' with 1 meaning YES and 0 meaning NO.}
#' \item{acutemi}{Acute
#' myocardial infarction within the previous 7 days; numeric, with 1 meaning
#' YES and 0 meaning NO.}
#' \item{ejecfrac}{Left ejection fraction;
#' numeric value from 0 percent to 90 percent.}
#' \item{ves1proc}{Number
#' of vessels involved in the patient's initial PCI procedure; numeric integer
#' from 0 to 5.}
#' }
#'
#' @source Package USPS, by R. L. Obenchain.
#' @keywords datasets
NULL





#' Graphical Analysis of Variance
#'
#' A collection of functions that primarily produce graphics to aid in a
#' Propensity Score Analysis (PSA).  Functions include: cat.psa and box.psa to
#' test balance within strata of categorical and quantitative covariates,
#' circ.psa for a representation of the estimated effect size by stratum,
#' loess.psa that provides a graphic and loess based effect size estimate, and
#' various balance functions that provide measures of the balance achieved via
#' a PSA in a categorical covariate.
#'
#' \tabular{ll}{ Package: \tab PSAgraphics\cr Version: \tab 2.0\cr License:
#' \tab GPL (>= 2)\cr }
#'
#' @name granova-package
#' @aliases PSAgraphics-package PSAgraphics
#' @docType package
#' @author James E. Helmreich <James.Helmreich@@Marist.edu> and Robert M. Pruzek <RMPruzek@@yahoo.com>
#'
#' Maintainer: Jason Bryer <jason@bryer.org>
#' @seealso
#'
#' \code{\link{box.psa}} \code{\link{cat.psa}} \code{\link{circ.psa}}
#' \code{\link{loess.psa}} \code{\link{bal.ks.psa}} \code{\link{bal.ms.psa}}
#'
#' \code{\link{bal.fe.psa}} \code{\link{cstrata.psa}}
#' \code{\link{cv.trans.psa}} \code{\link{cv.bal.psa}}
#' @importFrom grDevices terrain.colors
#' @importFrom graphics abline axis box boxplot hist legend lines par points rect rug segments symbols text title
#' @importFrom stats aggregate fisher.test ks.test loess na.omit qt quantile runif sd var
NULL


#' Data on 996 initial Percutaneous Coronary Interventions (PCIs) performed in
#' 1997 at the Lindner Center, Christ Hospital, Cincinnati.
#'
#' Data from an observational study of 996 patients receiving a PCI at Ohio
#' Heart Health in 1997 and followed for at least 6 months by the staff of the
#' Lindner Center. This is a landmark dataset in the literature on propensity
#' score adjustment for treatment selection bias due to practice of evidence
#' based medicine; patients receiving abciximab tended to be more severely
#' diseased than those who did not receive a IIb/IIIa cascade blocker.
#'
#'
#' @name lindner
#' @docType data
#' @format A data frame with 996 observations on the following 10 variables, no
#' NAs. \describe{ \item{list("lifepres")}{Mean life years preserved due to
#' survival for at least 6 months following PCI; numeric value of either 11.4
#' or 0.} \item{list("cardbill")}{Cardiac related costs incurred within 6
#' months of patient's initial PCI; numeric value in 1998 dollars; costs were
#' truncated by death for the 26 patients with lifepres == 0.}
#' \item{list("abcix")}{Numeric treatment selection indicator; 0 implies usual
#' PCI care alone; 1 implies usual PCI care deliberately augmented by either
#' planned or rescue treatment with abciximab.} \item{list("stent")}{Coronary
#' stent deployment; numeric, with 1 meaning YES and 0 meaning NO.}
#' \item{list("height")}{Height in centimeters; numeric integer from 108 to
#' 196.} \item{list("female")}{Female gender; numeric, with 1 meaning YES and 0
#' meaning NO.} \item{list("diabetic")}{Diabetes mellitus diagnosis; numeric,
#' with 1 meaning YES and 0 meaning NO.} \item{list("acutemi")}{Acute
#' myocardial infarction within the previous 7 days; numeric, with 1 meaning
#' YES and 0 meaning NO.} \item{list("ejecfrac")}{Left ejection fraction;
#' numeric value from 0 percent to 90 percent.} \item{list("ves1proc")}{Number
#' of vessels involved in the patient's initial PCI procedure; numeric integer
#' from 0 to 5.} }
#' @source Package USPS, by R. L. Obenchain.
#' @keywords datasets
NULL


#' Graphical Analysis of Variance
#'
#' A collection of functions that primarily produce graphics to aid in a
#' Propensity Score Analysis (PSA).  Functions include: cat.psa and box.psa to
#' test balance within strata of categorical and quantitative covariates,
#' circ.psa for a representation of the estimated effect size by stratum,
#' loess.psa that provides a graphic and loess based effect size estimate, and
#' various balance functions that provide measures of the balance achieved via
#' a PSA in a categorical covariate.
#'
#' \tabular{ll}{ Package: \tab PSAgraphics\cr Version: \tab 2.0\cr License:
#' \tab GPL (>= 2)\cr }
#'
#' @name granova-package
#' @aliases PSAgraphics-package PSAgraphics
#' @docType package
#' @author James E. Helmreich <James.Helmreich@@Marist.edu> and
#'
#' Robert M. Pruzek <RMPruzek@@yahoo.com>
#'
#' Maintainer: James E. Helmreich <James.Helmreich@@Marist.edu>
#' @seealso
#'
#' \code{\link{box.psa}} \code{\link{cat.psa}} \code{\link{circ.psa}}
#' \code{\link{loess.psa}} \code{\link{bal.ks.psa}} \code{\link{bal.ms.psa}}
#'
#' \code{\link{bal.fe.psa}} \code{\link{cstrata.psa}}
#' \code{\link{cv.trans.psa}} \code{\link{cv.bal.psa}}
#' @keywords hplot
NULL

