#' Multiple Covariate Balance Assessment Plot
#'
#' Provides a graphic that depicts covarite effect size differences between
#' treatment groups both before and after stratification.  Function will create
#' stata internally if desired, and returns numerical output used to create
#' graphic.
#'
#' Effect sizes between treatments for each covariate are presented in one
#' graphic, both before and after stratification.
#'
#' @param covariates Dataframe of covariates.  Factors should be recoded using
#' \code{cv.trans.psa}
#' @param treatment Binary vector or factor defining the two treatments
#' @param propensity Vector of same length as \code{treatment} containing
#' estimated propensity scores.
#' @param strata Either a vector of same length as \code{treatment} of
#' predefined stratum number, or one integer \code{n} used to assign rows to
#' \code{n} strata \code{propensity} scores, each of approximately the same
#' number of cases. If relatively few unique propensity scores have been
#' defined (as from a classification tree) then the logical \code{tree} should
#' be set equal to \code{TRUE}.
#' @param int Either a number \code{m} used to divide \code{[0,1]} into
#' \code{m} equal length subintervals, or a vector containing cut points
#' between 0 and 1 that define subintervals (perhaps as suggested by
#' loess.psa). In either case the subintervals define strata, for which sizes
#' can differ.
#' @param tree Logical, default \code{FALSE}.  If there are few unique
#' propensity scores, say from a recursively partitioned tree, then TRUE forces
#' strata to be defined by the unique propensity scores.
#' @param minsize Smallest allowable stratum-treatment size.  If violated, rows
#' in the stratum are removed.  User may wish to redefine strata.
#' @param universal.psd Logical, default = TRUE.  Forces standard deviations
#' used to be unadjusted for stratification.
#' @param trM Numeric, default = 0; passed to \code{mean} for trimming
#' purposes.
#' @param absolute.es Logical, default TRUE. If TRUE, graphic depicts absolute
#' values of all effect sizes.  Note that the adjusted effect size plotted is
#' the absolute value of weighted averages of the signed by-stratum effect size
#' values when \code{absolute.es} is TRUE.
#' @param trt.value Character string; if desired allows the name of an active
#' treatment to be given.  Should be a level (value) of the \code{treatment}
#' factor (vector).
#' @param use.trt.var Logical, default FALSE.  If TRUE, uses just active
#' treatment standard deviations for effect size, as per a suggestion of Rubin
#' and Stuart (see reference below).
#' @param verbose Logical, default FALSE.  Numerical output is returned
#' invisibly.
#' @param xlim Binary vector passed to plot for overriding default choices.
#' Default NULL.
#' @param plot.strata Logical, default TRUE.  Adds effect size values for
#' individual strata to graphic.
#' @param \dots Other graphical parameters passed to \code{plot}.
#' @return Graphic plots covariate balance before and after stratication on
#' propensity scores.  The default version (absolute.es = TRUE) plots the
#' absolute values of effect sizes for each stratum, though the overall
#' estimate is the weighted mean before taking the absolute values.  Numerical
#' output consists of seven addressable objects.  If \code{verbose} is FALSE
#' (default), output is not printed. \item{original.strata}{Matrix of
#' strata-treatment counts as originally input.} \item{strata.used}{Matrix of
#' strata-treatment counts used in effectsize calculations after any
#' \code{minsize} reductions.} \item{mean.diff.strata.wtd}{Matrix of strata by
#' covariate weighted (by strata size) average differences.}
#' \item{mean.diff.unadj}{Matrix of covariate effects sizes before
#' stratification.} \item{effect.sizes}{Matrix of effect sizes by covariate and
#' statum.} \item{treatment.levels}{Names of treatments.}
#' \item{effects.strata.treatment}{Matrix of standard deviations and
#' stratum-treatment covariate means used to calculate the
#' \code{effect.sizes}.}
#' @author Robert M. Pruzek \email{RMPruzek@@yahoo.com}
#'
#' James E. Helmreich \email{James.Helmreich@@Marist.edu}
#'
#' KuangNan Xiong \email{harryxkn@@yahoo.com}
#' @seealso \code{\link{cv.bal.psa}}, \code{\link{loess.psa}},
#' \code{\link{cstrata.psa}}, \code{\link{cv.trans.psa}}
#' @references ``Matching Methods for Causal Inference: A review and a look
#' forward." Forthcoming in Statistical Science.
#' @keywords hplot
#' @examples
#'
#' data(lindner)
#' attach(lindner)
#' lindner.ps <- glm(abcix ~ stent + height + female +
#'       diabetic + acutemi + ejecfrac + ves1proc,
#'       data = lindner, family = binomial)
#' ps<-lindner.ps$fitted
#' lindner.cv <- lindner[,4:10]
#' cv.bal.psa(lindner.cv, abcix, ps, strata = 5)
#' cv.bal.psa(lindner.cv, abcix, ps, strata = 10)
#' cv.bal.psa(lindner.cv, abcix, ps, int = c(.2, .5, .6, .75, .8))
#'
#' @export cv.bal.psa
cv.bal.psa <-
   function(covariates,
            treatment,
            propensity,
            strata = NULL,
            int = NULL,
            tree = FALSE,
            minsize = 2,
            universal.psd = TRUE,
            trM = 0,
            absolute.es = TRUE,
            trt.value = NULL,
            use.trt.var = FALSE,
            verbose = FALSE,
            xlim = NULL,
            plot.strata = TRUE,
            ...) {
      #covariates == dataframe of interest
      #treatment == binary vector of 0s and 1s (necessarily? what if character, or 1, 2?)
      #propensity == PS scores from some method or other.
      #strata == either a vector of strata number for each row of covariate, or one number n in which case it is attempted to group rows by ps scores into n strata of size approximately 1/n.  This does not seem to work well in the case of few specific propensity values, as from a tree.
      #int == either a number m used to divide [0,1] into m equal length subintervals, or a vector of cut points between 0 an    1 defining the subintervals (perhaps as suggested by loess.psa).  In either case these subintervals define strata, so strata can be of any size.
      #tree == logical, if unique ps scores are few, as from a recursively partitioned tree, then TRUE will force each ps value to define a stratum.
      #minsize == smallest allowable stratum-treatment size.  If violated, strata is removed.  Is this the best idea?
      #universal.psd == If 'TRUE', forces standard deviations used to be unadjusted for stratification.
      #trM == trimming proportion for mean calculations.
      #absolute.es == logical, if 'TRUE' routine uses absolute values of all effect sizes.
      #trt.value == allows user to specify which value is active treatment, if desired.
      #use.trt.var == logical, if true then Rubin-Stuart method using only treatment variance with be used in effect size calculations.
      #verbose == logical, controls output that is visibly returned.
      #xlim == binary vector, can be used to override default horizontal limits.  Is this really necessary?
      #Note: effect sizes are calculated as treatment 1 - treatment 0, or treatment B - treatment A.


      ########################################################
      X <- covariates

      treat.lev <- sort(unique(treatment))

      if (is.null(trt.value))
         trt.value = treat.lev[2]
      if (!is.null(trt.value))
      {
         if ((trt.value != treat.lev[2]) & (trt.value != treat.lev[1]))
         {
            stop("WARNING: trt.value as defined does not match a treatment level")
         }
         if (trt.value == treat.lev[1]) {
            treat.lev <- treat.lev[c(2, 1)]
         }
      }
      ######################################################## BEGIN C
      # Call cstrat.psa for definitions of strata
      cstrat <-
         cstrata.psa(
            treatment = treatment,
            propensity = propensity,
            strata = strata,
            int = int,
            tree = tree,
            minsize = minsize,
            graphic = FALSE
         )
      shom <- cstrat$Original.Strata
      som <- cstrat$Used.Strata
      psct <- cstrat$strata

      XX <- cbind(X, treatment, psct)

      names.cov <- colnames(X)
      n.cov <- length(names.cov)
      names.strata <- colnames(som)
      n.strata <- length(names.strata)
      ######################################################## END C

      ######################################################## BEGIN D
      #Calculating balance effect sizes for stratified covraiates
      #Initializing matrices used later.

      uess = matrix(0, nrow = n.cov, ncol = 2)
      effect.size.ji = matrix(0, nrow = n.cov, ncol = n.strata)
      effect.size.ji.adj = matrix(0, nrow = n.cov, ncol = n.strata)
      var.cov.by.strattreat = matrix(0, nrow = n.cov, ncol = 2 * n.strata)
      mean.diff <- matrix(0, nrow = n.cov, ncol = n.strata)
      mean.diff.adj = matrix(0, nrow = n.cov, ncol = n.strata)
      sd.adj <- matrix(0, nrow = n.cov, ncol = n.strata)
      sd.un <- rep(0, n.cov)
      mean.diff.unadj <- rep(0, nrow = n.cov)

      #mean.diff is the mean difference between tr/cr for each covariate across strata
      #sd.adj is the PS-adjusted pooled standard deviation for each covariate across strata
      #mean.diff.adj is the mean differences adjusted by strata size
      #var.cov.by.strattreat is the variance for each cell in strata*tr/ct table
      #################################################################Z
      for (j in 1:n.cov)
      {
         for (i in 1:n.strata)
         {
            #ha  is (size ith strata by 2)-matrix that picks off values of jth covariate and treatment in ith stratum
            ha = XX[XX[, n.cov + 2] == names.strata[i], c(j, n.cov + 1)]
            mean.diff[j, i] = (mean(ha[ha[, 2] == treat.lev[2], 1], trim = trM) -
                                  mean(ha[ha[, 2] == treat.lev[1] , 1], trim = trM))
            mean.diff.adj[j, i] = mean.diff[j, i] * som[3, i] / sum(som[3,])
            if (use.trt.var)
            {
               var.cov.by.strattreat[j, i]            = var(ha[ha[, 2] == treat.lev[2], 1])
               var.cov.by.strattreat[j, i + n.strata] = var(ha[ha[, 2] == treat.lev[2], 1])
               sd.adj[j, i] = sd(ha[ha[, 2] == treat.lev[2], 1])
            } else
            {
               var.cov.by.strattreat[j, i] =              var(ha[ha[, 2] == treat.lev[1], 1])
               var.cov.by.strattreat[j, i + n.strata] =  var(ha[ha[, 2] == treat.lev[2], 1])
               sd.adj[j, i] = sqrt((var.cov.by.strattreat[j, i] + var.cov.by.strattreat[j, i + n.strata]) /
                                      2)
            }
         }
         # uess[j,1] contains unadjusted ES for jth covariate by direct calculation
         # mean.diff.unadj and sd.un are mean.diff and sd.adj for each covariate without propensity score adjustment.
         mean.diff.unadj[j] = (mean(XX[XX[, n.cov + 1] == treat.lev[2], j], trim = trM) -
                                  mean(XX[XX[, n.cov + 1] == treat.lev[1], j], trim = trM))
         if (use.trt.var) {
            sd.un[j] = sd(XX[XX[, n.cov + 1] == treat.lev[2] , j])
         } else{
            sd.un[j] = sqrt((var(XX[XX[, n.cov + 1] == treat.lev[1], j]) +
                                var(XX[XX[, n.cov + 1] == treat.lev[2], j])) /
                               2)
         }
         uess[j, 1] = if (sd.un[j] > 0) {
            mean.diff.unadj[j] / sd.un[j]
         } else{
            0
         }

         #effect.size.ji provides the effect size mean.diff/sd.adj for each stratum for each covariate;
         #these are shown as letters on graphic.
         #effect.size.ji.adj is the middle step to calculate uess[j,2] which is the sum of mean.diff/psd in all strata
         #weighted by strata sizes.
         #uess[j,2] is going to be shown as weighted-avg of above dots from effect.sizeji in the final graphic.

         if (universal.psd == TRUE) {
            sd.adj[j, ] = sd.un[j]
         }
         for (i in 1:n.strata)
         {
            effect.size.ji[j, i] = if (sd.adj[j, i] > 0) {
               mean.diff[j, i] / sd.adj[j, i]
            } else{
               0
            }
            effect.size.ji.adj[j, i] = if (sd.adj[j, i] > 0) {
               mean.diff.adj[j, i] / sd.adj[j, i]
            } else{
               0
            }
         }
         uess[j, 2] = sum(effect.size.ji.adj[j, ])
      }
      #####################################################################Z

      #Name dimensions of everything.

      n.strata2 = n.strata * 2
      sd.un <- matrix(sd.un, ncol = 1)
      rownames(sd.un) <- names.cov
      dimnames(uess) = list(names.cov, c("stES_unadj", "stES_adj"))
      dimnames(mean.diff) = list(names.cov, names.strata)
      dimnames(mean.diff.adj) = list(names.cov, names.strata)
      dimnames(effect.size.ji) = list(names.cov, names.strata)
      mean.diff.unadj <- matrix(mean.diff.unadj, ncol = 1)
      dimnames(mean.diff.unadj) = list(names.cov, "mean.diff_unadj")
      dimnames(var.cov.by.strattreat) = list(names.cov, paste("cellvar", 1:n.strata2))

      #when absolute.es is set as true, take absolute values.
      if (absolute.es == TRUE)
      {
         effect.size.ji = abs(effect.size.ji)
         mean.diff = abs(mean.diff)
         mean.diff.adj = abs(mean.diff.adj)
         mean.diff.unadj = abs(mean.diff.unadj)
         uess = abs(uess)
      }

      se = order(uess[, 1])
      se2 = order(uess[, 1], decreasing = TRUE)
      sd.un = as.matrix(sd.un[se2,])
      colnames(sd.un) <- "st.dev.unadj"
      ord.uess = uess[se, ]
      ord.uess.2 = uess[se2, ]
      #matrix is ordered according to unadjusted ES values
      effect.size.ji1 = effect.size.ji[se,]
      effect.size.ji2 = effect.size.ji[se2,]
      mean.diff.adj = mean.diff.adj[se2,]
      mean.diff.unadj = mean.diff.unadj[se2,]
      #sd.adj.3 = sd.adj.3[se2, ]
      var.cov.by.strattreat = var.cov.by.strattreat[se2,]
      mean.diff = mean.diff[se2,]
      colnames(effect.size.ji2) =  letters[1:n.strata]
      colnames(som) = letters[1:n.strata]
      colnames(mean.diff) = letters[1:n.strata]
      colnames(mean.diff.adj) = letters[1:n.strata]
      colnames(var.cov.by.strattreat) = c(paste(letters[1:n.strata], "_", treat.lev[1], sep = ""),
                                          paste(letters[1:n.strata], "_", treat.lev[2], sep = ""))
      # final matrix contains all the final ES values as well as stratum-specific values for plotting.
      final = cbind(ord.uess, effect.size.ji1)
      ######################################################## END D

      ###################### PLOTTING ##################################
      if (is.null(xlim))
      {
         xlim <- 1.1 * range(final)
         if (xlim[1] > 0)
            xlim[1] <- 0
      }

      #Define the y-coordinates for the plot
      y.coord <- matrix(rep(1:n.cov, (2 + n.strata)), nrow = n.cov)

      main <- "Standardized Covariate Effect Sizes
w/ & w/o PS adjustment"
      if (absolute.es == TRUE) {
         main <- paste("Absolute", main)
      }

      plot(
         final,
         y.coord,
         xlim = xlim,
         ylim = NULL,
         axes = FALSE,
         sub = 'Open circles are stES-unadj; Closed circles are stES-adj; Letters represent strata' ,
         xlab = paste(
            "Standardized Effect Sizes: treatment",
            treat.lev[2],
            "- treatment ",
            treat.lev[1]
         ),
         ylab = ' ',
         main = main,
         font = 2,
         cex = 0,
         ...
      )

      axis(1, font = 2, las = 1)
      axis(
         2,
         1:n.cov,
         dimnames(ord.uess)[[1]],
         font = 2,
         las = 1,
         tick = FALSE,
         cex.axis = .79
      )

      #J: Plotting the red open circle line and the solid blue circle line
      points(final[, 1], y.coord[, 1], col = "dark red", pch = 21)
      lines(final[, 1], y.coord[, 1],  col = "dark red", lwd = 1.5)
      points(final[, 2], y.coord[, 2], col = "blue", pch = 19)
      lines(final[, 2], y.coord[, 2],  col = "blue", lwd = 1.5)

      #J: loop below plots blue letters for each covariate/stratum
      #pch = 94 + m,
      if (plot.strata) {
         for (m in 3:(2 + n.strata)) {
            points(
               final[, m],
               y.coord[, m],
               col = "blue",
               pch = 94 + m,
               cex = .66
            )
         }
      }
      abline(v = 0, lty = 2, lwd = 1.2)
      abline(v = seq(round(xlim[1] - 1, 0), round(xlim[2] + 1, 0), .5),
             lty = 3,
             lwd = .7)
      box()
      ####################### END PLOTTING ################################

      ####################### OUTPUT ################################
      sd.ESs = apply(effect.size.ji2, 1, sd)
      final2 = round(cbind(ord.uess.2, effect.size.ji2, sd.ESs), 2)
      out <-
         list(
            shom,
            som,
            round(mean.diff.adj, 2),
            round(mean.diff.unadj, 2),
            final2,
            treat.lev,
            round(cbind(sd.un, var.cov.by.strattreat), 2)
         )
      names(out) <-
         c(
            "original.strata",
            "strata.used",
            "mean.diff.strata.wtd",
            "mean.diff.unadj",
            "effect.sizes",
            "treatment.levels",
            "effects.strata.treatment"
         )
      if (verbose) {
         return(out)
      } else{
         return(invisible(out))
      }
   }
