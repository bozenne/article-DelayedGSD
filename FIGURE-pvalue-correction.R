### FIGURE-pvalue-correction.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 10 2023 (13:16) 
## Version: 
## Last-Updated: jun  5 2024 (14:49) 
##           By: Brice Ozenne
##     Update #: 50
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(DelayedGSD)
library(RColorBrewer)
DelayedGSD.options(FCT.p_value = "FinalPvalue2")
export <- TRUE

## * 2 stages

## ** bounds
bound2stage.m12b <- CalcBoundaries(kMax = 2,alpha = 0.025, beta = 0.2, InfoR.i = c(0.58, 1.00), InfoR.d = c(0.68, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                             method = 1, cNotBelowFixedc = FALSE, bindingFutility = TRUE)
bound2stage.m12b.fixC <- CalcBoundaries(kMax = 2, alpha = 0.025, beta = 0.2, InfoR.i = c(0.58, 1.00), InfoR.d = c(0.68, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                                  method = 1, cNotBelowFixedc = TRUE, bindingFutility = TRUE)
bound2stage.m3b <- CalcBoundaries(kMax = 2, alpha = 0.025, beta = 0.2, InfoR.i = c(0.58, 1.00), InfoR.d = c(0.68, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                            method = 3, cNotBelowFixedc = TRUE, bindingFutility = TRUE)
bound2stage.m12nb <- CalcBoundaries(kMax = 2, alpha = 0.025, beta = 0.2, InfoR.i = c(0.58, 1.00), InfoR.d = c(0.68, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                              method = 1, cNotBelowFixedc = FALSE, bindingFutility = FALSE)
bound2stage.m12nb.fixC <- CalcBoundaries(kMax = 2, alpha = 0.025, beta = 0.2, InfoR.i = c(0.58, 1.00), InfoR.d = c(0.68, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                                   method = 1, cNotBelowFixedc = TRUE, bindingFutility = FALSE)
bound2stage.m3nb <- CalcBoundaries(kMax = 2, alpha = 0.025, beta = 0.2, InfoR.i = c(0.58, 1.00), InfoR.d = c(0.68, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                             method = 3, cNotBelowFixedc = TRUE, bindingFutility = FALSE)

## alpha-spent
##     0.00841
##     0.02500

## ** p-values
DelayedGSD.options(FCT.p_value = "FinalPvalue2")
pval2stage.m12b <- gridFinalPvalue(bound2stage.m12b)
pval2stage.m12b.fixC <- gridFinalPvalue(bound2stage.m12b.fixC)
pval2stage.m3b <- gridFinalPvalue(bound2stage.m3b)
pval2stage.m12nb <- gridFinalPvalue(bound2stage.m12nb)
pval2stage.m12nb.fixC <- gridFinalPvalue(bound2stage.m12nb.fixC)
pval2stage.m3nb <- gridFinalPvalue(bound2stage.m3nb)

attr(pval2stage.m3nb,"terms")
## $`stage=1: z=-3`
##         efficacy  reversal  continue
## [1,] 0.009738676 0.7459376 0.9902613

## $`stage=1: z=0`
##         efficacy  reversal  continue
## [1,] 0.009738676 0.2500008 0.9902613

## $`stage=1: z=1.95995398454005`
##         efficacy     reversal  continue
## [1,] 0.008410045 7.994106e-06 0.9902613

## $`stage=2: z=-3`
##         efficacy reversal continue
## [1,] 0.008409995        0        0
## [2,] 0.988911426        0        0

## ** graphical display
if(export){
    pdf("figures/illustration-pvalue-2stage.pdf", width = 11, height = 10)
}

par(mfrow = c(2,3), mar = rep(2.5,4))
plot(pval2stage.m12b, xlim = c(0.7,2.3), title = "method 1 or 2")
plot(pval2stage.m12b.fixC, xlim = c(0.7,2.3), title = "method 1 or 2 (ck>1.96)")
plot(pval2stage.m3b, xlim = c(0.7,2.3), title = "method 3")
plot(pval2stage.m12nb, xlim = c(0.7,2.3), title = "method 1 or 2 (non-binding)")
plot(pval2stage.m12nb.fixC, xlim = c(0.7,2.3), title = "method 1 or 2 (ck>1.96, non-binding)")
plot(pval2stage.m3nb, xlim = c(0.7,2.3), title = "method 3 (non-binding)")

if(export){
    dev.off()
}

## * 3 stages

## ** bounds
bound3stage.m12b <- CalcBoundaries(kMax = 3,alpha = 0.025, beta = 0.2, InfoR.i = c(0.40,0.65,1), InfoR.d = c(0.50,0.75, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                                   method = 1, cNotBelowFixedc = FALSE, bindingFutility = TRUE)
bound3stage.m12b.fixC <- CalcBoundaries(kMax = 3,alpha = 0.025, beta = 0.2, InfoR.i = c(0.40,0.65,1), InfoR.d = c(0.50,0.75, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                                        method = 1, cNotBelowFixedc = TRUE, bindingFutility = TRUE)
bound3stage.m3b <- CalcBoundaries(kMax = 3,alpha = 0.025, beta = 0.2, InfoR.i = c(0.40,0.65,1), InfoR.d = c(0.50,0.75, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                                  method = 3, cNotBelowFixedc = TRUE, bindingFutility = TRUE)
bound3stage.m12nb <- CalcBoundaries(kMax = 3,alpha = 0.025, beta = 0.2, InfoR.i = c(0.40,0.65,1), InfoR.d = c(0.50,0.75, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                                    method = 1, cNotBelowFixedc = FALSE, bindingFutility = FALSE)
bound3stage.m12nb.fixC <- CalcBoundaries(kMax = 3,alpha = 0.025, beta = 0.2, InfoR.i = c(0.40,0.65,1), InfoR.d = c(0.50,0.75, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                                         method = 1, cNotBelowFixedc = TRUE, bindingFutility = FALSE)
bound3stage.m3nb <- CalcBoundaries(kMax = 3,alpha = 0.025, beta = 0.2, InfoR.i = c(0.40,0.65,1), InfoR.d = c(0.50,0.75, 1), rho_alpha = 2, rho_beta = 2, delta = 0.6,
                                   method = 3, cNotBelowFixedc = TRUE, bindingFutility = FALSE)

## alpha-spent
##     0.00400
##     0.01056
##     0.02500

## ** p-values
DelayedGSD.options(FCT.p_value = "FinalPvalue2")
pval3stage.m12b <- gridFinalPvalue(bound3stage.m12b)
pval3stage.m12b.fixC <- gridFinalPvalue(bound3stage.m12b.fixC)
pval3stage.m3b <- gridFinalPvalue(bound3stage.m3b)
pval3stage.m12nb <- gridFinalPvalue(bound3stage.m12nb)
pval3stage.m12nb.fixC <- gridFinalPvalue(bound3stage.m12nb.fixC)
pval3stage.m3nb <- gridFinalPvalue(bound3stage.m3nb)

attr(pval3stage.m3nb,"terms")
## $`stage=1: z=-3`
##         efficacy  reversal  continue
## [1,] 0.004418901 0.4873815 0.9955811

## $`stage=1: z=0`
##         efficacy   reversal  continue
## [1,] 0.004418901 0.06828444 0.9955811

## $`stage=1: z=1.95995398454005`
##         efficacy     reversal  continue
## [1,] 0.004000049 1.839197e-07 0.9955811

## $`stage=2: z=-3`
##         efficacy  reversal  continue
## [1,] 0.004000033 0.0000000 0.0000000
## [2,] 0.007516245 0.8080949 0.9880661

## ** graphical display
if(export){
    pdf("figures/illustration-pvalue-3stage.pdf", width = 11, height = 10)
}

par(mfrow = c(2,3), mar = rep(2.5,4))
plot(pval3stage.m12b, xlim = c(0.7,3.3), title = "method 1 or 2")
plot(pval3stage.m12b.fixC, xlim = c(0.7,3.3), title = "method 1 or 2 (ck>1.96)")
plot(pval3stage.m3b, xlim = c(0.7,3.3), title = "method 3")
plot(pval3stage.m12nb, xlim = c(0.7,3.3), title = "method 1 or 2 (non-binding)")
plot(pval3stage.m12nb.fixC, xlim = c(0.7,3.3), title = "method 1 or 2 (ck>1.96, non-binding)")
plot(pval3stage.m3nb, xlim = c(0.7,3.3), title = "method 3 (non-binding)")

if(export){
    dev.off()
}




##----------------------------------------------------------------------
### FIGURE-pvalue-correction.R ends here
