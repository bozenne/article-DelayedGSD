### FIGURE-illustration-trial.R --- g
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 28 2024 (15:41) 
## Version: 
## Last-Updated: jun 28 2024 (15:44) 
##           By: Brice Ozenne
##     Update #: 4
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(DelayedGSD)
data(simLP3trial, package = "DelayedGSD")

GSD.bound <- CalcBoundaries(kMax=2,  ## max number of analyses (including final)
                            alpha=0.025,  ## type I error
                            beta=0.1,  ## type II error
                            InfoR.i=56/100,  ## planned or observed information rates
                            rho_alpha=2,  ## rho parameter for alpha error spending function
                            rho_beta=2,  ## rho parameter for beta error spending function
                            method=3,  ## use method 1 or 2 from paper H&J
                            delta=1,  ## effect that the study is powered for
                            InfoR.d=c(65,100)/100,
                            cNotBelowFixedc=TRUE)

resI <- analyzeData(SelectData(simLP3trial,453.6779))
GSD.interim <- update(GSD.bound,delta=resI)

resD <- analyzeData(simLP3trial[1:166,])
GSD.decision <- update(GSD.interim,delta=resD)

pdf("figures/figure-article-illustration-trial.pdf", width = 7.5, height = 5)
plot(GSD.decision)
dev.off()

cairo_ps("figures/figure-article-illustration-trial.eps", width = 7.5, height = 5)
plot(GSD.decision)
dev.off()

##----------------------------------------------------------------------
### FIGURE-illustration-trial.R ends here
