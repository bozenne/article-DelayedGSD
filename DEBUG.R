### DEBUG.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun  5 2024 (13:50) 
## Version: 
## Last-Updated: sep 26 2024 (09:38) 
##           By: Brice Ozenne
##     Update #: 19
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)
library(DelayedGSD)
library(ggplot2)

MyMissProb <- matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538),  nrow = 2, ncol = 2,
                     dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")))

df.method <- data.frame(method = 1:3, binding = TRUE, fixC = c(FALSE,FALSE,TRUE))

 
args.GenData <- list(rand.block = c(1, 1, 0, 0),
                     allsd = c(2.5, 2.1, 2.4),
                     mean0 = c(10, 0, 0),
                     delta = c(0, 0.5, 1)*0,
                     ar = 30,
                     cor.01.1 = -0.15,
                     cor.ij.1 = 0.68,
                     cor.0j.1 = -0.27,
                     MissProb = MyMissProb,
                     DigitsOutcome = 2,
                     TimeFactor = 42,
                     DigitsTime = 0)
##### 2 stages #### 
debug2 <- operatingDelayedGSD(n.sim = 10, 
                              method = df.method,
                              args.GenData = args.GenData,
                              kMax = 2, InfoR.i = c(0.56,1), InfoR.d = c(0.65, 1), delta = 1,
                              PropForInterim = 0.5, lag = 21,
                              seed = 1:10)

##### 3 stages ####
res3stage[scenario == 4 & method == 3 & seed == 407971177]

debug3 <- operatingDelayedGSD(n.sim = 1, 
                              method = df.method[1:3,,drop=FALSE],
                              args.GenData = args.GenData,
                              kMax = 3, InfoR.i = c(0.40,0.65,1), InfoR.d = c(0.50,0.75, 1), delta = 1,
                              PropForInterim = c(0.35,0.6), lag = 21,
                              seed = 407971177)

debug3$results[debug3$results$method == 1,]
debug3$results[debug3$results$method == 3,]

res3 <- FinalPvalue2(Info.d = c(9.51425417, 12.00746544),  
                     Info.i = c(5.35450661, 9.31932803),  
                     ck = c(1.95996398, 1.97239755),
                     ck.unrestricted = c(1.74527964, 1.97239755),   
                     lk = c(0.18582186, -Inf),  
                     uk = c(2.38602074, Inf),
                     reason.interim = c("no boundary crossed", "Imax reached"),
                     kMax = 3, 
                     delta = 0,  
                     statistic = 1.974,
                     method = 3,  
                     bindingFutility = TRUE,
                     cNotBelowFixedc = TRUE,
                     continuity.correction = 1)

pmvnorm2(lower = c(2.38602074, 1.74527964),  
         upper = c(Inf, Inf),
         mean = c(0,0),
         sigma = cbind(c(1, 0.75019187), 
                       c(0.75019187, 1)))

pmvnorm2(lower = c(2.38602074, 1.96),  
         upper = c(Inf, Inf),
         mean = c(0,0),
         sigma = cbind(c(1, 0.75019187), 
                       c(0.75019187, 1)))

pmvnorm2(lower = c(0.18582186, 1.974), 
         upper = c(2.386, Inf),
         mean = c(0,0),
         sigma = cbind(c(1, 0.667781), 
                       c(0.667781, 1))
         )

##----------------------------------------------------------------------
### DEBUG.R ends here
