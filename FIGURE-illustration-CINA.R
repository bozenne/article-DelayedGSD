### FIGURE-illustration-CINA.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 25 2024 (16:08) 
## Version: 
## Last-Updated: sep 26 2024 (10:32) 
##           By: Brice Ozenne
##     Update #: 5
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


## * Cannot compute CI with non-binding futility rule 

## ** find problematic case
res2stage <- readRDS(file.path("Results-built","res2stage.rds"))
res2stage[binding==FALSE & is.na(lower_MUE) & !is.na(p.value_MUE),.(minp = min(p.value_MUE),medianp=median(p.value_MUE),maxp = max(p.value_MUE)),by="method"]
##    method      minp medianp  maxp
##     <int>     <num>   <num> <num>
## 1:      1 1.0000000       1     1
## 2:      2 1.0000000       1     1
## 3:      3 0.9782508       1     1
res2stage[hypo == "typeI" & method == 1 & ar == 2 & binding==FALSE & fixC==FALSE & is.na(lower_MUE) & !is.na(p.value_MUE),.(seed,p.value_MUE)]
## 551131307
res2stage[hypo == "typeI" & method == 3 & ar == 2 & binding==FALSE & fixC==FALSE & is.na(lower_MUE) & !is.na(p.value_MUE),.(seed,p.value_MUE)]
## 962697766

## ** warper to estimate rejection rate under various estimate
rejRate <- function(object, delta = NULL){
    object.GSD <- tail(object$delayedGSD,1)[[1]]
    if(is.null(delta)){delta <- seq(-1,1,by=0.001)}

    out <- do.call(rbind,lapply(delta, function(iDelta){
        iRes <- FinalPvalue2(Info.d = object.GSD$Info.d[1:object.GSD$stage$k],  
                             Info.i = object.GSD$Info.i[1:object.GSD$stage$k],  
                             ck = object.GSD$ck[1:object.GSD$stage$k],
                             ck.unrestricted = object.GSD$ck.unrestricted[1:object.GSD$stage$k], 
                             lk = object.GSD$lk[1:object.GSD$stage$k],  
                             uk = object.GSD$uk[1:object.GSD$stage$k],
                             reason.interim = object.GSD$conclusion["reason.interim",1:object.GSD$stage$k],
                             kMax = object.GSD$kMax, 
                             delta = iDelta,  
                             statistic = tail(object.GSD$delta[object.GSD$delta$method=="ML","statistic"],1),
                             method = object.GSD$method,  
                             bindingFutility = object.GSD$bindingFutility,
                             cNotBelowFixedc = object.GSD$cNotBelowFixedc,
                             continuity.correction = 1)
        return(c(delta = iDelta, p.value = iRes, attr(iRes,"terms")[1,]))
    }))
    return(as.data.frame(out))
}


## ** normal case
resSeed.normal <- operatingDelayedGSD(n.sim = 1, 
                                      method = data.frame(method = 1, binding = TRUE, fixC = FALSE),
                                      args.GenData = list(rand.block = c(1, 1, 0, 0),
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
                                                          DigitsTime = 0),
                                      kMax = 2, InfoR.i = c(0.56,1), InfoR.d = c(0.65, 1), delta = 1,
                                      PropForInterim = 0.5, lag = 21,
                                      seed = 551131307)
resSeed.normal$results[,c("stage","type.stage","info","infoPC","uk","lk","statistic","p.value_MUE","lower_MUE","upper_MUE")]
##   stage type.stage     info    infoPC       uk        lk   statistic p.value_MUE  lower_MUE upper_MUE
## 1     1    interim 5.984031 0.5419574 2.440034 0.5561203 -0.02018354          NA         NA        NA
## 2     1   decision 8.254085 0.7475499 1.672540 1.6725401  1.27212681   0.2952138 -0.5748515 0.9710938

resRate.normal <- rejRate(resSeed.normal)
range(resRate.normal$p.value)
## [1] 0.001340989 0.979118818
resRate.normal[c(which.min(abs(resRate.normal$p.value-0.025)),which.min(abs(resRate.normal$p.value-0.975))),]
##       delta   p.value     efficacy     reversal   continue
## 426  -0.575 0.0249787 5.021477e-05 0.0001479188 0.02478056
## 1972  0.971 0.9749408 4.740971e-01 0.0094746539 0.49136907

## ** problematic case (p.value=1)
resSeed.pb1 <- operatingDelayedGSD(n.sim = 1, 
                               method = data.frame(method = 1, binding = FALSE, fixC = FALSE),
                               args.GenData = list(rand.block = c(1, 1, 0, 0),
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
                                                   DigitsTime = 0),
                               kMax = 2, InfoR.i = c(0.56,1), InfoR.d = c(0.65, 1), delta = 1,
                               PropForInterim = 0.5, lag = 21,
                               seed = 551131307)
resSeed.pb1$results[,c("stage","type.stage","info","infoPC","uk","lk","statistic","p.value_MUE","lower_MUE","upper_MUE")]
##   stage type.stage      info    infoPC       uk        lk statistic p.value_MUE lower_MUE upper_MUE
## 1     1    interim  7.180679 0.6418053 2.315315 0.9426542  2.432482          NA        NA        NA
## 2     1   decision 10.059519 0.8991146 1.793917 1.7939172  1.499768           1        NA        NA

resRate.pb1 <- rejRate(resSeed.pb1)
range(resRate.pb1$p.value)
## [1] 1 1

## ** problematic case (p.value<1)
resSeed.pb2 <- operatingDelayedGSD(n.sim = 1, 
                               method = data.frame(method = 3, binding = FALSE, fixC = TRUE),
                               args.GenData = list(rand.block = c(1, 1, 0, 0),
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
                                                   DigitsTime = 0),
                               kMax = 2, InfoR.i = c(0.56,1), InfoR.d = c(0.65, 1), delta = 1,
                               PropForInterim = 0.5, lag = 21,
                               seed = 962697766)
resSeed.pb2$results[,c("stage","type.stage","info","infoPC","uk","lk","statistic","p.value_MUE","lower_MUE","upper_MUE")]
##   stage type.stage     info    infoPC       uk        lk statistic p.value_MUE lower_MUE upper_MUE
## 1     1    interim 6.405767 0.5719671 2.263137 0.5795951  2.666548          NA        NA        NA
## 2     1   decision 7.805370 0.6969368 2.096479 2.0964792  2.050658   0.9967586        NA 0.6585938
resRate.pb2 <- rejRate(resSeed.pb2)
range(resRate.pb2$p.value)
## [1] 0.9835125 0.9999994

## ** figure
resRateL.normal <- reshape(cbind(resRate.normal, sum = resRate.normal$efficacy + resRate.normal$reversal + resRate.normal$continue), direction = "long",
                        varying = c("efficacy","reversal","continue","sum"), idvar = "delta",
                        timevar = "term", v.names = "value", times = c("efficacy","reversal","continue","sum"))
resRateL.pb1 <- reshape(cbind(resRate.pb1, sum = resRate.pb1$efficacy + resRate.pb1$reversal + resRate.pb1$continue), direction = "long",
                        varying = c("efficacy","reversal","continue","sum"), idvar = "delta",
                        timevar = "term", v.names = "value", times = c("efficacy","reversal","continue","sum"))
resRateL.pb2 <- reshape(cbind(resRate.pb2, sum = resRate.pb2$efficacy + resRate.pb2$reversal + resRate.pb2$continue), direction = "long",
                        varying = c("efficacy","reversal","continue","sum"), idvar = "delta",
                        timevar = "term", v.names = "value", times = c("efficacy","reversal","continue","sum"))

resRateL.all <- rbind(cbind(resRateL.normal, type = paste0("CI=[",round(resSeed.normal$results$lower_MUE[2],3),
                                                           ";",round(resSeed.normal$results$upper_MUE[2],3),"] \t p=",
                                                           round(resSeed.normal$results$p.value_MUE[2],3))),
                      cbind(resRateL.pb1, type = paste0("CI=[",round(resSeed.pb1$results$lower_MUE[2],3),
                                                        ";",round(resSeed.pb1$results$upper_MUE[2],3),"] \t p=",
                                                        round(resSeed.pb1$results$p.value_MUE[2],3))),
                      cbind(resRateL.pb2, type = paste0("CI=[",round(resSeed.pb2$results$lower_MUE[2],3),
                                                           ";",round(resSeed.pb2$results$upper_MUE[2],3),"] \t p=",
                                                           round(resSeed.pb2$results$p.value_MUE[2],3))))
resRateL.all$type <- factor(resRateL.all$type, levels = unique(resRateL.all$type))

figure.CINA <- ggplot(resRateL.all, aes(x=delta,y=value, group=term,color=term))
figure.CINA <- figure.CINA + geom_hline(yintercept = c(0.025,0.975), color = "black", linetype = 2)
figure.CINA <- figure.CINA + facet_wrap(~type) + labs(x = "Effect size", y = "Probability of a more extreme results") 
figure.CINA <- figure.CINA + geom_line(aes(linewidth = term))
figure.CINA <- figure.CINA + scale_linewidth_manual(values = c("continue" = 1, 
                                                               "efficacy" = 1,
                                                               "reversal" = 1,
                                                               "sum" = 2))
figure.CINA

ggsave(filename = "figures/review-CINA.pdf", figure.CINA, width = 10, height = 5)


##----------------------------------------------------------------------
### FIGURE-illustration-CINA.R ends here
