### TEXT-simulation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun  5 2024 (13:55) 
## Version: 
## Last-Updated: jan 21 2025 (10:20) 
##           By: Brice Ozenne
##     Update #: 142
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

if(Sys.info()["login"] == "bozenne"){
}else if(Sys.info()["login"] == "hpl802"){
  setwd("x:/DelayedGSD/")
}

options(width = 120, digits = 4)

library(data.table)
library(ggplot2)
library(xtable)

## * formating functions
##' @description Average and nicely format as percentage (fixed number of digits, e.g. 0.010%)
mean2pc <- function(x, digits = 2){
    if(length(x)==0){return(NA)}
    mean.x <- mean(x, na.rm = TRUE)
    out <- paste0(formatC(100*mean.x, format='f', digits = digits), "%")
    if(round(100*mean.x,digits)==0){
        out[round(100*mean.x,digits)==0] <- paste0("<0.",rep(0,digits-1),"1%")
    }
    if(abs(mean.x)<1e-12){
        out[abs(mean.x)<1e-12] <- "0"
    }
    if(any(is.na(x))){
        out <- paste0(out," (NA: ",mean2pc(is.na(x), digits = digits),")")
    }
    return(out)
}

##' @description Average and nicely format as numeric (fixed number of digits, e.g. 0.010%)
mean2num <- function(x, digits = 3){
    if(length(x)==0){return(NA)}
    out <- formatC(mean(x, na.rm = TRUE), format='f', digits = digits)
    if(any(is.na(x))){
        out <- paste0(out," (NA: ",mean2pc(mean(is.na(x)), digits = digits),")")
    }
    return(out)
}

## * load data
res2stage <- readRDS(file.path("Results-built","res2stage.rds"))
res2stage[, method.char := paste0("method ",method, c(""," fixC")[fixC+1])]
res2stage[, stage.char := factor(stage, 1:2, c("interim","final"))]
res2stage[, truth := ifelse(hypo=="power",1,0)]

res3stage <- readRDS(file.path("Results-built","res3stage.rds"))
res3stage[, method.char := paste0("method ",method, c(""," fixC")[fixC+1])]
res3stage[, stage.char := factor(stage, 1:3, c("interim1","interim2","final"))]
res3stage[, truth := ifelse(hypo=="power",1,0)]

## * Main text
## ** sample size
res2stage[type == "final", .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                             n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                             n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]")), by = c("stage","type")]
##    stage     type   n.patients    n.outcome  n.missing
## 2:     2    final 268[237;268] 240[237;240]   28[0;28]

res3stage[type=="final", .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                           n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                           n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]")), by = c("stage","type")]
## +    stage   type   n.patients    n.outcome n.missing
## 1:     3  final 274[242;274] 245[241;245]  29[0;29]


## ** Type 1 error & power

## monte carlo error
quantile(100*sapply(1:20000, function(i){mean(rbinom(2e4, size = 1, prob = 0.025))}), probs = c(0.025,0.975))
##  2.5% 97.5% 
## 2.285 2.715 
## quantile(100*sapply(1:10000, function(i){mean(rbinom(1e5, size = 1, prob = 0.025))}), probs = c(0.025,0.975))
##  2.5% 97.5% 
## 2.402 2.597 

## *** 2 stages
## For each run, create a binary indicator for rejection for efficacy
res2stage.rejection <- res2stage[type != "interim" & infoBias == 0,.(n.stage = .N, rejection = "efficacy" %in% decision),
                                 by = c("method.char","seed","scenario","missing","binding","fixC","ar","hypo")]

## Average over runs and method within scenario
res2stageS.rejection <- res2stage.rejection[,.(n.sim = .N, rejectionRate = 100*mean(rejection)),
                                            by=c("method.char","scenario","binding","missing","fixC","ar","hypo")]


res2stageS.rejection[hypo == "typeI", range(rejectionRate)]
## [1] 2.350 2.645

res2stageS.rejection[missing == FALSE & hypo == "power" & method.char == "method 3", rejectionRate]
## [1] 89.59

res2stageS.rejection[missing & ar == 1 & hypo == "power", range(rejectionRate)]
## [1] 89.92 90.72

range(res2stageS.rejection[hypo == "power" & method.char == "method 1", rejectionRate]- res2stageS.rejection[hypo == "power" & method.char == "method 2", rejectionRate])
## [1] 0.01 0.10

range(res2stageS.rejection[hypo == "power" & method.char == "method 1", rejectionRate]- res2stageS.rejection[hypo == "power" & method.char == "method 3", rejectionRate])
## [1] 0.28 0.69

## *** 3 stages
## For each run, create a binary indicator for rejection for efficacy
res3stage.rejection <- res3stage[type != "interim",.(n.stage = .N, rejection = "efficacy" %in% decision),
                                 by = c("method.char","seed","scenario","missing","binding","fixC","ar","hypo")]

## Average over runs and method within scenario
res3stageS.rejection <- res3stage.rejection[,.(n.sim = .N, rejectionRate = 100*mean(rejection)),
                                            by=c("method.char","scenario","binding","missing","fixC","ar","hypo")]


res3stageS.rejection[hypo == "typeI", range(rejectionRate)]
## [1] 2.030 2.625

res3stageS.rejection[missing == FALSE & hypo == "power" & method.char == "method 3", rejectionRate]
## [1] 89.09

res3stageS.rejection[missing & hypo == "power", range(rejectionRate)]
## [1] 89.63 90.98

range(res3stageS.rejection[hypo == "power" & method.char == "method 1", rejectionRate]- res3stageS.rejection[hypo == "power" & method.char == "method 2", rejectionRate])
## [1] 0.045 0.075

range(res3stageS.rejection[hypo == "power" & method.char == "method 1", rejectionRate]- res3stageS.rejection[hypo == "power" & method.char == "method 3", rejectionRate])
## [1] 0.645 0.865

## ** coverage
res2stage.coverage <- res2stage[decision %in% c("futility","efficacy"),
                                .(N = .N,
                                  lowerNA = 100*mean(is.na(lower_MUE)),
                                  upperNA = 100*mean(is.na(upper_MUE)),
                                  coverage = 100*mean( (lower_MUE <= truth) & (truth <= upper_MUE), na.rm=TRUE)),
                                by = c("hypo","method.char","scenario","binding","fixC","missing","ar")]

res2stage.coverage[,.(lowerNA = paste0("[",min(lowerNA),"; ",max(lowerNA),"]"), upperNA = paste0("[",min(upperNA),"; ",max(upperNA),"]")), by = c("binding")]
##    binding        lowerNA       upperNA
##     <lgcl>         <char>        <char>
## 1:    TRUE         [0; 0]        [0; 0]
## 2:   FALSE [0.045; 3.695] [0.045; 3.65]
range(res2stage.coverage$coverage,na.rm=TRUE)
## [1] 91.84 97.17
range(res2stage.coverage[hypo=="power",coverage],na.rm=TRUE)
## [1] 94.54 97.17
res2stage[type != "interim" & is.na(lower_MUE), .(.N, minp.value_MUE = min(p.value_MUE)), by = c("binding","decision","stage")]
##    binding decision stage     N minp.value_MUE
##     <lgcl>   <char> <num> <int>          <num>
## 1:   FALSE futility     1 12315         0.9783
res2stage[type != "interim" & is.na(upper_MUE), .(.N, minp.value_MUE = min(p.value_MUE)), by = c("binding","decision","stage")]
## Empty data.table (0 rows and 5 cols): binding,decision,stage,N,minp.value_MUE

res3stage.coverage <- res3stage[decision %in% c("futility","efficacy"),
                                .(N = .N,
                                  lowerNA = 100*mean(is.na(lower_MUE)),
                                  upperNA = 100*mean(is.na(upper_MUE)),
                                  coverage = 100*mean( (lower_MUE <= truth) & (truth <= upper_MUE), na.rm=TRUE)),
                                by = c("hypo","method.char","scenario","binding","fixC","missing","ar")]

res3stage.coverage[,.(lowerNA = paste0("[",min(lowerNA),"; ",max(lowerNA),"]"), upperNA = paste0("[",min(upperNA),"; ",max(upperNA),"]")), by = c("binding")]
##    binding       lowerNA       upperNA
##     <lgcl>        <char>        <char>
## 1:    TRUE        [0; 0]        [0; 0]
## 2:   FALSE [0.07; 5.125] [0.07; 5.085]
range(res3stage.coverage$coverage,na.rm=TRUE)
## [1] 89.52 97.89
range(res3stage.coverage[hypo=="power",coverage],na.rm=TRUE)
## [1] 94.50 97.89
res3stage[type != "interim" & is.na(lower_MUE), .(.N, minp.value_MUE = min(p.value_MUE)), by = c("binding","decision","stage")]
##    binding decision stage     N minp.value_MUE
##     <lgcl>   <char> <num> <int>          <num>
## 1:   FALSE futility     2  6523         0.9831
## 2:   FALSE futility     1  5065         0.9874
res3stage[type != "interim" & is.na(upper_MUE), .(.N, minp.value_MUE = min(p.value_MUE)), by = c("binding","decision","stage")]
##    binding decision stage     N minp.value_MUE
##     <lgcl>   <char> <num> <int>          <num>
## 1:   FALSE futility     2  6439         0.9913
## 2:   FALSE futility     1  5051         0.9954
res3stage[type != "interim", table(lower = is.na(lower_MUE),upper = is.na(upper_MUE))]
##        upper
## lower     FALSE    TRUE
##   FALSE 1068412       0
##   TRUE       98   11490

## ** reversal

## *** 2 stages
res2stage[,decision2 := ifelse(type=="interim",reason, decision)]
res2stage.reversal <- dcast(res2stage, seed+method+scenario+missing+binding+fixC+ar+hypo ~ type, value.var = "decision2")

res2stage.reversal[, reversal.fu2eff := (interim == "futility")*(decision == "efficacy")]
res2stage.reversal[is.na(reversal.fu2eff), reversal.fu2eff:=0]
res2stage.reversal[, reversal.eff2fu := (interim == "efficacy")*(decision == "futility")]
res2stage.reversal[is.na(reversal.eff2fu), reversal.eff2fu:=0]
res2stage.reversalS <- res2stage.reversal[, .(.N,fu2eff = 100*mean(reversal.fu2eff), eff2fu = 100*mean(reversal.eff2fu)), by = c("method","scenario","missing","binding","fixC","ar","hypo")]

res2stage.reversalS[,.(fu2eff = paste0("min=",min(fu2eff),", median = ", median(fu2eff),", max = ",max(fu2eff)),
                       eff2fu = paste0("min=",min(eff2fu),", median = ", median(eff2fu),", max = ",max(eff2fu))),
                       by = "method"]
##    method                              fu2eff                                  eff2fu
##     <int>                              <char>                                  <char>
## 1:      1  min=0, median = 0.1225, max = 0.75 min=0.045, median = 0.1425, max = 0.495
## 2:      2 min=0, median = 0.1075, max = 0.875  min=0.045, median = 0.145, max = 0.485
## 3:      3          min=0, median = 0, max = 0   min=0.2, median = 0.5675, max = 1.015

## *** 3 stages
res3stage[,decision2 := ifelse(type=="interim",reason, decision)]
res3stage.reversal <- dcast(res3stage, seed+method+scenario+missing+binding+fixC+ar+hypo ~ paste0(type,stage), value.var = "decision2")

res3stage.reversal[, c("reversal.fu2eff","reversal.eff2fu"):=0]
res3stage.reversal[(interim1 == "futility")*(decision1 == "efficacy") | (interim2 == "futility")*(decision2 == "efficacy"), reversal.fu2eff := 1]
res3stage.reversal[(interim1 == "efficacy")*(decision1 == "futility") | (interim2 == "efficacy")*(decision2 == "futility"), reversal.eff2fu := 1]
res3stage.reversalS <- res3stage.reversal[, .(.N,fu2eff = 100*mean(reversal.fu2eff), eff2fu = 100*mean(reversal.eff2fu)), by = c("method","scenario","missing","binding","fixC","ar","hypo")]

res3stage.reversalS[,.(fu2eff = paste0("min=",min(fu2eff),", median = ", median(fu2eff),", max = ",max(fu2eff)),
                       eff2fu = paste0("min=",min(eff2fu),", median = ", median(eff2fu),", max = ",max(eff2fu))),
                    by = "method"]
##    method                            fu2eff                               eff2fu
##     <int>                            <char>                               <char>
## 1:      1 min=0, median = 0.145, max = 0.875 min=0.07, median = 0.2475, max = 0.71
## 2:      2  min=0, median = 0.145, max = 1.06 min=0.07, median = 0.2725, max = 0.72
## 3:      3         min=0, median = 0, max = 0 min=0.285, median = 0.835, max = 1.13

## ** Rejection below 1.96
table2stage.below196 <- res2stage[type %in% c("decision","final"),
                                  .(.N, rejection = mean2pc(decision=="efficacy"), rejectionBelow196 = mean2pc((statistic<qnorm(0.975))*(decision=="efficacy"))),
                                  by = c("scenario","missing","method","binding","fixC","ar","hypo")]
table2stage.below196[method %in% 1:2 & fixC == FALSE,range(rejectionBelow196)]
## [1] "0.07%" "0.59%"


table3stage.below196 <- res3stage[type %in% c("decision","final"),
                            .(.N, rejection = mean2pc(decision=="efficacy"), rejectionBelow196 = mean2pc((statistic<qnorm(0.975))*(decision=="efficacy"))),
                            by = c("scenario","missing","method","binding","fixC","ar","hypo")]
table3stage.below196[method %in% 1:2 & fixC == FALSE,range(rejectionBelow196)]
## [1]"0.07%" "0.73%"

## ** Acceptance above 1.96
table2stage.above196 <- res2stage[type %in% c("decision","final"),
                                  .(.N, rejection = mean2pc(decision=="futility"), rejectionAbove196 = mean2pc((statistic>qnorm(0.975))*(decision=="futility"))),
                                  by = c("scenario","missing","method","binding","fixC","ar","hypo")]
table2stage.above196[method == 3 | fixC,range(rejectionAbove196)]
## [1] "0.17%" "0.98%"


table3stage.above196 <- res3stage[type %in% c("decision","final"),
                            .(.N, rejection = mean2pc(decision=="futility"), rejectionAbove196 = mean2pc((statistic>qnorm(0.975))*(decision=="futility"))),
                            by = c("scenario","missing","method","binding","fixC","ar","hypo")]
table3stage.above196[method == 3 | fixC,range(rejectionAbove196)]
## [1] "0.18%" "1.36%"


## ** Special cases
res2stage[, reasonNA := ifelse(is.na(reason),"NA",reason)]
res3stage[, reasonNA := ifelse(is.na(reason),"NA",reason)]
normal.case <- c("efficacy","futility","no boundary crossed","NA")

table2stage.special <- ftable(reason = res2stage[reasonNA %in% normal.case == FALSE,reason],
                              method = res2stage[reasonNA %in% normal.case == FALSE,method],
                              scenario = res2stage[reasonNA %in% normal.case == FALSE,scenario])
table2stage.special
##                                     scenario   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
## reason                       method                                                                                         
## decreasing information       1                31  46   0   0  31  46   0   0  41   3   0   0  41  41  43   3   0   0  24  32
##                              2                31  46   0   0  31  45   0   0  41   3   0   0  41  41  43   3   0   0  24  32
##                              3                32  45   0   0  32  45   0   0  42   3   0   0  42  42  45   3   0   0  25  31
## Imax reached                 1                 8   8 245 245   8   8 245 245   7   7 203 203   7   0  41   7 203 203   9   9
##                              2                 8   8 250 250   6   6 221 221   5   5 182 182   7   0  41   7 203 203   9   9
##                              3                 7   7 233 233   7   7 233 233   7   7 200 200   7   0  40   7 200 200   9   9
## stop for futility at interim 1                 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
##                              2                 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
##                              3                 8   1  62   9   8   1  62   9   4   0  59   0   4  10   1   0  59   0  11   1
100*range(as.double(table2stage.special[1:3,]))/res2stage[scenario == 1 & method == 1 & type == "interim",.N]
## [1]  0.00 0.23
100*range(as.double(table2stage.special[4:6,]))/res2stage[scenario == 1 & method == 1 & type == "interim",.N]
## [1] 0.00 1.25

table3stage.special <- ftable(reason = res3stage[reasonNA %in% normal.case == FALSE,reason],
                              method = res3stage[reasonNA %in% normal.case == FALSE,method],
                              scenario = res3stage[reasonNA %in% normal.case == FALSE,scenario])
table3stage.special
## reason                       method                                                                                                   
## decreasing information       1                 45   47    0    0   45   47    0    0   31    9    0    0   31    9    0    0   32   23
##                              2                 45   47    0    0   45   44    0    0   33    7    0    0   31    9    0    0   32   23
##                              3                 42   43    0    0   42   43    0    0   32   10    0    0   32   10    0    0   32   22
## Imax reached                 1                 60   30  727  500   60   30  727  500   40   97  641 1279   40   97  641 1279   33   23
##                              2                 60   31  742  485   48   30  639  451   36   85  573 1139   40   99  643 1293   34   23
##                              3                 56   32  650  522   56   32  650  522   36   96  584 1276   36   96  584 1276   31   24
## stop for futility at interim 1                  0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
##                              2                  0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
##                              3                 17    2   66   12   17    2   66   12   18    0   80    0   18    0   80    0   24    2
100*range(as.double(table3stage.special[1:3,]))/res3stage[scenario == 1 & method == 1 & stage == 1 & type == "interim",.N]
## [1] 0.000 0.235
100*range(as.double(table3stage.special[4:6,]))/res3stage[scenario == 1 & method == 1 & stage == 1  & type == "interim",.N]
## [1] 0.115 6.465

## ** coherence
## *** 2 stages
res2stage.PmismatchEFF <- res2stage[type != "interim",.(N = .N,
                                                           mismatchP = 100*mean(decision=="efficacy" & p.value_MUE>0.025),
                                                           mismatchCI = 100*mean(decision=="efficacy" & lower_MUE<0, na.rm=TRUE)),
                                    by = c("method.char","scenario","missing","binding","fixC","ar","hypo")]
res2stage.PmismatchEFF[,.(max(mismatchP),max(mismatchCI))]
## 1:  0  0

res2stage.PmismatchFU <- res2stage[type != "interim",.(N = .N,
                                                       mismatchP = 100*mean(decision=="futility" & p.value_MUE<0.025),
                                                       mismatchCI = 100*mean(decision=="futility" & lower_MUE>0, na.rm=TRUE)),
                                   by = c("method.char","scenario","missing","binding","fixC","ar","hypo")]
res2stage.PmismatchFU[,.(max(mismatchP),max(mismatchCI))]
## 1:     0     0

## *** 3 stages
res3stage.mismatch <- res3stage[(decision=="efficacy" & (p.value_MUE>0.025 | lower_MUE<0)) | (decision=="futility" & (p.value_MUE<0.025 | lower_MUE>0))]
NROW(res3stage.mismatch)
## [1] 24
print(res3stage.mismatch[type == "final", .(statistic = statistic-ck, p.value = p.value_MUE, lower = lower_MUE, upper = upper_MUE)],digits = 5)
##      statistic  p.value       lower  upper
##          <num>    <num>       <num>  <num>
## 1:  2.4897e-05 0.025005  1.9565e-05 1.1936
## 2:  6.6121e-05 0.025004 -1.3442e-05 1.2674
## 3:  2.4897e-05 0.025002 -1.1246e-06 1.1936
## 4: -1.0330e-04 0.025002  1.5304e-05 1.1901
## 5: -1.0330e-04 0.025002  1.5304e-05 1.1901
## 6:  8.0310e-05 0.025001  1.2907e-05 1.2013

res3stage.mismatch[type != "final", .(.N, method = unique(method), stage = unique(stage),
                                      statVSbound = paste("[",round(min(statistic-ck),5),";",round(max(statistic-ck),5),"]"),
                                      p.value = paste("[",round(min(p.value_MUE),5),";",round(max(p.value_MUE),5),"]")
                                      )]
##        N method stage          statVSbound               p.value
##    <int>  <int> <num>               <char>                <char>
## 1:    18      3     2 [ 0.00067 ; 0.0283 ] [ 0.02508 ; 0.02628 ]

res3stage.mismatch3 <- res3stage[method == 3][paste(file,seed,sep="_:_") %in% res3stage.mismatch[type != "final", paste(file,seed,sep="_:_")]]
NROW(res3stage.mismatch3)/18 ## number of line per run
## [1] 3
res3stage.mismatch3[stage==2 & type == "interim",unique(reason)]
## [1] "Imax reached"

## ** bias
true_eff <- 1

## *** two stage
res2stage[, truth := c(0,true_eff)[(hypo=="power")+1]]
res2stage.bias <- res2stage[decision %in% c("futility","efficacy"),
                            .(N = .N,
                              bias_MLE = estimate_ML-truth,
                              bias_MUE = estimate_MUE-truth,
                              mbias_MLE = (estimate_ML>truth) - 0.5,
                              mbias_MUE = (estimate_MUE>truth) - 0.5),
                            by = c("method","scenario","seed","missing","binding","fixC","ar","hypo")]
all(res2stage.bias$N==1)

res2stageS.bias <- res2stage.bias[,.(N = .N,
                                     bias_MLE = mean(bias_MLE, na.rm = TRUE),
                                     bias_MUE = mean(bias_MUE, na.rm = TRUE),
                                     mbias_MLE = mean(mbias_MLE, na.rm = TRUE),
                                     mbias_MUE = mean(mbias_MUE, na.rm = TRUE)),
                                  by=c("method","scenario","missing","binding","fixC","ar","hypo")]

res2stageS.biasS <- res2stageS.bias[method %in% 1:2 & binding == TRUE & fixC == FALSE,
                                    .(N, missing, ar, hypo, bias_MLE,bias_MUE,bias_ratio=abs(bias_MLE)/abs(bias_MUE),mbias_MLE,mbias_MUE,mbias_ratio=abs(mbias_MLE)/abs(mbias_MUE)),
                                    by = c("method","scenario")]
res2stageS.biasS


rbind(bias_MLE = quantile(abs(res2stageS.biasS$bias_MLE)),
      bias_MUE = quantile(abs(res2stageS.biasS$bias_MUE)),
      ratio = quantile(res2stageS.biasS$bias_ratio),
      mbias_MLE = 100*quantile(abs(res2stageS.biasS$mbias_MLE)),
      mbias_MUE = 100*quantile(abs(res2stageS.biasS$mbias_MUE)))
##                 0%      25%     50%     75%    100%
## bias_MLE  0.017139 0.018943 0.03560 0.03644 0.03993
## bias_MUE  0.002606 0.004265 0.01428 0.01533 0.01623
## ratio     2.166379 2.475683 2.60432 4.08294 7.31607
## mbias_MLE 1.340000 2.168750 3.25750 3.89125 4.20000
## mbias_MUE 0.335000 0.343750 0.67500 0.90125 1.03500

res2stage[type != "interim" & is.na(estimate_MUE), .(.N, minp.value_MUE = min(p.value_MUE)), by = c("binding","decision","stage")]
##    binding decision stage     N minp.value_MUE
##     <lgcl>   <char> <num> <int>          <num>
## 1:   FALSE futility     1 12315         0.9783

## *** three stage
res3stage[, truth := c(0,true_eff)[(hypo=="power")+1]]
res3stage.bias <- res3stage[decision %in% c("futility","efficacy"),
                            .(N = .N,
                              bias_MLE = estimate_ML-truth,
                              bias_MUE = estimate_MUE-truth,
                              mbias_MLE = (estimate_ML>truth) - 0.5,
                              mbias_MUE = (estimate_MUE>truth) - 0.5),
                            by = c("method","scenario","seed","missing","binding","fixC","ar","hypo")]
all(res3stage.bias$N==1)

res3stageS.bias <- res3stage.bias[,.(N = .N,
                                     bias_MLE = mean(bias_MLE, na.rm = TRUE),
                                     bias_MUE = mean(bias_MUE, na.rm = TRUE),
                                     mbias_MLE = mean(mbias_MLE, na.rm = TRUE),
                                     mbias_MUE = mean(mbias_MUE, na.rm = TRUE)),
                                  by=c("method","scenario","missing","binding","fixC","ar","hypo")]

res3stageS.biasS <- res3stageS.bias[method %in% 1:2 & binding == TRUE & fixC == FALSE,
                                    .(N, missing, ar, hypo, bias_MLE,bias_MUE,bias_ratio=abs(bias_MLE)/abs(bias_MUE),mbias_MLE,mbias_MUE,mbias_ratio=abs(mbias_MLE)/abs(mbias_MUE)),
                                    by = c("method","scenario")]

rbind(bias_MLE = quantile(abs(res3stageS.biasS$bias_MLE)),
      bias_MUE = quantile(abs(res3stageS.biasS$bias_MUE)),
      ratio = quantile(res3stageS.biasS$bias_ratio),
      mbias_MLE = 100*quantile(abs(res3stageS.biasS$mbias_MLE)),
      mbias_MUE = 100*quantile(abs(res3stageS.biasS$mbias_MUE)))
##                0%     25%     50%     75%    100%
## bias_MLE  0.02932 0.04256 0.05279 0.06484 0.06981
## bias_MUE  0.01560 0.02064 0.02743 0.03417 0.03643
## ratio     1.85978 1.87851 1.90747 1.98972 2.07678
## mbias_MLE 2.74000 3.54250 4.88000 5.05875 5.42500
## mbias_MUE 0.08500 0.33875 0.36750 0.43125 0.45000

res3stage[type != "interim" & is.na(estimate_MUE), .(.N, minp.value_MUE = min(p.value_MUE)), by = c("binding","decision","stage")]
##    binding decision stage     N minp.value_MUE
##     <lgcl>   <char> <num> <int>          <num>
## 1:   FALSE futility     2  6523         0.9831
## 2:   FALSE futility     1  5065         0.9874

res3stage[type != "interim", table(is.na(estimate_MUE), is.na(lower_MUE))]
##         FALSE    TRUE
## FALSE 1068412       0
## TRUE        0   11588

## * Appendix

## ** information
## *** 2 stages
res2stage.info1 <- res2stage[ar == 1 & infoBias == 0 & missing == TRUE, .(n.interim = sum(type=="interim"), infoPC.interim = .SD[type=="interim", mean(infoPC)],
                                                          n.decision = sum(type=="decision"), infoPC.decision = .SD[type=="decision", mean(infoPC)],
                                                          n.final = sum(type=="final"), infoPC.final = .SD[type=="final", mean(infoPC)]),
                             by = c("method","binding","fixC","ar","scenario")]
rbind(interim = quantile(res2stage.info1$infoPC.interim),
      decision = quantile(res2stage.info1$infoPC.decision),
      final = quantile(res2stage.info1$infoPC.final))
##              0%    25%    50%    75%   100%
## interim  0.5404 0.5429 0.5435 0.5450 0.5467
## decision 0.6720 0.6774 0.6781 0.6827 0.7163
## final    1.0114 1.0177 1.0208 1.0217 1.0359

res2stage.info2 <- res2stage[ar == 2 & infoBias == 0 & missing == TRUE, .(n.interim = sum(type=="interim"), infoPC.interim = .SD[type=="interim", mean(infoPC)],
                                                          n.decision = sum(type=="decision"), infoPC.decision = .SD[type=="decision", mean(infoPC)],
                                                          n.final = sum(type=="final"), infoPC.final = .SD[type=="final", mean(infoPC)]),
                             by = c("method","binding","fixC","ar","scenario")]
rbind(interim = quantile(res2stage.info2$infoPC.interim),
      decision = quantile(res2stage.info2$infoPC.decision),
      final = quantile(res2stage.info2$infoPC.final))
##              0%    25%    50%    75%   100%
## interim  0.5617 0.5643 0.5650 0.5669 0.5683
## decision 0.8176 0.8219 0.8233 0.8492 0.9195
## final    1.0090 1.0148 1.0178 1.0220 1.0342

## *** 3 stages
res3stage.info1 <- res3stage[ar == 1 & missing == TRUE, .(n.interim1 = sum(stage == 1 & type=="interim"), infoPC.interim1 = .SD[stage == 1 & type=="interim", mean(infoPC)],
                                                          n.decision1 = sum(stage == 1 & type=="decision"), infoPC.decision1 = .SD[stage == 1 & type=="decision", mean(infoPC)],
                                                          n.interim2 = sum(stage == 2 & type=="interim"), infoPC.interim2 = .SD[stage == 2 & type=="interim", mean(infoPC)],
                                                          n.decision2 = sum(stage == 2 & type=="decision"), infoPC.decision2 = .SD[stage == 2 & type=="decision", mean(infoPC)],
                                                          n.final = sum(type=="final"), infoPC.final = .SD[type=="final", mean(infoPC)]),
                             by = c("method","binding","fixC","ar","scenario")]
rbind(interim1 = quantile(res3stage.info1$infoPC.interim1),
      decision1 = quantile(res3stage.info1$infoPC.decision1),
      interim2 = quantile(res3stage.info1$infoPC.interim2),
      decision2 = quantile(res3stage.info1$infoPC.decision2),
      final = quantile(res3stage.info1$infoPC.final))
##               0%    25%    50%    75%   100%
## interim1  0.3846 0.3872 0.3874 0.3895 0.3906
## decision1 0.5185 0.5236 0.5261 0.5317 0.5449
## interim2  0.6304 0.6352 0.6381 0.6401 0.6473
## decision2 0.7568 0.7633 0.7679 0.7983 0.8887
## final     1.0006 1.0087 1.0114 1.0180 1.0329

res3stage.info2 <- res3stage[ar == 2 & missing == TRUE, .(n.interim1 = sum(stage == 1 & type=="interim"), infoPC.interim1 = .SD[stage == 1 & type=="interim", mean(infoPC)],
                                                          n.decision1 = sum(stage == 1 & type=="decision"), infoPC.decision1 = .SD[stage == 1 & type=="decision", mean(infoPC)],
                                                          n.interim2 = sum(stage == 2 & type=="interim"), infoPC.interim2 = .SD[stage == 2 & type=="interim", mean(infoPC)],
                                                          n.decision2 = sum(stage == 2 & type=="decision"), infoPC.decision2 = .SD[stage == 2 & type=="decision", mean(infoPC)],
                                                          n.final = sum(type=="final"), infoPC.final = .SD[type=="final", mean(infoPC)]),
                             by = c("method","binding","fixC","ar","scenario")]
rbind(interim1 = quantile(res3stage.info2$infoPC.interim1),
      decision1 = quantile(res3stage.info2$infoPC.decision1),
      interim2 = quantile(res3stage.info2$infoPC.interim2),
      decision2 = quantile(res3stage.info2$infoPC.decision2),
      final = quantile(res3stage.info2$infoPC.final))
##               0%    25%    50%    75%   100%
## interim1  0.4039 0.4067 0.4069 0.4093 0.4105
## decision1 0.6599 0.6648 0.6658 0.6749 0.6947
## interim2  0.6508 0.6562 0.6585 0.6613 0.6686
## decision2 0.9011 0.9084 0.9117 0.9498 1.0482
## final     0.9935 1.0030 1.0051 1.0108 1.0241


## * Review
## ** ck in table 1

res2stage[method %in% 1:2 & fixC == FALSE & missing == TRUE & ar == 1 & type == "decision" & infoBias == 0,
          .(n.sim = .N, ckBelow1.96 = 100*mean(ck<qnorm(0.975)), statisticBelow1.96 = 100*mean(statistic < qnorm(0.975))),
          by = c("method.char","hypo","binding")]
##    method.char   hypo binding n.sim ckBelow1.96 statisticBelow1.96
##         <char> <char>  <lgcl> <int>       <num>              <num>
## 1:    method 1  power    TRUE 10625       99.92              7.501
## 2:    method 2  power    TRUE 10632       99.90              7.543
## 3:    method 1  typeI    TRUE 14264       99.94             98.829
## 4:    method 2  typeI    TRUE 14293       99.92             98.832
## 5:    method 1  power   FALSE 10596       99.93              6.663
## 6:    method 2  power   FALSE 10605       99.90              6.723
## 7:    method 1  typeI   FALSE   201       96.52             15.423
## 8:    method 2  typeI   FALSE   201       96.02             15.423

res2stage[method == 1 & fixC == FALSE & missing == TRUE & ar == 1 & type == "decision" & hypo == "typeI" & binding == TRUE & infoBias == 0,hist(ck)]
res2stage[method == 1 & fixC == FALSE & missing == TRUE & ar == 1 & type == "decision" & hypo == "typeI" & binding == TRUE & infoBias == 0,quantile(ck, c(0,0.01,0.25,0.5,0.75,1), na.rm=TRUE)]
##     0%     1%    25%    50%    75%   100% 
## 0.8417 1.4545 1.5726 1.6249 1.6788 1.9600 

res2stage[method == 1 & fixC == FALSE & missing == TRUE & ar == 1 & type == "decision" & hypo == "typeI" & binding == FALSE & infoBias == 0, hist(ck)]
res2stage[method == 1 & fixC == FALSE & missing == TRUE & ar == 1 & type == "decision" & hypo == "typeI" & binding == FALSE & infoBias == 0,quantile(ck, c(0,0.01,0.25,0.5,0.75,1))]
##    0%    1%   25%   50%   75%  100% 
## 1.468 1.500 1.609 1.664 1.730 1.960 


res2stage[method %in% 1:2 & fixC == FALSE & missing == TRUE & ar == 1 & type != "interim" & infoBias == 0,
          .(n.sim = .N, ckBelow1.96 = 100*mean(ck<qnorm(0.975)), statisticBelow1.96 = 100*mean(statistic < qnorm(0.975))),
          by = c("method.char","hypo","binding")]
##    method.char   hypo binding n.sim ckBelow1.96 statisticBelow1.96
##         <char> <char>  <lgcl> <int>       <num>              <num>
## 1:    method 1  power    TRUE 20000      53.085              9.650
## 2:    method 2  power    TRUE 20000      53.120              9.655
## 3:    method 1  typeI    TRUE 20000      71.280             97.295
## 4:    method 2  typeI    TRUE 20000      71.430             97.295
## 5:    method 1  power   FALSE 20000      52.945              8.870
## 6:    method 2  power   FALSE 20000      52.970              8.880
## 7:    method 1  typeI   FALSE 20000       0.970             97.110
## 8:    method 2  typeI   FALSE 20000       0.965             97.110

res2stage[method == 1 & fixC == FALSE & missing == TRUE & ar == 1 & type != "interim" & hypo == "typeI" & binding == TRUE & infoBias == 0,hist(ck)]
res2stage[method == 1 & fixC == FALSE & missing == TRUE & ar == 1 & type != "interim" & hypo == "typeI" & binding == TRUE & infoBias == 0,quantile(ck, c(0,0.01,0.25,0.5,0.75,1))]
##     0%     1%    25%    50%    75%   100% 
## 0.8417 1.4638 1.5949 1.6671 1.9943 2.0299 

res2stage[method == 1 & fixC == FALSE & missing == TRUE & ar == 1 & type != "interim" & hypo == "typeI" & binding == FALSE & infoBias == 0, hist(ck)]
res2stage[method == 1 & fixC == FALSE & missing == TRUE & ar == 1 & type != "interim" & hypo == "typeI" & binding == FALSE & infoBias == 0,quantile(ck, c(0,0.01,0.25,0.5,0.75,1))]
##    0%    1%   25%   50%   75%  100% 
## 1.468 1.960 2.015 2.028 2.043 2.250

##----------------------------------------------------------------------
### TEXT-simulation.R ends here
