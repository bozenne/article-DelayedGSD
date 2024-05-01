### table1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 19 2023 (10:24) 
## Version: 
## Last-Updated: maj  1 2024 (11:09) 
##           By: Brice Ozenne
##     Update #: 104
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

##' @description Generate table
createTableResSim <- function(data, xtable, sep = " / "){

    ## prepare data
    ## select the last interim
    data.interim <- data[type == "interim", .(decision.interim=paste0(decision[which.max(stage)]," (",reason[which.max(stage)],")")),
                         by = c("scenario","seed","method")]
    data.decision <- merge(x = data[decision %in% c("futility","efficacy")],
                           y = data.interim, by = c("scenario","method","seed"),
                           all = TRUE)
    data.decision[, method.char := factor(method.char, levels = c("method 1","method 1 fixC","method 2","method 2 fixC","method 3"))]
    setkeyv(data.decision, "method.char")

    ls.sumstat <- list()
    ls.sumstat$type1 <- data.decision[hypo=="typeI",.(statistic = "type 1 error",
                                                      n = .N,
                                                      ## freq = sum(decision=="efficacy"),
                                                      value = mean2pc(decision=="efficacy")),
                                      by="method.char"]
    ls.sumstat$power <- data.decision[hypo=="power",.(statistic = "power",
                                                      n = .N,
                                                      ## freq = sum(decision=="efficacy"),
                                                      value = mean2pc(decision=="efficacy")),
                                      by="method.char"]
    ls.sumstat$CINA <- data.decision[hypo=="power",.(statistic = "CI [NA,NA]",
                                                     n = .N,
                                                     ## freq = sum(is.na(lower_MUE) | is.na(upper_MUE)),
                                                     value = mean2pc(is.na(lower_MUE) | is.na(upper_MUE))),
                                     by="method.char"]
    ls.sumstat$coverage <- data.decision[hypo=="power" & !is.na(lower_MUE) & !is.na(upper_MUE),.(statistic = "coverage",
                                                                                                 n = .N,
                                                                                                 ## freq = sum(lower_MUE <= truth & truth <= upper_MUE),
                                                                                                 value = mean2pc(lower_MUE <= truth & truth <= upper_MUE)),
                                         by="method.char"]
    ls.sumstat$reversal <- data.decision[hypo=="power",.(statistic = "reversal",
                                                         n = .N,
                                                         ## freq = paste(sum(decision.interim=="stop (futility)" & decision=="efficacy"),
                                                         ##              sum(decision.interim=="stop (efficacy)" & decision=="futility"),
                                                         ##              sep="/"),
                                                         value = paste(mean2pc(decision.interim=="stop (futility)" & decision=="efficacy"),
                                                                       mean2pc(decision.interim=="stop (efficacy)" & decision=="futility"),
                                                                       sep="/")
                                                         ),by="method.char"]
    ls.sumstat$abnormal <- data.decision[hypo=="power",.(statistic = "abnormal",
                                                         n = .N,
                                                         ## freq = paste(sum(decision=="efficacy" & statistic < 1.96),
                                                         ##              sum(decision=="futility" & statistic >= ck),
                                                         ##              sep="/"),
                                                         value = paste(mean2pc(decision=="efficacy" & statistic < 1.96),
                                                                       mean2pc(decision=="futility" & statistic >= ck),
                                                                       sep="/")),
                                         by="method.char"]
    ls.sumstat$LMMEmeanB <- data.decision[!is.na(estimate_MUE),.(statistic = "mean bias LMME",
                                                                 n = .N,
                                                                 value = paste(mean2num(.SD[hypo=="typeI",estimate_ML]-.SD[hypo=="typeI",truth]),
                                                                               mean2num(.SD[hypo=="power",estimate_ML]-.SD[hypo=="power",truth]),
                                                                               sep=sep)),
                                          by="method.char"]
    ls.sumstat$MUEmeanB <- data.decision[!is.na(estimate_MUE),.(statistic = "mean bias MUE",
                                                                n = .N,
                                                                value = paste(mean2num(.SD[hypo=="typeI",estimate_MUE]-.SD[hypo=="typeI",truth]),
                                                                              mean2num(.SD[hypo=="power",estimate_MUE]-.SD[hypo=="power",truth]),
                                                                              sep=sep)),
                                         by="method.char"]
    ls.sumstat$LMMEmedianB <- data.decision[!is.na(estimate_MUE),.(statistic = "median bias LMME",
                                                                   n = .N,
                                                                   value = paste(mean2pc((.SD[hypo=="typeI",estimate_ML]>.SD[hypo=="typeI",truth])-0.5),
                                                                                 mean2pc((.SD[hypo=="power",estimate_ML]>.SD[hypo=="power",truth])-0.5),
                                                                                 sep=sep)),
                                            by="method.char"]
    ls.sumstat$MUEmedianB <- data.decision[!is.na(estimate_MUE),.(statistic = "median bias MUE",
                                                                  n = .N,
                                                                  value = paste(mean2pc((.SD[hypo=="typeI",estimate_MUE]>.SD[hypo=="typeI",truth])-0.5),
                                                                                mean2pc((.SD[hypo=="power",estimate_MUE]>.SD[hypo=="power",truth])-0.5),
                                                                                sep=sep)),
                                           by="method.char"]

    dtL.sumstat <- do.call(rbind,ls.sumstat)
    dtL.sumstat[, statistic := factor(statistic, unique(statistic))]
    dtW.sumstat <- dcast(dtL.sumstat, formula = statistic~method.char, value.var = "value")


    if(xtable){
        ## % -> \\%
        dtW.sumstat[[2]] <- gsub("%","\\%",dtW.sumstat[[2]], fixed = TRUE)
        dtW.sumstat[[3]] <- gsub("%","\\%",dtW.sumstat[[3]], fixed = TRUE)
        dtW.sumstat[[4]] <- gsub("%","\\%",dtW.sumstat[[4]], fixed = TRUE)
        dtW.sumstat[[5]] <- gsub("%","\\%",dtW.sumstat[[5]], fixed = TRUE)
        dtW.sumstat[[6]] <- gsub("%","\\%",dtW.sumstat[[6]], fixed = TRUE)
        
        add <- "\\hspace{3mm}"
        dtW.sumstat$statistic <- paste0(add,c("Type 1 error",
                                              "Power",
                                              "CI=[NA;NA]",
                                              "Coverage",
                                              "Reversal(F\\(\\rightarrow\\)E/E\\(\\rightarrow\\)F)",
                                              "E\\(\\left(\\tilde{Z}_k<1.96\\right)\\)/F\\(\\left(\\tilde{Z}\\geq c_k\\right)\\)",
                                              "Mean bias (\\(\\mathcal{H}_0\\)/\\(\\mathcal{H}_1\\)): LMME",
                                              "\\hphantom{Mean bias (\\(\\mathcal{H}_0\\)/\\(\\mathcal{H}_1\\))}: MUE",
                                              "Median bias (\\(\\mathcal{H}_0\\)/\\(\\mathcal{H}_1\\)): LMME",
                                              "\\hphantom{Median bias (\\(\\mathcal{H}_0\\)/\\(\\mathcal{H}_1\\))}: MUE"))
        xtable.sumstat <- xtable(cbind(dtW.sumstat[,1:3],"x"="",dtW.sumstat[,4:5],"y"="",dtW.sumstat[,6]),
                                 type='latex')        
        add.to.row <- list(pos = list(4,6), command = c("[2mm]","[2mm]"))
        print(xtable.sumstat, include.rownames = FALSE, sanitize.text.function=identity, add.to.row = add.to.row)
    }else{
        return(dtW.sumstat)
    }

}

## * 2 stages

## ** Load data
res2stage <- readRDS(file.path("Results-built","res2stage.rds"))
res2stage[, method.char := paste0("method ",method, c(""," fixC")[fixC+1])]
res2stage[, stage.char := factor(stage, 1:2, c("interim","final"))]
res2stage[, truth := ifelse(hypo=="power",1,0)]

## method 3 fixC same as method 3
res2stage.red <- res2stage[method.char != "method 3 fixC" & missing==TRUE]
## unique(res2stage.red$scenario)


## ** Generate table
keep.col <- c("scenario", "hypo", "method", "stage", "type", "statistic", "ck",
              "nX1.interim", "nX3.interim", "sample.size",
              "estimate_ML", "se_ML", "p.value_ML", "lower_ML", "upper_ML",
              "estimate_MUE", "p.value_MUE", "lower_MUE", "upper_MUE",
              "decision", "reason", "method.char", "stage.char", "truth","seed")       

## *** ar 1 binding
res2stage.ar1binding <- res2stage.red[infoBias==0 & binding==TRUE & ar==1 & !is.na(decision),.SD,.SDcols=keep.col]

## sample size
unique(unlist(res2stage.ar1binding$sample.size))
## [1] 265

## number of simulations
res2stage.ar1binding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000

## number of patients: included/with final outcome value
res2stage.ar1binding[, .(included = paste0("[",min(nX1.interim),";", median(nX1.interim), ";", max(nX1.interim),"]"),
                            incomplete = paste0("[",min(nX1.interim-nX3.interim),";",median(nX1.interim-nX3.interim),";", max(nX1.interim-nX3.interim),"]"),
                            outcome = paste0("[",min(nX3.interim),";", median(nX3.interim),";", max(nX3.interim),"]")), by = "hypo"]
##     hypo      included incomplete       outcome
## 1: power [153;163;179] [29;44;60] [109;119;128]
## 2: typeI [153;163;179] [29;44;60] [109;119;128]

## table
table2stage.ar1binding <- createTableResSim(res2stage.ar1binding, xtable = FALSE)
table2stage.ar1binding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC        method 3
##  1:     type 1 error          2.60%           2.46%          2.61%           2.47%           2.54%
##  2:            power         91.00%          90.53%         90.98%          90.67%          90.53%
##  3:       CI [NA,NA]              0               0              0               0               0
##  4:         coverage         94.66%          95.68%         94.65%          95.79%          95.26%
##  5:         reversal    0.16%/0.06%     0.06%/0.43%    0.15%/0.06%     0.05%/0.41%         0/0.75%
##  6:         abnormal        0.47%/0             0/0        0.46%/0             0/0         0/0.02%
##  7:   mean bias LMME -0.036 / 0.041  -0.036 / 0.041 -0.036 / 0.041  -0.036 / 0.041  -0.036 / 0.042
##  8:    mean bias MUE -0.011 / 0.018 -0.092 / -0.031 -0.011 / 0.018 -0.096 / -0.032 -0.038 / -0.001
##  9: median bias LMME -3.00% / 4.78%  -3.00% / 4.78% -2.98% / 4.78%  -3.16% / 4.79%  -3.21% / 4.67%
## 10:  median bias MUE 0.01% / -0.65% -8.82% / -3.70% 0.03% / -0.65% -8.75% / -3.79% -3.21% / -1.69%
createTableResSim(res2stage.ar1binding, xtable = TRUE)

res2stage.ar1binding[type %in% c("decision","final") & method.char == "method 1", .(.N, prob.below = mean((statistic<=1.96)*(decision=="efficacy"))), by = c("hypo")]
##     hypo     N prob.below
## 1: power 10000     0.5327
## 2: typeI 10000     0.7114

## *** ar 2 binding
res2stage.ar2binding <- res2stage.red[infoBias==0 & binding==TRUE & ar==2 & !is.na(decision),.SD,.SDcols=keep.col]

## sample size
unique(unlist(res2stage.ar2binding$sample.size))
## [1] 265

## number of simulations
res2stage.ar2binding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000

## number of patients: included/with final outcome value
res2stage.ar2binding[, .(included = paste0("[",min(nX1.interim),";", median(nX1.interim), ";", max(nX1.interim),"]"),
                         incomplete = paste0("[",min(nX1.interim-nX3.interim),";",median(nX1.interim-nX3.interim),";", max(nX1.interim-nX3.interim),"]"),
                         outcome = paste0("[",min(nX3.interim),";", median(nX3.interim),";", max(nX3.interim),"]")), by = "hypo"]
##     hypo      included incomplete       outcome
## 1: power [178;193;211] [58;74;93] [109;119;129]
## 2: typeI [178;193;211] [58;74;93] [109;119;129]

## table
table2stage.ar2binding <- createTableResSim(res2stage.ar2binding, xtable = FALSE)
table2stage.ar2binding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC        method 3
##  1:     type 1 error          2.56%           2.31%          2.51%           2.33%           2.41%
##  2:            power         91.40%          90.86%         91.29%          90.90%          90.74%
##  3:       CI [NA,NA]              0               0              0               0               0
##  4:         coverage         94.95%          95.52%         94.93%          95.58%          95.26%
##  5:         reversal    0.64%/0.18%     0.41%/0.49%    0.78%/0.20%     0.40%/0.50%         0/0.88%
##  6:         abnormal        0.54%/0             0/0        0.49%/0             0/0         0/0.29%
##  7:   mean bias LMME -0.015 / 0.021  -0.015 / 0.021 -0.015 / 0.021  -0.016 / 0.022  -0.016 / 0.021
##  8:    mean bias MUE  0.001 / 0.007 -0.053 / -0.027  0.001 / 0.007 -0.055 / -0.028 -0.024 / -0.014
##  9: median bias LMME -0.45% / 2.94%  -0.45% / 2.94% -0.35% / 2.93%  -0.63% / 2.98%  -0.81% / 2.75%
## 10:  median bias MUE 0.59% / -0.47% -6.89% / -2.94% 0.59% / -0.46% -6.77% / -2.98% -3.15% / -1.40%
createTableResSim(res2stage.ar2binding, xtable = TRUE)

power2stage.ar2binding <- as.numeric(gsub("%","",unlist(table2stage.ar2binding[2,.SD,.SDcols = names(table2stage.ar2binding)[-1]]),fixed=TRUE))
power2stage.ar2binding[1:2]-power2stage.ar2binding[3:4]
power2stage.ar2binding[1:2]-power2stage.ar2binding[5]
## [1]  0.11 -0.04
## [1] 0.66 0.12

## *** ar 1 non-binding
res2stage.ar1nonbinding <- res2stage.red[infoBias==0 & binding==FALSE & ar==1 & !is.na(decision),.SD,.SDcols=keep.col]

## sample size
unique(unlist(res2stage.ar1nonbinding$sample.size))
## [1] 268

## number of simulations
res2stage.ar1nonbinding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000

## number of patients: included/with final outcome value
res2stage.ar1nonbinding[, .(included = paste0("[",min(nX1.interim),";", median(nX1.interim), ";", max(nX1.interim),"]"),
                            incomplete = paste0("[",min(nX1.interim-nX3.interim),";",median(nX1.interim-nX3.interim),";", max(nX1.interim-nX3.interim),"]"),
                            outcome = paste0("[",min(nX3.interim),";", median(nX3.interim),";", max(nX3.interim),"]")), by = "hypo"]
##     hypo      included incomplete       outcome
## 1: power [154;164;178] [30;44;60] [111;120;128]
## 2: typeI [154;164;178] [30;44;60] [111;120;128]

## table
table2stage.ar1nonbinding <- createTableResSim(res2stage.ar1nonbinding, xtable = FALSE)
table2stage.ar1nonbinding
##            statistic      method 1   method 1 fixC      method 2   method 2 fixC       method 3
##  1:     type 1 error         2.57%           2.53%         2.57%           2.54%          2.60%
##  2:            power        91.20%          90.82%        91.20%          90.99%         90.90%
##  3:       CI [NA,NA]         3.12%           3.50%         3.15%           3.28%          3.48%
##  4:         coverage        95.82%          97.23%        95.82%          97.27%         96.56%
##  5:         reversal   0.21%/0.09%     0.07%/0.33%   0.22%/0.08%     0.07%/0.33%        0/0.61%
##  6:         abnormal       0.38%/0             0/0       0.39%/0             0/0        0/0.03%
##  7:   mean bias LMME 0.006 / 0.066   0.006 / 0.068 0.006 / 0.066   0.006 / 0.067  0.005 / 0.067
##  8:    mean bias MUE 0.006 / 0.045   0.005 / 0.000 0.006 / 0.045  0.005 / -0.003  0.006 / 0.029
##  9: median bias LMME 0.01% / 6.07%  -0.02% / 6.29% 0.01% / 6.09%  -0.04% / 6.18% -0.09% / 6.25%
## 10:  median bias MUE 0.01% / 0.64% -0.02% / -2.59% 0.01% / 0.65% -0.04% / -2.99%  0.03% / 0.02%
createTableResSim(res2stage.ar1nonbinding, xtable = TRUE)

power2stage.ar1nonbinding <- as.numeric(gsub("%","",unlist(table2stage.ar1nonbinding[2,.SD,.SDcols = names(table2stage.ar1nonbinding)[-1]]),fixed=TRUE))
power2stage.ar1nonbinding[1:2]-power2stage.ar1nonbinding[3:4]
power2stage.ar1nonbinding[1:2]-power2stage.ar1nonbinding[5]
## [1]  0.01 -0.11
## [1]  0.31 -0.13

## *** ar 2 non-binding
res2stage.ar2nonbinding <- res2stage.red[infoBias==0 & binding==FALSE & ar==2 & !is.na(decision),.SD,.SDcols=keep.col]

## sample size
unique(unlist(res2stage.ar2nonbinding$sample.size))
## [1] 268

## number of simulations
res2stage.ar2nonbinding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000

## number of patients: included/with final outcome value
res2stage.ar2nonbinding[, .(included = paste0("[",min(nX1.interim),";", median(nX1.interim), ";", max(nX1.interim),"]"),
                            incomplete = paste0("[",min(nX1.interim-nX3.interim),";",median(nX1.interim-nX3.interim),";", max(nX1.interim-nX3.interim),"]"),
                            outcome = paste0("[",min(nX3.interim),";", median(nX3.interim),";", max(nX3.interim),"]")), by = "hypo"]
##     hypo      included  incomplete       outcome
## 1: power [179;194;215] [55;74;100] [111;120;129]
## 2: typeI [179;194;215] [55;74;100] [111;120;129]

## table
table2stage.ar2nonbinding <- createTableResSim(res2stage.ar2nonbinding, xtable = FALSE)
table2stage.ar2nonbinding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC       method 3
##  1:     type 1 error          2.47%           2.36%          2.46%           2.37%          2.45%
##  2:            power         91.55%          91.14%         91.42%          91.22%         91.25%
##  3:       CI [NA,NA]          2.84%           3.25%          3.16%           3.07%          3.25%
##  4:         coverage         96.18%          97.14%         96.29%          97.11%         96.83%
##  5:         reversal    0.69%/0.17%     0.46%/0.35%    0.77%/0.17%     0.43%/0.35%        0/0.61%
##  6:         abnormal        0.41%/0             0/0        0.36%/0             0/0        0/0.29%
##  7:   mean bias LMME  0.004 / 0.041   0.003 / 0.042  0.004 / 0.042   0.003 / 0.042  0.002 / 0.041
##  8:    mean bias MUE  0.004 / 0.028  0.003 / -0.002  0.004 / 0.030  0.003 / -0.004  0.005 / 0.012
##  9: median bias LMME -0.08% / 3.41%  -0.13% / 3.63% -0.08% / 3.57%  -0.14% / 3.57% -0.25% / 3.28%
## 10:  median bias MUE -0.08% / 0.57% -0.13% / -1.74% -0.08% / 0.75% -0.14% / -1.85% 0.13% / -0.60%
createTableResSim(res2stage.ar2nonbinding, xtable = TRUE)

power2stage.ar2nonbinding <- as.numeric(gsub("%","",unlist(table2stage.ar2nonbinding[2,.SD,.SDcols = names(table2stage.ar2nonbinding)[-1]]),fixed=TRUE))
power2stage.ar2nonbinding[1:2]-power2stage.ar2nonbinding[3:4]
power2stage.ar2nonbinding[1:2]-power2stage.ar2nonbinding[5]
## [1]  0.13 -0.08
## [1]  0.30 -0.11

## *** information bias
res2stage.infoBias.pos <- res2stage.red[infoBias==1 & !is.na(decision),.SD,.SDcols=keep.col]

## number of simulations
res2stage.infoBias.pos[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000

## table (NA in bias because no simulation under the null)
res2stage.infoBias.pos <- createTableResSim(res2stage.infoBias.pos, xtable = FALSE)
res2stage.infoBias.pos
##           statistic    method 1    method 2   method 3
## 1:            power      91.20%      91.15%     90.71%
## 2:       CI [NA,NA]       3.10%       3.20%      3.69%
## 3:         coverage      95.82%      95.81%     96.27%
## 4:         reversal 0.21%/0.07% 0.20%/0.08%    0/1.01%
## 5:         abnormal     0.38%/0     0.35%/0        0/0
## 6:   mean bias LMME  NA / 0.066  NA / 0.067 NA / 0.068
## 7:    mean bias MUE  NA / 0.045  NA / 0.045 NA / 0.043
## 8: median bias LMME  NA / 6.06%  NA / 6.12% NA / 6.36%
## 9:  median bias MUE  NA / 0.63%  NA / 0.67% NA / 1.10%

res2stage.infoBias.neg <- res2stage.red[infoBias==-1 & !is.na(decision),.SD,.SDcols=keep.col]

## number of simulations
res2stage.infoBias.neg[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000

## table (NA in bias because no simulation under the null)
res2stage.infoBias.neg <- createTableResSim(res2stage.infoBias.neg, xtable = FALSE)
res2stage.infoBias.neg
##           statistic    method 1    method 2    method 3
## 1:            power      91.20%      91.20%      90.79%
## 2:       CI [NA,NA]       3.13%       3.14%       3.56%
## 3:         coverage      95.83%      95.83%      97.37%
## 4:         reversal 0.21%/0.09% 0.22%/0.09%     0/0.42%
## 5:         abnormal     0.38%/0     0.39%/0     0/0.07%
## 6:   mean bias LMME  NA / 0.066  NA / 0.066  NA / 0.068
## 7:    mean bias MUE  NA / 0.045  NA / 0.045 NA / -0.005
## 8: median bias LMME  NA / 6.08%  NA / 6.08%  NA / 6.30%
## 9:  median bias MUE  NA / 0.65%  NA / 0.64% NA / -3.07%

## * 3 stages

## ** Load data
res3stage <- readRDS(file.path("Results-built","res3stage.rds"))
res3stage[, method.char := paste0("method ",method, c(""," fixC")[fixC+1])]
res3stage[, stage.char := factor(stage, 1:3, c("interim1","interim2","final"))]
res3stage[, truth := ifelse(hypo=="power",1,0)]

## method 3 fixC same as method 3
res3stage.red <- res3stage[method.char != "method 3 fixC" & missing==TRUE]

## ** Generate table
keep.col <- c("scenario", "hypo", "method", "stage", "type", "statistic", "ck",
              "estimate_ML", "se_ML", "p.value_ML", "lower_ML", "upper_ML",
              "estimate_MUE", "p.value_MUE", "lower_MUE", "upper_MUE",
              "decision", "reason", "method.char", "stage.char", "truth","seed")       

## *** ar 1 binding
res3stage.ar1binding <- res3stage.red[binding==TRUE & ar==1 & !is.na(decision),.SD,.SDcols=keep.col]


## number of simulations
res3stage.ar1binding[type!="interim",.N,by=c("method.char","hypo")][,unique(N)]
## [1] 10000

table3stage.ar1binding <- createTableResSim(res3stage.ar1binding, xtable = FALSE)
table3stage.ar1binding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC        method 3
##  1:     type 1 error          2.70%           2.52%          2.71%           2.52%           2.57%
##  2:            power         90.68%          90.09%         90.59%          90.31%          89.78%
##  3:       CI [NA,NA]          0.02%           0.02%          0.02%           0.03%           0.02%
##  4:         coverage         94.77%          95.91%         94.79%          95.91%          94.93%
##  5:         reversal    0.39%/0.11%     0.16%/0.47%    0.38%/0.13%     0.13%/0.46%         0/1.06%
##  6:         abnormal        0.59%/0             0/0        0.55%/0             0/0         0/0.10%
##  7:   mean bias LMME -0.066 / 0.054  -0.066 / 0.054 -0.066 / 0.054  -0.067 / 0.055  -0.067 / 0.055
##  8:    mean bias MUE -0.033 / 0.029 -0.130 / -0.023 -0.032 / 0.028 -0.132 / -0.024  -0.063 / 0.005
##  9: median bias LMME -4.83% / 5.66%  -4.83% / 5.66% -4.82% / 5.66%  -5.03% / 5.73%  -5.07% / 5.49%
## 10:  median bias MUE 0.75% / -0.25% -4.73% / -4.66% 0.74% / -0.26% -4.98% / -4.89% -1.46% / -2.47%
createTableResSim(res3stage.ar1binding, xtable = TRUE)

## *** ar 2 binding
res3stage.ar2binding <- res3stage.red[binding==TRUE & ar==2 & !is.na(decision),.SD,.SDcols=keep.col]


## number of simulations
res3stage.ar2binding[type!="interim",.N,by=c("method.char","hypo")][,unique(N)]
## [1] 10000

table3stage.ar2binding <- createTableResSim(res3stage.ar2binding, xtable = FALSE)
table3stage.ar2binding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC        method 3
##  1:     type 1 error          2.65%           2.29%          2.63%           2.27%           2.39%
##  2:            power         91.04%          90.36%         90.95%          90.48%          90.41%
##  3:       CI [NA,NA]              0               0              0               0               0
##  4:         coverage         94.65%          95.50%         94.61%          95.61%          95.19%
##  5:         reversal    0.72%/0.22%     0.51%/0.69%    0.92%/0.27%     0.52%/0.71%         0/1.17%
##  6:         abnormal        0.68%/0             0/0        0.64%/0             0/0         0/0.30%
##  7:   mean bias LMME -0.039 / 0.031  -0.039 / 0.031 -0.039 / 0.030  -0.039 / 0.032  -0.039 / 0.033
##  8:    mean bias MUE -0.018 / 0.017 -0.084 / -0.019 -0.018 / 0.017 -0.087 / -0.020 -0.049 / -0.005
##  9: median bias LMME -3.17% / 2.95%  -3.17% / 2.96% -3.20% / 2.94%  -3.24% / 2.99%  -3.39% / 3.02%
## 10:  median bias MUE 0.49% / -0.28% -3.06% / -2.89% 0.53% / -0.32% -3.22% / -2.99% -1.44% / -2.67%
createTableResSim(res3stage.ar2binding, xtable = TRUE)

## *** ar 1 non-binding
res3stage.ar1nonbinding <- res3stage.red[binding==FALSE & ar==1 & !is.na(decision),.SD,.SDcols=keep.col]

## number of simulations
res3stage.ar1nonbinding[type!="interim",.N,by=c("method.char","hypo")][,unique(N)]
## [1] 10000

table3stage.ar1nonbinding <- createTableResSim(res3stage.ar1nonbinding, xtable = FALSE)
table3stage.ar1nonbinding
##            statistic      method 1  method 1 fixC      method 2  method 2 fixC      method 3
##  1:     type 1 error         2.26%          2.21%         2.26%          2.20%         2.27%
##  2:            power        91.08%         90.51%        91.02%         90.67%        90.58%
##  3:       CI [NA,NA]         4.09%          4.66%         4.19%          4.45%         4.79%
##  4:         coverage        96.33%         98.09%        96.35%         98.14%        97.15%
##  5:         reversal   0.20%/0.16%    0.09%/0.62%   0.20%/0.16%    0.08%/0.61%       0/0.99%
##  6:         abnormal       0.57%/0            0/0       0.56%/0            0/0       0/0.07%
##  7:   mean bias LMME 0.007 / 0.096  0.007 / 0.098 0.007 / 0.096  0.007 / 0.097 0.007 / 0.099
##  8:    mean bias MUE 0.007 / 0.069  0.006 / 0.023 0.007 / 0.069  0.006 / 0.020 0.008 / 0.052
##  9: median bias LMME 0.64% / 8.51%  0.61% / 8.86% 0.64% / 8.57%  0.60% / 8.69% 0.56% / 8.98%
## 10:  median bias MUE 0.61% / 2.52% 0.58% / -1.68% 0.62% / 2.57% 0.60% / -2.18% 0.68% / 0.42%
createTableResSim(res3stage.ar1nonbinding, xtable = TRUE)

## *** ar 2 non-binding
res3stage.ar2nonbinding <- res3stage.red[binding==FALSE & ar==2 & !is.na(decision),.SD,.SDcols=keep.col]

## number of simulations
res3stage.ar2nonbinding[type!="interim",.N,by=c("method.char","hypo")][,unique(N)]
## [2] 10000

table3stage.ar2nonbinding <- createTableResSim(res3stage.ar2nonbinding, xtable = FALSE)
table3stage.ar2nonbinding
##            statistic      method 1 method 1 fixC      method 2 method 2 fixC      method 3
##  1:     type 1 error         2.17%         2.08%         2.11%         2.07%         2.33%
##  2:            power        91.43%        90.64%        91.35%        90.78%        90.69%
##  3:       CI [NA,NA]         3.69%         4.48%         3.95%         4.26%         4.72%
##  4:         coverage        96.67%        98.08%        96.73%        98.11%        97.61%
##  5:         reversal   0.85%/0.29%   0.60%/0.83%   1.02%/0.34%   0.67%/0.83%       0/1.25%
##  6:         abnormal       0.79%/0           0/0       0.76%/0           0/0       0/0.41%
##  7:   mean bias LMME 0.005 / 0.064 0.004 / 0.067 0.005 / 0.065 0.004 / 0.066 0.004 / 0.067
##  8:    mean bias MUE 0.005 / 0.051 0.004 / 0.021 0.004 / 0.053 0.004 / 0.019 0.007 / 0.035
##  9: median bias LMME 0.63% / 5.72% 0.59% / 6.18% 0.61% / 5.87% 0.57% / 6.07% 0.46% / 6.15%
## 10:  median bias MUE 0.61% / 2.39% 0.57% / 0.19% 0.59% / 2.54% 0.54% / 0.03% 0.81% / 0.50%
createTableResSim(res3stage.ar2nonbinding, xtable = TRUE)


## * text article

## ** Boundaries below 1.96
table2stage.Ck <- res2stage[decision %in% c("efficacy","futility"),
                        .(.N, rejection = mean2pc(decision=="efficacy"), rejectionBelow196 = mean2pc((statistic<qnorm(0.975))*(decision=="efficacy"))), by = c("scenario","missing","method","binding","fixC","ar","hypo")]
table2stage.Ck[method %in% 1:2 & fixC == FALSE & ar==10 & hypo == "power",
               range(rejectionBelow196)]
## [1] "0.64%" "0.85%"

table3stage.Ck <- res3stage[decision %in% c("efficacy","futility"),
                        .(.N, rejection = mean2pc(decision=="efficacy"), rejectionBelow196 = mean2pc((statistic<qnorm(0.975))*(decision=="efficacy"))), by = c("scenario","missing","method","binding","fixC","ar","hypo")]
table3stage.Ck[method %in% 1:2 & fixC == FALSE & ar==10 & hypo == "power",
               range(rejectionBelow196)]
## [1] "0.92%" "1.03%"

## ** Special cases
res2stage[, reasonNA := ifelse(is.na(reason),"NA",reason)]
res3stage[, reasonNA := ifelse(is.na(reason),"NA",reason)]
normal.case <- c("efficacy","futility","no boundary crossed","NA")

table2stage.special <- ftable(reason = res2stage[reasonNA %in% normal.case == FALSE,reason],
                              method = res2stage[reasonNA %in% normal.case == FALSE,method],
                              scenario = res2stage[reasonNA %in% normal.case == FALSE,scenario])
table2stage.special
##                                     scenario  1  2  3  4  5  6  7  8  9 13 17
## reason                       method                                          
## decreasing information       1                0  0  1  1  0  0  1  1  0  0  0
##                              2                0  0  1  1  0  0  1  1  0  0  0
##                              3                0  0  1  1  0  0  1  1  0  0  0
## Imax reached                 1                1  1  0  0  1  1  0  0  0  0  0
##                              2                1  1  0  0  1  1  0  0  0  0  0
##                              3                1  1  0  0  1  1  0  0  0  0  0
## stop for futility at interim 1                0  0  0  0  0  0  0  0  0  0  0
##                              2                0  0  0  0  0  0  0  0  0  0  0
##                              3               11  1  2  0 11  1  2  0  8  8  1

table3stage.special <- ftable(reason = res3stage[reasonNA %in% normal.case == FALSE,reason],
                              method = res3stage[reasonNA %in% normal.case == FALSE,method],
                              scenario = res3stage[reasonNA %in% normal.case == FALSE,scenario])
table3stage.special
##                                     scenario  1  2  3  4  5  6  7  8  9 10 11 13 14 15 17
## reason                       method                                                      
## decreasing information       1                0  0  0  1  0  0  0  1  0  0  1  0  0  1  0
##                              2                0  0  0  1  0  0  0  1  0  0  1  0  0  1  0
##                              3                0  0  0  1  0  0  0  1  0  0  1  0  0  1  0
## Imax reached                 1               28 13  0  0 28 13  0  0 24 43  0 24 43  0  0
##                              2               31 13  0  0 23 10  0  0 18 29  0 25 44  0  0
##                              3               25 15  0  0 25 15  0  0 24 44  0 24 44  0  0
## stop for futility at interim 1                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##                              2                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##                              3               28  2  1  0 28  2  1  0 26  0  1 26  0  1  3

## ** Type 1 error
## *** 2 stages
## For each run, create a binary indicator for rejection for efficacy
res2stage.rejection <- res2stage[,.(n.stage = .N, rejection = "efficacy" %in% na.omit(decision)),
                                 by = c("method.char","seed","scenario","missing","binding","fixC","ar","hypo")]

## Average over runs and method within scenario
res2stageS.rejection <- res2stage.rejection[,.(n.sim = .N, rejectionRate = 100*mean(rejection)),
                                            by=c("method.char","scenario","binding","missing","fixC","ar","hypo")]


res2stageS.rejection[missing & ar == 10 & hypo == "typeI", range(rejectionRate)]
## [1] 2.22 2.57
res2stageS.rejection[missing & ar == 10 & hypo == "power", range(rejectionRate)]
## [1] 79.86 81.00

res2stageS.rejection[missing & ar == 10 & hypo == "power" & method.char == "method 1",
                     rejectionRate]- res2stageS.rejection[missing & ar == 10 & hypo == "power" & method.char == "method 2",
                     rejectionRate]
## [1] 0.07 0.06

res2stageS.rejection[missing & ar == 10 & hypo == "power" & method.char == "method 1",
                     rejectionRate]- res2stageS.rejection[missing & ar == 10 & hypo == "power" & method.char == "method 3",
                     rejectionRate]
## [1] 0.57 0.24

## *** 3 stages
## For each run, create a binary indicator for rejection for efficacy
res3stage.rejection <- res3stage[,.(n.stage = .N, rejection = "efficacy" %in% na.omit(decision)),
                                 by = c("method.char","seed","scenario","missing","binding","fixC","ar","hypo")]

## Average over runs and method within scenario
res3stageS.rejection <- res3stage.rejection[,.(n.sim = .N, rejectionRate = 100*mean(rejection)),
                                            by=c("method.char","scenario","binding","missing","fixC","ar","hypo")]


res3stageS.rejection[missing & hypo == "typeI", range(rejectionRate)]
## [1] 2.31 2.61
res3stageS.rejection[missing & hypo == "power", range(rejectionRate)]
## [1] 79.58 80.87

res3stageS.rejection[missing & hypo == "power" & method.char == "method 1",
                     rejectionRate]- res3stageS.rejection[missing & hypo == "power" & method.char == "method 2",
                     rejectionRate]
## [1] 0.08 0.01 0.08 0.02

res3stageS.rejection[missing & hypo == "power" & method.char == "method 1",
                     rejectionRate]- res3stageS.rejection[missing & hypo == "power" & method.char == "method 3",
                     rejectionRate]
## [1] 0.73 0.47 0.46 0.48

## ** reversal
## *** 2 stages
res2stageS.reversal <- res2stage.conclusion[missing & hypo == "power" & ar == 10, .(N = .N,
                                                fu2eff = 100*mean(stage.interim == stage.conclusion & decision.interim == "stop" & reason.interim == "futility" & conclusion == "efficacy"),
                                                eff2fu = 100*mean(stage.interim == stage.conclusion & decision.interim == "stop" & reason.interim == "efficacy" & conclusion == "futility")),
                                            by = c("method","scenario","missing","binding","fixC","ar","hypo")]
res2stageS.reversalW <- dcast(res2stageS.reversal, scenario + N + hypo + missing + ar + binding + fixC ~ method, value.var = c("fu2eff","eff2fu"))

range(c(res2stageS.reversalW$fu2eff_1,res2stageS.reversalW$fu2eff_2))
## [1] 0.11 0.61
range(c(res2stageS.reversalW$eff2fu_1,res2stageS.reversalW$eff2fu_2))
## [1] 0.17 0.67
range(res2stageS.reversalW$eff2fu_3)
## [1] 1.04 1.07

## *** 3 stages
res3stageS.reversal <- res3stage.conclusion[missing == TRUE, .(N = .N,
                                                fu2eff = 100*mean(stage.interim == stage.conclusion & decision.interim == "stop" & reason.interim == "futility" & conclusion == "efficacy"),
                                                eff2fu = 100*mean(stage.interim == stage.conclusion & decision.interim == "stop" & reason.interim == "efficacy" & conclusion == "futility")),
                                            by = c("method","scenario","missing","binding","fixC","ar","hypo")]
res3stageS.reversalW <- dcast(res3stageS.reversal, N + scenario + hypo + missing + ar + binding + fixC ~ method, value.var = c("fu2eff","eff2fu"))

range(c(res3stageS.reversalW$fu2eff_1,res3stageS.reversalW$fu2eff_2))
## [1] 0.00 0.69
range(c(res3stageS.reversalW$eff2fu_1,res3stageS.reversalW$eff2fu_2))
## [1] 0.01 1.10
range(res3stageS.reversalW$eff2fu_3)
## [1] 0.12 1.60

## ** coverage
res2stage.coverage <- res2stage[missing==TRUE & hypo=="power" & ar == 10 & decision %in% c("futility","efficacy"),
                                .(N = .N,
                                  "NA" = 100*mean(is.na(lower_MUE) | is.na(upper_MUE)),
                                  coverage = 100*mean( (lower_MUE <= truth) & (truth <= upper_MUE), na.rm=TRUE)),
                                by = c("method.char","missing","binding","fixC")]

range(res2stage.coverage$"NA",na.rm=TRUE)
## [1] 0.00 6.57
range(res2stage.coverage$coverage,na.rm=TRUE)
## [1] 94.83 97.47

res3stage.coverage <- res3stage[missing==TRUE & hypo=="power" & ar == 10 & decision %in% c("futility","efficacy"),
                                .(N = .N,
                                  "NA" = 100*mean(is.na(lower_MUE) | is.na(upper_MUE)),
                                  coverage = 100*mean( (lower_MUE <= truth) & (truth <= upper_MUE), na.rm=TRUE)),
                                by = c("method.char","missing","binding","fixC")]

range(res3stage.coverage$"NA",na.rm=TRUE)
## [1] 0.01 9.05
range(res3stage.coverage$coverage,na.rm=TRUE)
## [1] 95.29 98.09


## ** coherence
res2stage.PmismatchFU <- res2stage[decision=="futility",.(N = .N,
                                                          mismatchP = 100*mean(p.value_MUE<0.025),
                                                          mismatchCI = 100*mean(lower_MUE>0, na.rm=TRUE)),
                                  by = c("method.char","scenario","missing","binding","fixC","ar","hypo")]
res2stage.PmismatchFU[,.(max(mismatchP),max(mismatchCI))]
## 1:  0  0

res2stage.PmismatchEFF <- res2stage[decision=="efficacy",.(N = .N,
                                                           mismatchP = 100*mean(p.value_MUE>0.025),
                                                           mismatchCI = 100*mean(lower_MUE<0, na.rm=TRUE)),
                                  by = c("method.char","scenario","missing","binding","fixC","ar","hypo")]
res2stage.PmismatchEFF[,.(max(mismatchP),max(mismatchCI))]
## 1:  0  0

res3stage.PmismatchFU <- res3stage[decision=="futility",.(N = .N,
                                                          mismatchP = 100*mean(p.value_MUE<0.025),
                                                          mismatchCI = 100*mean(lower_MUE>0, na.rm=TRUE)),
                                  by = c("method.char","scenario","missing","binding","fixC","ar","hypo")]
res3stage.PmismatchFU[,.(max(mismatchP),max(mismatchCI))]
## 1: 0.04985  0
print(res3stage[decision=="futility" & p.value_MUE<0.025,min(p.value_MUE)], digits = 10)
## 0.02500576079

res3stage.PmismatchEFF <- res3stage[decision=="efficacy",.(N = .N,
                                                           mismatchP = 100*mean(p.value_MUE>0.025),
                                                           mismatchCI = 100*mean(lower_MUE<0, na.rm=TRUE)),
                                  by = c("method.char","scenario","missing","binding","fixC","ar","hypo")]
res3stage.PmismatchEFF[,.(max(mismatchP),max(mismatchCI))]
## 1: 0.4115 0.01249
print(res3stage[decision=="efficacy" & p.value_MUE>0.025,max(p.value_MUE)], digits = 10)
## [1] 0.02500576079
print(res3stage[decision=="efficacy",min(lower_MUE)], digits = 10)

## ** bias

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

res2stageS.bias[missing == TRUE & binding == TRUE & fixC == FALSE & ar == 10,.(bias_MLE,bias_MUE,ratio=abs(bias_MLE)/abs(bias_MUE)), by = c("method","scenario")]
##    method scenario bias_MLE  bias_MUE ratio
## 1:      1        1  0.01345  0.005985 2.247
## 2:      2        1  0.01315  0.005661 2.323
## 3:      3        1  0.01468 -0.003299 4.450
## 4:      1        2 -0.01794 -0.004533 3.957
## 5:      2        2 -0.01784 -0.004483 3.981
## 6:      3        2 -0.01856 -0.016750 1.108
res2stageS.bias[missing == TRUE & binding == TRUE & fixC == FALSE & ar == 10,.(mbias_MLE,mbias_MUE,ratio=abs(mbias_MLE)/abs(mbias_MUE)), by = c("method","scenario")]
##    method scenario mbias_MLE  mbias_MUE   ratio
## 1:      1        1    0.0261 -0.0024000 10.8750
## 2:      2        1    0.0260 -0.0025000 10.4000
## 3:      3        1    0.0301 -0.0054505  5.5224
## 4:      1        2   -0.0173  0.0010002 17.2965
## 5:      2        2   -0.0170  0.0007502 22.6599
## 6:      3        2   -0.0202 -0.0232523  0.8687

res2stageS.bias[missing == TRUE & (binding == FALSE | fixC == TRUE) & ar == 10,.(mbias_MLE,mbias_MUE,ratio=abs(mbias_MLE)/abs(mbias_MUE)), by = c("method","scenario")]


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

res3stageS.bias[missing == TRUE & binding == TRUE & fixC == FALSE & ar == 10,.(bias_MLE,bias_MUE,ratio=abs(bias_MLE)/abs(bias_MUE)), by = c("method","scenario")]
##    method scenario bias_MLE  bias_MUE ratio
## 1:      1        1  0.02280  0.016230 1.405
## 2:      2        1  0.02262  0.016051 1.410
## 3:      3        1  0.02481  0.005772 4.299
## 4:      1        2 -0.03396 -0.014852 2.286
## 5:      2        2 -0.03380 -0.014699 2.300
## 6:      3        2 -0.03400 -0.028405 1.197
res3stageS.bias[missing == TRUE & binding == TRUE & fixC == FALSE & ar == 10,.(mbias_MLE,mbias_MUE,ratio=abs(mbias_MLE)/abs(mbias_MUE)), by = c("method","scenario")]
##    method scenario mbias_MLE mbias_MUE  ratio
## 1:      1        1    0.0359 -0.003801  9.445
## 2:      2        1    0.0358 -0.004001  8.948
## 3:      3        1    0.0374 -0.009511  3.932
## 4:      1        2   -0.0367  0.010353  3.545
## 5:      2        2   -0.0365  0.010204  3.577
## 6:      3        2   -0.0380 -0.002351 16.162
##----------------------------------------------------------------------
### table1.R ends here
