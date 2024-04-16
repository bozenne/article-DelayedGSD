### table1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 19 2023 (10:24) 
## Version: 
## Last-Updated: nov  9 2023 (16:24) 
##           By: Brice Ozenne
##     Update #: 60
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
  mean.x <- mean(x, na.rm = TRUE)
  out <- paste0(formatC(100*mean.x, format='f', digits = digits), "%")
  if(round(100*mean.x,digits)==0){
    out[round(100*mean.x,digits)==0] <- paste0("<0.",rep(0,digits-1),"1%")
  }
  if(abs(mean.x)<1e-12){
    out[abs(mean.x)<1e-12] <- "0"
  }
  if(any(is.na(x))){
    out <- paste0(out," (NA: ",dec2pc(mean(is.na(x)), digits = digits),")")
  }
  return(out)
}

##' @description Average and nicely format as numeric (fixed number of digits, e.g. 0.010%)
mean2num <- function(x, digits = 3){
  out <- formatC(mean(x, na.rm = TRUE), format='f', digits = digits)
  if(any(is.na(x))){
    out <- paste0(out," (NA: ",dec2pc(mean(is.na(x)), digits = digits),")")
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
                                                      value = mean2pc(decision=="efficacy")),
                                      by="method.char"]
    ls.sumstat$power <- data.decision[hypo=="power",.(statistic = "power",
                                                      value = mean2pc(decision=="efficacy")),
                                      by="method.char"]
    ls.sumstat$CINA <- data.decision[hypo=="power",.(statistic = "CI [NA,NA]",
                                                     value = mean2pc(is.na(lower_MUE) | is.na(upper_MUE))),
                                     by="method.char"]
    ls.sumstat$coverage <- data.decision[hypo=="power" & !is.na(lower_MUE) & !is.na(upper_MUE),.(statistic = "coverage",
                                                                                                 value = mean2pc(lower_MUE <= truth & truth <= upper_MUE)),
                                         by="method.char"]
    ls.sumstat$reversal <- data.decision[hypo=="power",.(statistic = "reversal",
                                                         value = paste(mean2pc(decision.interim=="stop (futility)" & decision=="efficacy"),
                                                                       mean2pc(decision.interim=="stop (efficacy)" & decision=="futility"),
                                                                       sep="/")),
                                         by="method.char"]
    ls.sumstat$abnormal <- data.decision[hypo=="power",.(statistic = "abnormal",
                                                         value = paste(mean2pc(decision=="efficacy" & statistic < 1.96),
                                                                       mean2pc(decision=="futility" & statistic >= ck),
                                                                       sep="/")),
                                         by="method.char"]
    ls.sumstat$LMMEmeanB <- data.decision[!is.na(estimate_MUE),.(statistic = "mean bias LMME",
                                                                 value = paste(mean2num(.SD[hypo=="typeI",estimate_ML]-.SD[hypo=="typeI",truth]),
                                                                               mean2num(.SD[hypo=="power",estimate_ML]-.SD[hypo=="power",truth]),
                                                                               sep=sep)),
                                          by="method.char"]
    ls.sumstat$MUEmeanB <- data.decision[!is.na(estimate_MUE),.(statistic = "mean bias MUE",
                                                                value = paste(mean2num(.SD[hypo=="typeI",estimate_MUE]-.SD[hypo=="typeI",truth]),
                                                                              mean2num(.SD[hypo=="power",estimate_MUE]-.SD[hypo=="power",truth]),
                                                                              sep=sep)),
                                         by="method.char"]
    ls.sumstat$LMMEmedianB <- data.decision[!is.na(estimate_MUE),.(statistic = "median bias LMME",
                                                                   value = paste(mean2pc((.SD[hypo=="typeI",estimate_ML]>.SD[hypo=="typeI",truth])-0.5),
                                                                                 mean2pc((.SD[hypo=="power",estimate_ML]>.SD[hypo=="power",truth])-0.5),
                                                                                 sep=sep)),
                                            by="method.char"]
    ls.sumstat$MUEmedianB <- data.decision[!is.na(estimate_MUE),.(statistic = "median bias MUE",
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
res2stage[, truth := ifelse(hypo=="power",0.6,0)]

## ** Generate table
keep.col <- c("scenario", "hypo", "method", "stage", "type", "statistic", "ck",
              "estimate_ML", "se_ML", "p.value_ML", "lower_ML", "upper_ML",
              "estimate_MUE", "p.value_MUE", "lower_MUE", "upper_MUE",
              "decision", "reason", "method.char", "stage.char", "truth","seed")       

## *** ar 5 binding
res2stage.ar5binding <- res2stage[method.char != "method 3 fixC" & missing==TRUE & binding==TRUE & ar==5 & !is.na(decision),.SD,.SDcols=keep.col]
res2stage.ar5binding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000
table2stage.ar5binding <- createTableResSim(res2stage.ar5binding, xtable = FALSE)
table2stage.ar5binding
##           statistic     method 1 method 1 fixC     method 2 method 2 fixC     method 3
##  1:     type 1 error        2.40%         2.32%        2.40%         2.31%        2.35%
##  2:            power       80.53%        80.08%       80.53%        80.20%       80.14%
##  3:       CI [NA,NA]        0.01%         0.01%        0.01%         0.01%        0.01%
##  4:         coverage       94.74%        96.29%       94.74%        96.33%       95.14%
##  5:         reversal  0.08%/0.07%   0.02%/0.46%  0.08%/0.07%   0.02%/0.45%      0/0.67%
##  6:         abnormal      0.45%/0           0/0      0.45%/0           0/0      0/0.02%
##  7:   mean bias LMME -0.030/0.023  -0.030/0.023 -0.030/0.023  -0.031/0.023 -0.031/0.024
##  8:    mean bias MUE -0.012/0.010 -0.012/-0.015 -0.012/0.010 -0.013/-0.015 -0.012/0.004
##  9: median bias LMME -3.29%/4.05%  -3.29%/4.05% -3.28%/4.05%  -3.45%/4.07% -3.44%/4.32%
## 10:  median bias MUE 0.06%/-0.35%  0.06%/-0.87% 0.06%/-0.34%  0.08%/-0.76% 0.07%/-0.55%
createTableResSim(res2stage.ar5binding, xtable = TRUE)

res2stage.ar5binding[type %in% c("decision","final") & hypo=="power" & method == 1 & method.char == "method 1", mean(estimate_MUE<truth, na.rm = TRUE)]
res2stage.ar5binding[type %in% c("decision","final") & hypo=="power" & method == 1 & method.char == "method 1 fixC", mean(estimate_MUE<truth, na.rm = TRUE)]
res2stage.ar5binding[type %in% c("decision","final") & hypo=="typeI" & method == 1 & method.char == "method 1", mean(estimate_MUE<truth, na.rm = TRUE)]
res2stage.ar5binding[type %in% c("decision","final") & hypo=="typeI" & method == 1 & method.char == "method 1 fixC", mean(estimate_MUE<truth, na.rm = TRUE)]

## [1] 0.7006
debug2stage <- res2stage.ar5binding[type %in% c("decision","final") & hypo=="typeI" & method == 1]

debug2stage[method.char == "method 1", mean(ck<=1.95)]
debug2stage[method.char == "method 1 fixC", mean(ck<=1.95)]

debug2stage[method.char == "method 1", estimate_MUE] - debug2stage[method.char == "method 1 fixC", estimate_MUE]
debug2stage[method.char == "method 1", estimate_MUE] - debug2stage[method.char == "method 1 fixC", estimate_MUE]

debug2stage[method.char == "method 1"]

res2stage.ar5binding[type %in% c("decision","final") & hypo=="typeI" & method == 1 & method.char == "method 1"]
res2stage.ar5binding[type %in% c("decision","final") & hypo=="typeI" & method == 1 & method.char == "method 1 fixC", hist(ck)]


## res2stage.ar5binding[type %in% c("decision","final") & hypo=="typeI" & method == 1 & seed == 402321297]
## res2stage[type %in% "decision" & hypo=="typeI" & method == 1 & seed == 402321297 & ar == 5 & binding == TRUE & missing == TRUE]

power2stage.ar5binding <- as.numeric(gsub("%","",unlist(table2stage.ar5binding[2,.SD,.SDcols = names(table2stage.ar5binding)[-1]]),fixed=TRUE))
power2stage.ar5binding[1:2]-power2stage.ar5binding[3:4]
power2stage.ar5binding[1:2]-power2stage.ar5binding[5]
## > [1]  0.00 -0.12
## > [1]  0.39 -0.06

## *** ar 10 binding
res2stage.ar10binding <- res2stage[method.char != "method 3 fixC" & missing==TRUE & binding==TRUE & ar==10 & !is.na(decision),.SD,.SDcols=keep.col]
res2stage.ar10binding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]

res2stage.ar10binding[hypo=="power",.(statistic = "reversal",
                                      value = paste(mean2pc(decision.interim=="stop (futility)" & decision=="efficacy"),
                                                    mean2pc(decision.interim=="stop (efficacy)" & decision=="futility"),
                                                    sep="/")),
                      by="method.char"]

## [1] 10000
table2stage.ar10binding <- createTableResSim(res2stage.ar10binding, xtable = FALSE)
table2stage.ar10binding
##           statistic     method 1 method 1 fixC     method 2 method 2 fixC      method 3
##  1:     type 1 error        2.42%         2.24%        2.39%         2.22%         2.37%
##  2:            power       81.00%        80.15%       80.93%        80.35%        80.43%
##  3:       CI [NA,NA]        0.02%         0.02%        0.02%         0.02%         0.02%
##  4:         coverage       94.84%        96.26%       94.82%        96.31%        95.34%
##  5:         reversal  0.57%/0.17%   0.22%/0.67%  0.61%/0.20%   0.16%/0.65%       0/1.07%
##  6:         abnormal      0.85%/0           0/0      0.84%/0           0/0       0/0.11%
##  7:   mean bias LMME -0.018/0.013  -0.018/0.013 -0.018/0.013  -0.019/0.014  -0.019/0.015
##  8:    mean bias MUE -0.005/0.006 -0.006/-0.015 -0.004/0.006 -0.006/-0.015 -0.005/-0.003
##  9: median bias LMME -1.72%/2.61%  -1.72%/2.61% -1.69%/2.60%  -1.97%/2.65%  -2.02%/3.01%
## 10:  median bias MUE 0.10%/-0.24%  0.10%/-1.11% 0.08%/-0.25% -0.07%/-1.05% -0.02%/-0.55%
createTableResSim(res2stage.ar10binding, xtable = TRUE)

power2stage.ar10binding <- as.numeric(gsub("%","",unlist(table2stage.ar10binding[2,.SD,.SDcols = names(table2stage.ar10binding)[-1]]),fixed=TRUE))
power2stage.ar10binding[1:2]-power2stage.ar10binding[3:4]
power2stage.ar10binding[1:2]-power2stage.ar10binding[5]
## [1]  0.07 -0.20
## [1]  0.57 -0.28

## *** ar 5 non-binding
res2stage.ar5nonbinding <- res2stage[method.char != "method 3 fixC" & missing==TRUE & binding==FALSE & ar==5 & !is.na(decision),.SD,.SDcols=keep.col]
res2stage.ar5nonbinding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000
table2stage.ar5nonbinding <- createTableResSim(res2stage.ar5nonbinding, xtable = FALSE)
table2stage.ar5nonbinding
##            statistic     method 1 method 1 fixC     method 2 method 2 fixC     method 3
##  1:     type 1 error        2.68%         2.63%        2.68%         2.64%        2.66%
##  2:            power       80.37%        79.93%       80.36%        80.04%       80.06%
##  3:       CI [NA,NA]        5.72%         6.16%        5.76%         5.86%        6.05%
##  4:         coverage       95.86%        97.77%       95.86%        97.76%       96.55%
##  5:         reversal  0.03%/0.04%   0.01%/0.46%  0.03%/0.04%   0.01%/0.44%      0/0.60%
##  6:         abnormal      0.44%/0           0/0      0.44%/0           0/0          0/0
##  7:   mean bias LMME  0.001/0.054   0.001/0.055  0.001/0.054   0.001/0.054  0.001/0.055
##  8:    mean bias MUE  0.001/0.042   0.000/0.017  0.001/0.042   0.000/0.015  0.001/0.037
##  9: median bias LMME -0.10%/7.92%  -0.12%/8.19% -0.10%/7.97%  -0.13%/7.97% -0.14%/8.40%
## 10:  median bias MUE -0.11%/2.83%  -0.13%/2.57% -0.11%/2.87%  -0.14%/2.47% -0.02%/2.80%
createTableResSim(res2stage.ar5nonbinding, xtable = TRUE)

power2stage.ar5nonbinding <- as.numeric(gsub("%","",unlist(table2stage.ar5nonbinding[2,.SD,.SDcols = names(table2stage.ar5nonbinding)[-1]]),fixed=TRUE))
power2stage.ar5nonbinding[1:2]-power2stage.ar5nonbinding[3:4]
power2stage.ar5nonbinding[1:2]-power2stage.ar5nonbinding[5]
## [1]  0.01 -0.11
## [1]  0.31 -0.13

## *** ar 10 non-binding
res2stage.ar10nonbinding <- res2stage[method.char != "method 3 fixC" & missing==TRUE & binding==FALSE & ar==10 & !is.na(decision),.SD,.SDcols=keep.col]
res2stage.ar10nonbinding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000
table2stage.ar10nonbinding <- createTableResSim(res2stage.ar10nonbinding, xtable = FALSE)
table2stage.ar10nonbinding
##            statistic     method 1 method 1 fixC     method 2 method 2 fixC     method 3
##  1:     type 1 error        2.53%         2.45%        2.53%         2.47%        2.57%
##  2:            power       80.50%        79.86%       80.44%        80.12%       80.26%
##  3:       CI [NA,NA]        5.95%         6.59%        6.11%         6.06%        6.21%
##  4:         coverage       95.90%        97.38%       95.89%        97.45%       97.00%
##  5:         reversal  0.41%/0.21%   0.14%/0.58%  0.42%/0.22%   0.11%/0.55%      0/1.04%
##  6:         abnormal      0.64%/0           0/0      0.64%/0           0/0      0/0.08%
##  7:   mean bias LMME -0.000/0.042  -0.001/0.043 -0.000/0.042  -0.001/0.042 -0.001/0.042
##  8:    mean bias MUE -0.000/0.036  -0.001/0.017 -0.000/0.036  -0.001/0.015  0.001/0.029
##  9: median bias LMME -0.15%/6.61%  -0.19%/7.00% -0.15%/6.69%  -0.19%/6.74% -0.26%/6.79%
## 10:  median bias MUE -0.15%/3.09%  -0.19%/2.72% -0.14%/3.18%  -0.19%/2.47% -0.03%/2.80%
createTableResSim(res2stage.ar10nonbinding, xtable = TRUE)

power2stage.ar10nonbinding <- as.numeric(gsub("%","",unlist(table2stage.ar10nonbinding[2,.SD,.SDcols = names(table2stage.ar10nonbinding)[-1]]),fixed=TRUE))
power2stage.ar10nonbinding[1:2]-power2stage.ar10nonbinding[3:4]
power2stage.ar10nonbinding[1:2]-power2stage.ar10nonbinding[5]
## [1]  0.06 -0.26
## [1]  0.24 -0.40

## * 3 stages

## ** Load data
res3stage <- readRDS(file.path("Results-built","res3stage.rds"))
res3stage[, method.char := paste0("method ",method, c(""," fixC")[fixC+1])]
res3stage[, stage.char := factor(stage, 1:3, c("interim1","interim2","final"))]
res3stage[, truth := ifelse(hypo=="power",0.6,0)]

## ** Generate table
keep.col <- c("scenario", "hypo", "method", "stage", "type", "statistic", "ck",
              "estimate_ML", "se_ML", "p.value_ML", "lower_ML", "upper_ML",
              "estimate_MUE", "p.value_MUE", "lower_MUE", "upper_MUE",
              "decision", "reason", "method.char", "stage.char", "truth","seed")       

## *** ar 5 binding
res3stage.ar5binding <- res3stage[method.char != "method 3 fixC" & missing==TRUE & binding==TRUE & ar==5 & !is.na(decision),.SD,.SDcols=keep.col]
res3stage.ar5binding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000
table3stage.ar5binding <- createTableResSim(res3stage.ar5binding, xtable = FALSE)
table3stage.ar5binding
##            statistic       method 1  method 1 fixC       method 2  method 2 fixC       method 3
##  1:     type 1 error          2.61%          2.56%          2.61%          2.55%          2.59%
##  2:            power         74.35%         74.00%         74.35%         74.01%         73.99%
##  3:       CI [NA,NA]              0          0.08%              0          0.08%          0.11%
##  4:         coverage         94.63%         96.84%         94.63%         96.82%         95.68%
##  5:         reversal    0.01%/0.02%        0/0.36%    0.01%/0.02%        0/0.35%        0/0.45%
##  6:         abnormal        0.35%/0            0/0        0.35%/0            0/0            0/0
##  7:   mean bias LMME -0.056 / 0.034 -0.056 / 0.034 -0.056 / 0.034 -0.056 / 0.035 -0.056 / 0.035
##  8:    mean bias MUE -0.034 / 0.026 -0.034 / 0.005 -0.034 / 0.026 -0.034 / 0.006 -0.034 / 0.016
##  9: median bias LMME -6.59% / 3.91% -6.60% / 3.89% -6.57% / 3.90% -6.59% / 3.85% -6.60% / 3.97%
## 10:  median bias MUE -0.15% / 0.81% -0.16% / 0.50% -0.12% / 0.80% -0.12% / 0.53% -0.13% / 0.53%
createTableResSim(res3stage.ar5binding, xtable = TRUE)

## *** ar 10 binding
res3stage.ar10binding <- res3stage[method.char != "method 3 fixC" & missing==TRUE & binding==TRUE & ar==10 & !is.na(decision),.SD,.SDcols=keep.col]
res3stage.ar10binding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000
table3stage.ar10binding <- createTableResSim(res3stage.ar10binding, xtable = FALSE)
table3stage.ar10binding
##            statistic       method 1  method 1 fixC       method 2  method 2 fixC       method 3
##  1:     type 1 error          2.60%          2.47%          2.60%          2.50%          2.49%
##  2:            power         74.51%         73.84%         74.51%         73.91%         74.01%
##  3:       CI [NA,NA]          0.12%          0.20%          0.12%          0.20%          0.23%
##  4:         coverage         95.07%         96.87%         95.07%         96.91%         96.12%
##  5:         reversal    0.25%/0.08%    0.04%/0.54%    0.25%/0.08%    0.04%/0.50%        0/0.73%
##  6:         abnormal        0.67%/0            0/0        0.67%/0            0/0        0/0.03%
##  7:   mean bias LMME -0.035 / 0.021 -0.035 / 0.021 -0.035 / 0.021 -0.034 / 0.022 -0.034 / 0.023
##  8:    mean bias MUE -0.027 / 0.023 -0.028 / 0.008 -0.027 / 0.023 -0.028 / 0.008 -0.028 / 0.014
##  9: median bias LMME -3.89% / 2.84% -3.91% / 2.85% -3.88% / 2.84% -3.80% / 2.85% -3.80% / 3.07%
## 10:  median bias MUE -0.01% / 0.72% -0.05% / 0.13%  0.02% / 0.70%  0.08% / 0.19%  0.13% / 0.25%
createTableResSim(res3stage.ar10binding, xtable = TRUE)

## *** ar 5 non-binding
res3stage.ar5nonbinding <- res3stage[method.char != "method 3 fixC" & missing==TRUE & binding==FALSE & ar==5 & !is.na(decision),.SD,.SDcols=keep.col]
res3stage.ar5nonbinding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000
table3stage.ar5nonbinding <- createTableResSim(res3stage.ar5nonbinding, xtable = FALSE)
table3stage.ar5nonbinding
##            statistic       method 1  method 1 fixC       method 2  method 2 fixC       method 3
##  1:     type 1 error          2.53%          2.43%          2.53%          2.44%          2.44%
##  2:            power         74.61%         74.33%         74.59%         74.37%         74.43%
##  3:       CI [NA,NA]          4.85%          5.13%          4.86%          5.00%          5.09%
##  4:         coverage         95.87%         98.24%         95.89%         98.24%         97.00%
##  5:         reversal    0.01%/0.01%    0.01%/0.29%    0.01%/0.01%    0.01%/0.28%        0/0.32%
##  6:         abnormal        0.28%/0            0/0        0.28%/0            0/0        0/0.01%
##  7:   mean bias LMME  0.004 / 0.065  0.003 / 0.065  0.004 / 0.065  0.003 / 0.065  0.003 / 0.066
##  8:    mean bias MUE  0.004 / 0.057  0.003 / 0.038  0.004 / 0.057  0.003 / 0.037  0.004 / 0.049
##  9: median bias LMME -0.35% / 5.44% -0.40% / 5.54% -0.35% / 5.44% -0.40% / 5.43% -0.41% / 5.62%
## 10:  median bias MUE -0.34% / 2.55% -0.39% / 2.41% -0.35% / 2.55% -0.39% / 2.35% -0.39% / 2.52%
createTableResSim(res3stage.ar5nonbinding, xtable = TRUE)

## *** ar 10 non-binding
res3stage.ar10nonbinding <- res3stage[method.char != "method 3 fixC" & missing==TRUE & binding==FALSE & ar==10 & !is.na(decision),.SD,.SDcols=keep.col]
res3stage.ar10nonbinding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 10000
table3stage.ar10nonbinding <- createTableResSim(res3stage.ar10nonbinding, xtable = FALSE)
table3stage.ar10nonbinding
##            statistic       method 1  method 1 fixC       method 2  method 2 fixC       method 3
##  1:     type 1 error          2.49%          2.37%          2.49%          2.37%          2.42%
##  2:            power         74.71%         74.17%         74.71%         74.25%         74.45%
##  3:       CI [NA,NA]          4.78%          5.25%          4.79%          4.88%          4.99%
##  4:         coverage         96.02%         97.96%         96.01%         97.95%         97.41%
##  5:         reversal    0.21%/0.09%    0.05%/0.47%    0.21%/0.09%    0.04%/0.46%        0/0.57%
##  6:         abnormal        0.54%/0            0/0        0.54%/0            0/0        0/0.03%
##  7:   mean bias LMME  0.002 / 0.046  0.002 / 0.046  0.002 / 0.046  0.002 / 0.045  0.002 / 0.046
##  8:    mean bias MUE  0.003 / 0.051  0.002 / 0.037  0.003 / 0.051  0.002 / 0.036  0.003 / 0.042
##  9: median bias LMME -0.40% / 4.29% -0.45% / 4.60% -0.40% / 4.29% -0.45% / 4.37% -0.47% / 4.58%
## 10:  median bias MUE -0.40% / 2.60% -0.45% / 2.33% -0.40% / 2.59% -0.45% / 2.09% -0.39% / 2.21%
createTableResSim(res3stage.ar10nonbinding, xtable = TRUE)

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
