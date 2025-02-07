### table1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 19 2023 (10:24) 
## Version: 
## Last-Updated: feb  7 2025 (10:14) 
##           By: Brice Ozenne
##     Update #: 188
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
res2stage.red[, type := factor(type, c("interim","decision","final"))]
setkeyv(res2stage.red, c("scenario","method","type"))
## unique(res2stage.red$scenario)


## ** Generate table
keep.col <- c("scenario", "hypo", "method", "stage", "type", "statistic", "ck",
              "nX1", "nX2", "nX3", "infoPC",
              "estimate_ML", "se_ML", "p.value_ML", "lower_ML", "upper_ML",
              "estimate_MUE", "p.value_MUE", "lower_MUE", "upper_MUE",
              "decision", "reason", "method.char", "stage.char", "truth","seed")       

## *** ar 1 binding
res2stage.ar1binding <- res2stage.red[infoBias==0 & binding==TRUE & ar==1,.SD,.SDcols=keep.col]

## sample size
res2stage.ar1binding[, .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                         n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                         n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]"),
                         pc.info = mean(infoPC)), by = "type"]
##        type   n.patients    n.outcome n.missing pc.info
##      <fctr>       <char>       <char>    <char>   <num>
## 1:  interim 163[153;179] 119[109;129] 44[29;60]  0.5455
## 2: decision 170[158;189] 152[139;170]  18[8;27]  0.6779
## 3:    final 265[265;265] 237[237;237] 28[28;28]  1.0183

## number of simulations
res2stage.ar1binding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 20000

## table
table2stage.ar1binding <- createTableResSim(res2stage.ar1binding, xtable = FALSE)
table2stage.ar1binding
##            statistic        method 1   method 1 fixC        method 2   method 2 fixC        method 3
##               <fctr>          <char>          <char>          <char>          <char>          <char>
##  1:     type 1 error           2.61%           2.51%           2.62%           2.52%           2.56%
##  2:            power          90.45%          89.92%          90.43%          90.08%          89.95%
##  3:       CI [NA,NA]               0               0               0               0               0
##  4:         coverage          94.54%          95.44%          94.55%          95.57%          95.06%
##  5:         reversal     0.20%/0.08%     0.08%/0.50%     0.19%/0.08%     0.07%/0.48%         0/0.78%
##  6:         abnormal         0.53%/0             0/0         0.52%/0             0/0         0/0.04%
##  7:   mean bias LMME  -0.040 / 0.036  -0.040 / 0.036  -0.040 / 0.036  -0.040 / 0.037  -0.040 / 0.037
##  8:    mean bias MUE  -0.015 / 0.015 -0.097 / -0.035  -0.015 / 0.015 -0.100 / -0.036 -0.043 / -0.005
##  9: median bias LMME  -3.31% / 4.20%  -3.31% / 4.20%  -3.29% / 4.20%  -3.44% / 4.20%  -3.50% / 4.06%
## 10:  median bias MUE -0.34% / -0.90% -9.44% / -4.18% -0.34% / -0.90% -9.34% / -4.28% -3.42% / -1.93%
createTableResSim(res2stage.ar1binding, xtable = TRUE)

## *** ar 1 non-binding
res2stage.ar1nonbinding <- res2stage.red[infoBias==0 & binding==FALSE & ar==1,.SD,.SDcols=keep.col]

## sample size
res2stage.ar1nonbinding[, .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                            n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                            n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]"),
                         pc.info = mean(infoPC)), by = "type"]
##        type   n.patients    n.outcome n.missing pc.info
##      <fctr>       <char>       <char>    <char>   <num>
## 1:  interim 164[153;178] 120[111;130] 44[28;60]  0.5428
## 2: decision 171[159;187] 153[139;171]  18[8;26]  0.6779
## 3:    final 268[268;268] 240[240;240] 28[28;28]  1.0294

## number of simulations
res2stage.ar1nonbinding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 20000

## table
table2stage.ar1nonbinding <- createTableResSim(res2stage.ar1nonbinding, xtable = FALSE)
table2stage.ar1nonbinding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC        method 3
##               <fctr>         <char>          <char>         <char>          <char>          <char>
##  1:     type 1 error          2.62%           2.55%          2.62%           2.56%           2.62%
##  2:            power         90.72%          90.31%         90.70%          90.48%          90.44%
##  3:       CI [NA,NA]          3.12%           3.52%          3.16%           3.31%           3.50%
##  4:         coverage         95.56%          97.00%         95.57%          97.06%          96.36%
##  5:         reversal    0.17%/0.08%     0.05%/0.36%    0.17%/0.07%     0.04%/0.35%         0/0.63%
##  6:         abnormal        0.40%/0             0/0        0.40%/0             0/0         0/0.02%
##  7:   mean bias LMME  0.002 / 0.062   0.002 / 0.064  0.002 / 0.062   0.002 / 0.063   0.002 / 0.064
##  8:    mean bias MUE  0.003 / 0.041  0.001 / -0.004  0.002 / 0.041  0.001 / -0.007   0.003 / 0.025
##  9: median bias LMME -0.38% / 5.73%  -0.42% / 5.96% -0.38% / 5.75%  -0.43% / 5.85%  -0.47% / 6.00%
## 10:  median bias MUE -0.38% / 0.43% -0.42% / -2.68% -0.38% / 0.45% -0.43% / -3.06% -0.33% / -0.25%
createTableResSim(res2stage.ar1nonbinding, xtable = TRUE)

## *** ar 2 binding
res2stage.ar2binding <- res2stage.red[infoBias==0 & binding==TRUE & ar==2,.SD,.SDcols=keep.col]

## sample size
res2stage.ar2binding[, .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                         n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                         n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]"),
                         pc.info = mean(infoPC)), by = "type"]
##        type   n.patients    n.outcome n.missing pc.info
##      <fctr>       <char>       <char>    <char>   <num>
## 1:  interim 193[178;213] 119[109;129] 74[56;93]  0.5670
## 2: decision 208[190;230] 186[169;206] 22[13;28]  0.8239
## 3:    final 265[265;265] 237[237;237] 28[28;28]  1.0163

## number of simulations
res2stage.ar2binding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 20000

## table
table2stage.ar2binding <- createTableResSim(res2stage.ar2binding, xtable = FALSE)
table2stage.ar2binding
##            statistic        method 1   method 1 fixC        method 2   method 2 fixC        method 3
##               <fctr>          <char>          <char>          <char>          <char>          <char>
##  1:     type 1 error           2.64%           2.39%           2.63%           2.42%           2.49%
##  2:            power          90.89%          90.30%          90.80%          90.38%          90.20%
##  3:       CI [NA,NA]               0               0               0               0               0
##  4:         coverage          94.81%          95.53%          94.78%          95.53%          95.23%
##  5:         reversal     0.75%/0.15%     0.47%/0.46%     0.88%/0.16%     0.44%/0.46%         0/0.85%
##  6:         abnormal         0.59%/0             0/0         0.56%/0             0/0         0/0.31%
##  7:   mean bias LMME  -0.019 / 0.017  -0.019 / 0.017  -0.019 / 0.017  -0.020 / 0.018  -0.020 / 0.017
##  8:    mean bias MUE  -0.003 / 0.004 -0.056 / -0.030  -0.003 / 0.004 -0.059 / -0.031 -0.028 / -0.017
##  9: median bias LMME  -1.40% / 2.17%  -1.40% / 2.17%  -1.34% / 2.16%  -1.54% / 2.20%  -1.65% / 1.98%
## 10:  median bias MUE -0.46% / -1.03% -7.46% / -3.52% -0.47% / -1.02% -7.33% / -3.57% -3.70% / -2.15%
createTableResSim(res2stage.ar2binding, xtable = TRUE)

## *** ar 2 non-binding
res2stage.ar2nonbinding <- res2stage.red[infoBias==0 & binding==FALSE & ar==2,.SD,.SDcols=keep.col]

## sample size
res2stage.ar2nonbinding[, .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                            n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                            n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]"),
                            pc.info = mean(infoPC)), by = "type"]
##        type   n.patients    n.outcome  n.missing pc.info
##      <fctr>       <char>       <char>     <char>   <num>
## 1:  interim 194[178;215] 120[111;131] 74[55;100]  0.5641
## 2: decision 209[192;233] 187[171;209]  22[12;28]  0.8248
## 3:    final 268[268;268] 240[240;240]  28[28;28]  1.0278

## number of simulations
res2stage.ar2nonbinding[,.N,by=c("method.char","type","hypo")][type=="interim",unique(N)]
## [1] 20000

## table
table2stage.ar2nonbinding <- createTableResSim(res2stage.ar2nonbinding, xtable = FALSE)
table2stage.ar2nonbinding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC        method 3
##               <fctr>         <char>          <char>         <char>          <char>          <char>
##  1:     type 1 error          2.48%           2.36%          2.46%           2.37%           2.49%
##  2:            power         91.10%          90.66%         91.00%          90.75%          90.67%
##  3:       CI [NA,NA]          2.78%           3.22%          3.03%           3.03%           3.35%
##  4:         coverage         95.64%          96.67%         95.67%          96.71%          96.66%
##  5:         reversal    0.68%/0.15%     0.45%/0.37%    0.78%/0.16%     0.41%/0.37%         0/0.73%
##  6:         abnormal        0.44%/0             0/0        0.40%/0             0/0         0/0.29%
##  7:   mean bias LMME  0.000 / 0.037  -0.001 / 0.039  0.000 / 0.039  -0.001 / 0.038  -0.001 / 0.038
##  8:    mean bias MUE  0.000 / 0.025 -0.001 / -0.005 -0.000 / 0.026 -0.001 / -0.007   0.002 / 0.008
##  9: median bias LMME -0.46% / 3.13%  -0.52% / 3.37% -0.47% / 3.26%  -0.53% / 3.29%  -0.62% / 3.24%
## 10:  median bias MUE -0.46% / 0.30% -0.52% / -2.15% -0.47% / 0.44% -0.53% / -2.29% -0.22% / -0.92%
createTableResSim(res2stage.ar2nonbinding, xtable = TRUE)

## * 3 stages

## ** Load data
res3stage <- readRDS(file.path("Results-built","res3stage.rds"))
res3stage[, method.char := paste0("method ",method, c(""," fixC")[fixC+1])]
res3stage[, stage.char := factor(stage, 1:3, c("interim1","interim2","final"))]
res3stage[, truth := ifelse(hypo=="power",1,0)]

## method 3 fixC same as method 3
res3stage.red <- res3stage[method.char != "method 3 fixC" & missing==TRUE]
res3stage.red[, type := factor(type, c("interim","decision","final"))]
setkeyv(res3stage.red, c("scenario","method","stage","type"))

## ** Generate table
keep.col <- c("scenario", "hypo", "method", "stage", "type", "statistic", "ck",
              "nX1", "nX2", "nX3", "infoPC",
              "estimate_ML", "se_ML", "p.value_ML", "lower_ML", "upper_ML",
              "estimate_MUE", "p.value_MUE", "lower_MUE", "upper_MUE",
              "decision", "reason", "method.char", "stage.char", "truth","seed")       

## *** ar 1 binding
res3stage.ar1binding <- res3stage.red[binding==TRUE & ar==1,.SD,.SDcols=keep.col]

## sample size
res3stage.ar1binding[, .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                         n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                         n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]"),
                         pc.info = mean(infoPC)), by = c("stage","type")]
##    stage     type   n.patients    n.outcome n.missing pc.info
##    <num>   <fctr>       <char>       <char>    <char>   <num>
## 1:     1  interim 125[115;139]    85[75;93] 40[27;57]  0.3894
## 2:     1 decision 132[121;149] 118[105;136]  14[5;24]  0.5247
## 3:     2  interim 192[181;208] 145[136;156] 47[32;66]  0.6360
## 4:     2 decision 199[189;214] 178[165;195] 21[12;28]  0.7655
## 5:     3    final 270[270;270] 241[241;241] 29[29;29]  1.0105

## number of simulations
res3stage.ar1binding[type!="interim",.N,by=c("method.char","hypo")][,unique(N)]
## [1] 20000

table3stage.ar1binding <- createTableResSim(res3stage.ar1binding, xtable = FALSE)
table3stage.ar1binding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC        method 3
##               <fctr>         <char>          <char>         <char>          <char>          <char>
##  1:     type 1 error          2.56%           2.36%          2.57%           2.35%           2.39%
##  2:            power         90.40%          89.88%         90.33%          90.10%          89.63%
##  3:       CI [NA,NA]          0.01%           0.01%         <0.01%           0.01%          <0.01%
##  4:         coverage         94.71%          95.82%         94.72%          95.88%          94.98%
##  5:         reversal    0.32%/0.11%     0.16%/0.46%    0.32%/0.12%     0.14%/0.46%         0/0.98%
##  6:         abnormal        0.52%/0             0/0        0.50%/0             0/0         0/0.08%
##  7:   mean bias LMME -0.070 / 0.052  -0.070 / 0.052 -0.070 / 0.052  -0.070 / 0.053  -0.070 / 0.053
##  8:    mean bias MUE -0.036 / 0.026 -0.135 / -0.025 -0.036 / 0.026 -0.137 / -0.025  -0.068 / 0.002
##  9: median bias LMME -5.05% / 5.42%  -5.05% / 5.42% -5.03% / 5.42%  -5.15% / 5.46%  -5.18% / 5.30%
## 10:  median bias MUE 0.38% / -0.42% -5.06% / -4.72% 0.39% / -0.44% -5.24% / -4.90% -1.89% / -2.66%
createTableResSim(res3stage.ar1binding, xtable = TRUE)

## *** ar 1 non-binding
res3stage.ar1nonbinding <- res3stage.red[binding==FALSE & ar==1,.SD,.SDcols=keep.col]

## sample size
res3stage.ar1nonbinding[, .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                            n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                            n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]"),
                            pc.info = mean(infoPC)), by = c("stage","type")]
##    stage     type   n.patients    n.outcome n.missing pc.info
##    <num>   <fctr>       <char>       <char>    <char>   <num>
## 1:     1  interim 126[117;140]    86[76;95] 40[26;57]  0.3869
## 2:     1 decision 133[123;147] 119[107;134]  14[5;24]  0.5264
## 3:     2  interim 195[185;209] 148[138;158] 47[33;64]  0.6423
## 4:     2 decision 202[191;215] 181[168;198] 21[13;28]  0.7700
## 5:     3    final 274[274;274] 245[245;245] 29[29;29]  1.0251

## number of simulations
res3stage.ar1nonbinding[type!="interim",.N,by=c("method.char","hypo")][,unique(N)]
## [1] 20000

table3stage.ar1nonbinding <- createTableResSim(res3stage.ar1nonbinding, xtable = FALSE)
table3stage.ar1nonbinding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC       method 3
##               <fctr>         <char>          <char>         <char>          <char>         <char>
##  1:     type 1 error          2.31%           2.21%          2.31%           2.21%          2.24%
##  2:            power         90.55%          89.92%         90.50%          90.09%         89.80%
##  3:       CI [NA,NA]          4.29%           4.93%          4.37%           4.66%          5.12%
##  4:         coverage         96.06%          97.89%         96.05%          97.89%         96.97%
##  5:         reversal    0.28%/0.12%     0.13%/0.60%    0.28%/0.12%     0.12%/0.60%        0/1.08%
##  6:         abnormal        0.64%/0             0/0        0.62%/0             0/0        0/0.09%
##  7:   mean bias LMME  0.003 / 0.091   0.002 / 0.094  0.003 / 0.092   0.002 / 0.093  0.002 / 0.095
##  8:    mean bias MUE  0.003 / 0.065   0.002 / 0.020  0.003 / 0.065   0.002 / 0.017  0.003 / 0.048
##  9: median bias LMME -0.12% / 7.72%  -0.16% / 8.11% -0.12% / 7.77%  -0.17% / 7.91% -0.22% / 8.22%
## 10:  median bias MUE -0.11% / 2.15% -0.15% / -1.92% -0.13% / 2.17% -0.18% / -2.43% -0.11% / 0.17%
createTableResSim(res3stage.ar1nonbinding, xtable = TRUE)

## *** ar 2 binding
res3stage.ar2binding <- res3stage.red[binding==TRUE & ar==2,.SD,.SDcols=keep.col]

## sample size
res3stage.ar2binding[, .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                         n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                         n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]"),
                         pc.info = mean(infoPC)), by = c("stage","type")]
##    stage     type   n.patients    n.outcome  n.missing pc.info
##    <num>   <fctr>       <char>       <char>     <char>   <num>
## 1:     1  interim 155[139;176]    85[76;93]  70[53;93]  0.4092
## 2:     1 decision 170[153;191] 152[132;173]   18[9;27]  0.6661
## 3:     2  interim 222[207;243] 145[136;157] 77[59;102]  0.6569
## 4:     2 decision 237[220;256] 212[193;230]  26[17;29]  0.9109
## 5:     3    final 270[270;270] 241[241;241]  29[29;29]  1.0042

## number of simulations
res3stage.ar2binding[type!="interim",.N,by=c("method.char","hypo")][,unique(N)]
## [1] 20000

table3stage.ar2binding <- createTableResSim(res3stage.ar2binding, xtable = FALSE)
table3stage.ar2binding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC        method 3
##               <fctr>         <char>          <char>         <char>          <char>          <char>
##  1:     type 1 error          2.56%           2.25%          2.55%           2.25%           2.33%
##  2:            power         90.87%          90.18%         90.80%          90.31%          90.22%
##  3:       CI [NA,NA]              0               0              0               0               0
##  4:         coverage         94.66%          95.54%         94.64%          95.66%          95.22%
##  5:         reversal    0.84%/0.26%     0.60%/0.71%    1.06%/0.29%     0.57%/0.72%         0/1.12%
##  6:         abnormal        0.69%/0             0/0        0.67%/0             0/0         0/0.33%
##  7:   mean bias LMME -0.043 / 0.030  -0.043 / 0.030 -0.042 / 0.029  -0.043 / 0.030  -0.043 / 0.031
##  8:    mean bias MUE -0.021 / 0.016 -0.086 / -0.020 -0.020 / 0.016 -0.089 / -0.021 -0.051 / -0.007
##  9: median bias LMME -3.55% / 2.75%  -3.55% / 2.75% -3.52% / 2.74%  -3.63% / 2.77%  -3.69% / 2.84%
## 10:  median bias MUE 0.32% / -0.36% -3.40% / -3.04% 0.34% / -0.38% -3.65% / -3.13% -1.90% / -2.64%
createTableResSim(res3stage.ar2binding, xtable = TRUE)


## *** ar 2 non-binding
res3stage.ar2nonbinding <- res3stage.red[binding==FALSE & ar==2,.SD,.SDcols=keep.col]

## sample size
res3stage.ar2nonbinding[, .(n.patients = paste0(median(nX1),"[",min(nX1),";",max(nX1),"]"),
                            n.outcome = paste0(median(nX3),"[",min(nX3),";",max(nX3),"]"),
                            n.missing = paste0(median(nX1-nX3),"[",min(nX1-nX3),";",max(nX1-nX3),"]"),
                            pc.info = mean(infoPC)), by = c("stage","type")]
##    stage     type   n.patients    n.outcome n.missing pc.info
##    <num>   <fctr>       <char>       <char>    <char>   <num>
## 1:     1  interim 156[141;178]    86[76;95] 70[51;91]  0.4062
## 2:     1 decision 171[155;193] 153[135;176]  18[9;26]  0.6647
## 3:     2  interim 225[206;244] 148[138;158] 77[57;99]  0.6635
## 4:     2 decision 240[223;267] 215[196;240] 26[18;29]  0.9304
## 5:     3    final 274[274;274] 245[245;245] 29[29;29]  1.0180

## number of simulations
res3stage.ar2nonbinding[type!="interim",.N,by=c("method.char","hypo")][,unique(N)]
## [1] 20000

table3stage.ar2nonbinding <- createTableResSim(res3stage.ar2nonbinding, xtable = FALSE)
table3stage.ar2nonbinding
##            statistic       method 1   method 1 fixC       method 2   method 2 fixC       method 3
##               <fctr>         <char>          <char>         <char>          <char>         <char>
##  1:     type 1 error          2.14%           2.04%          2.10%           2.03%          2.25%
##  2:            power         90.98%          90.25%         90.92%          90.38%         90.33%
##  3:       CI [NA,NA]          3.77%           4.50%          4.01%           4.25%          4.70%
##  4:         coverage         96.07%          97.49%         96.08%          97.53%         97.27%
##  5:         reversal    0.88%/0.26%     0.60%/0.71%    1.06%/0.29%     0.65%/0.71%        0/1.13%
##  6:         abnormal        0.73%/0             0/0        0.71%/0             0/0        0/0.40%
##  7:   mean bias LMME  0.000 / 0.059  -0.000 / 0.063  0.000 / 0.060  -0.000 / 0.061 -0.001 / 0.062
##  8:    mean bias MUE  0.000 / 0.047  -0.001 / 0.018 -0.000 / 0.049  -0.001 / 0.015  0.003 / 0.032
##  9: median bias LMME -0.10% / 5.08%  -0.15% / 5.50% -0.10% / 5.21%  -0.14% / 5.38% -0.26% / 5.44%
## 10:  median bias MUE -0.11% / 2.00% -0.16% / -0.10% -0.12% / 2.14% -0.15% / -0.29%  0.08% / 0.03%
createTableResSim(res3stage.ar2nonbinding, xtable = TRUE)


##----------------------------------------------------------------------
### table1.R ends here
