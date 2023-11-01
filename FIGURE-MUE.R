### FIGURE-MUE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 19 2023 (16:01) 
## Version: 
## Last-Updated: okt 31 2023 (16:37) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)
library(ggplot2)
library(xtable)

## * 2 stages

## ** Load data
res2stage <- readRDS(file.path("Results-built","res2stage.rds"))

res2stage[, method.char := paste0("method ",method, c(""," (constrained c\\u2085)")[fixC+1])]
res2stage$hypo.char <- factor(res2stage$hypo,
                              levels = c("typeI","power"),
                              labels = c(expression(paste("Under",~~H[0])),
                                         expression(paste("Under",~~H[1]))))
res2stage$method.num <-  res2stage$method+0.5*res2stage$fixC
res2stage$method.char <-  factor(res2stage$method.num,
                                 levels = c(1,1.5,2,2.5,3,3.5),
                                 labels=c(expression(paste("method 1")),
                                          expression(paste("method 1 (constrained ",c[k],")")),
                                          expression(paste("method 2")),
                                          expression(paste("method 2 (constrained ",c[k],")")),
                                          expression(paste("method 3")),
                                          expression(paste("method 3 (constrained ",c[k],")")))
                                 )
res2stage[, truth := ifelse(hypo=="power",0.6,0)]
res2stage[, stage.char := factor(stage, 1:2, c("interim","final"))]

## ** Generate figure

## *** ar 5 binding
res2stage.ar5binding <- res2stage[method.num != 3.5 & missing==TRUE & binding==TRUE & ar==5 & decision %in% c("futility","efficacy")]
## res2stage.ar5binding[]

gg2stageBeta.ar5binding <- ggplot(res2stage.ar5binding, aes(x = estimate_MUE, fill = stage.char, group = stage.char))
gg2stageBeta.ar5binding <- gg2stageBeta.ar5binding + geom_density(alpha=0.25) + facet_grid(hypo.char~method.char, labeller = label_parsed)
gg2stageBeta.ar5binding <- gg2stageBeta.ar5binding + labs(x = "estimate", fill = "stage", y = "Density (arbitrary unit)")
gg2stageBeta.ar5binding <- gg2stageBeta.ar5binding + theme(text = element_text(size=20), 
                                                 axis.line = element_line(linewidth = 1.25),
                                                 axis.ticks = element_line(linewidth = 2),
                                                 axis.ticks.length=unit(.25, "cm"),
                                                 legend.key.size = unit(3,"line"),
                                                 legend.position = "bottom")
gg2stageBeta.ar5binding

## *** ar 10 binding
res2stage.ar10binding <- res2stage[method.num != 3.5 & missing==TRUE & binding==TRUE & ar==10 & decision %in% c("futility","efficacy")]
## res2stage.ar5binding[]

gg2stageBeta.ar10binding <- ggplot(res2stage.ar10binding, aes(x = estimate_MUE, fill = stage.char, group = stage.char))
gg2stageBeta.ar10binding <- gg2stageBeta.ar10binding + geom_density(alpha=0.25) + facet_grid(hypo~method.char, labeller = label_parsed)
gg2stageBeta.ar10binding <- gg2stageBeta.ar10binding + labs(x = "estimate", fill = "stage", y = "Density (arbitrary unit)")
gg2stageBeta.ar10binding <- gg2stageBeta.ar10binding + theme(text = element_text(size=15), 
                                                 axis.line = element_line(linewidth = 1.25),
                                                 axis.ticks = element_line(linewidth = 2),
                                                 axis.ticks.length=unit(.25, "cm"),
                                                 legend.key.size = unit(3,"line"),
                                                 legend.position = "bottom")
gg2stageBeta.ar10binding

## *** ar 5 nonbinding
res2stage.ar5nonbinding <- res2stage[method.num != 3.5 & missing==TRUE & binding==FALSE & ar==5 & decision %in% c("futility","efficacy")]
## res2stage.ar5nonbinding[]

gg2stageBeta.ar5nonbinding <- ggplot(res2stage.ar5nonbinding, aes(x = estimate_MUE, fill = stage.char, group = stage.char))
gg2stageBeta.ar5nonbinding <- gg2stageBeta.ar5nonbinding + geom_density(alpha=0.25) + facet_grid(hypo~method.char, labeller = label_parsed)
gg2stageBeta.ar5nonbinding <- gg2stageBeta.ar5nonbinding + labs(x = "estimate", fill = "stage", y = "Density (arbitrary unit)")
gg2stageBeta.ar5nonbinding <- gg2stageBeta.ar5nonbinding + theme(text = element_text(size=15), 
                                                 axis.line = element_line(linewidth = 1.25),
                                                 axis.ticks = element_line(linewidth = 2),
                                                 axis.ticks.length=unit(.25, "cm"),
                                                 legend.key.size = unit(3,"line"),
                                                 legend.position = "bottom")
gg2stageBeta.ar5nonbinding

## *** ar 10 nonbinding
res2stage.ar10nonbinding <- res2stage[method.num != 3.5 & missing==TRUE & binding==FALSE & ar==10 & decision %in% c("futility","efficacy")]
## res2stage.ar5nonbinding[]

gg2stageBeta.ar10nonbinding <- ggplot(res2stage.ar10nonbinding, aes(x = estimate_MUE, fill = stage.char, group = stage.char))
gg2stageBeta.ar10nonbinding <- gg2stageBeta.ar10nonbinding + geom_density(alpha=0.25) + facet_grid(hypo~method.char, labeller = label_parsed)
gg2stageBeta.ar10nonbinding <- gg2stageBeta.ar10nonbinding + labs(x = "estimate", fill = "stage", y = "Density (arbitrary unit)")
gg2stageBeta.ar10nonbinding <- gg2stageBeta.ar10nonbinding + theme(text = element_text(size=15), 
                                                 axis.line = element_line(linewidth = 1.25),
                                                 axis.ticks = element_line(linewidth = 2),
                                                 axis.ticks.length=unit(.25, "cm"),
                                                 legend.key.size = unit(3,"line"),
                                                 legend.position = "bottom")
gg2stageBeta.ar10nonbinding

## * export
ggsave(gg2stageBeta.ar5binding, filename = file.path("figures","figure-article-simulation-MUE.pdf"), width = 16, height = 7.5)

##----------------------------------------------------------------------
### FIGURE-MUE.R ends here
