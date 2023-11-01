### FIGURE-pipeline.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  7 2023 (14:59) 
## Version: 
## Last-Updated: okt 31 2023 (16:38) 
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

library(DelayedGSD)
library(ggplot2)

df.sim <- GenData(n = 549,
                  N.fw = 2,
                  rand.block = c(1,1,0,0),
                  allsd = c(2.5,2.1,2.4),
                  mean0 = c(10,0,0),
                  delta = c(0,0.3,0.6),
                  ar = (0.86*2)*2*10,
                  cor.01.1 = -0.15,
                  cor.ij.1 = 0.68,
                  cor.0j.1 = -0.27,
                  seed = 1,
                  MissProb = matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538), 
                                    nrow = 2, 
                                    ncol = 2, 
                                    dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")) 
                                    ),
                  DigitsOutcome = 2,
                  TimeFactor = 14,
                  DigitsTime = 0)$d


Utime <- sort(unique(c(df.sim$t1,df.sim$t2,df.sim$t3)))

df.pipeline <- do.call(rbind,lapply(Utime, function(iTime){
    rbind(data.frame(time = iTime,
                     timepoint = "pipeline data: baseline only",
                     n = sum(df.sim$t1 <= iTime),
                     pc = mean(df.sim$t1 <= iTime)),
          data.frame(time = iTime,
                     timepoint = "pipeline data: baseline and first follow-up",
                     n = sum(df.sim$t2 <= iTime),
                     pc = mean(df.sim$t2 <= iTime)),
          data.frame(time = iTime,
                     timepoint = "complete data",
                     n = sum(df.sim$t3 <= iTime),
                     pc = mean(df.sim$t3 <= iTime))
          )
}))
df.pipeline$timepoint <- factor(df.pipeline$timepoint,
                                levels = unique(df.pipeline$timepoint))

gg.pip <- ggplot(df.pipeline, aes(x = time, y = pc, group = timepoint, color = timepoint))
gg.pip <- gg.pip + geom_line(linewidth = 2) + scale_y_continuous(labels = scales::percent)
gg.pip <- gg.pip + labs(x = "time since start of the study (in days)", y = "Percentage of patients", color = "")
gg.pip <- gg.pip + guides(color = guide_legend(nrow = 2, byrow = TRUE))
gg.pip <- gg.pip + theme(legend.position="bottom",
                         text = element_text(size=15),
                         axis.title = element_text(size=12), 
                         axis.line = element_line(linewidth = 1.25),
                         axis.ticks = element_line(linewidth = 2),
                         axis.ticks.length=unit(.25, "cm"),
                         legend.key.size = unit(2,"line"))

ggsave(gg.pip, filename = "figures/illustration-pipeline.pdf",
       width = 7.25, height = 5)

##----------------------------------------------------------------------
### FIGURE-pipeline.R ends here


       
