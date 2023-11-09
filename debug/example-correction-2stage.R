### example-correction-2stage.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  9 2023 (14:42) 
## Version: 
## Last-Updated: nov  9 2023 (16:23) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(DelayedGSD)

## * User interface
iter_sim <- 19
n.iter_sim <- 100

missing <- TRUE 
binding <- TRUE 
cNotBelowFixedc <- TRUE 
ar.factor <- 5 
delta.factor <- 0
n.method <- 3 


## * Settings
nsim <- 100 # number of simulations
method <- 1:3 # methods used to compute the boundaries
                                        #--- to plan the trial ----
kMax <- 2  #max number of analyses (including final)
alpha <- 0.025  #type I error (one sided)
beta <- 0.2  #type II error
informationRates <- c(0.58,1)  #planned  information rates
rho_alpha <- 2  # rho parameter for alpha error spending function
rho_beta <- 2  # rho parameter for beta error spending function
## deltaPower <- 0.75 # just to try another value when Id > Imax
Id <- c(0.68)  #(expected) information rate at each decision analysis
                                        #
                                        #---- to generate data -----------
                                        #
block <- c(1,1,0,0) 
allsd <- c(2.5,2.1,2.4) # sd, first from baseline measurement, then the two changes from baseline
mean0 <- c(10,0,0) # mean placebo group (again, first is absolute value, then change from baseline)
delta <- c(0,0.5,1)*delta.factor # treatment effect
ar <- (0.86*2)*2*ar.factor # orginial accrual rate from data from Corine is 0.86 per week, hence we multiply by 2 for by 14 days. As too low, we further multiply by 2
cor011 <- -0.15 # ~ from data from Corine
corij1 <- 0.68  # ~ from data from Corine
cor0j1 <- -0.27  # ~ from data from Corine
if(missing){
    Miss11 <- 5/104 # miss both V1 and V2
    Miss12 <- 1/104 # miss V1 and but not V2
    Miss21 <- 6/104 # do not miss V1 and but miss V2
    Miss22 <- 92/104 # miss none
    MyMissProb <- matrix(c(Miss11,Miss12,Miss21,Miss22),ncol=2,nrow=2,byrow=TRUE, # to additionnally remove 1 more because some FASFL=N
                         dimnames = list(c("V1 missing","V1 not missing"), c("V2 missing","V2 not missing")))
}else{
    MyMissProb <- NULL
}

PropForInterim <- c(0.5) # Decide to have interim analysiz when PropForInterim % of all subjects have had the chance to have one follow-up measuement recorded in the data to be available for analysis.
theDelta.t <- 1.50001 # time lag to process the data and make them ready to analyze after collecting them (unit is time between two follow-up visits)
TimeFactor <- 14 ## number of days between two visits
                                        #
                                        #--- actually for both planing the trial  and generating data-----
                                        #
                                        #
deltaPower <- 0.6 # effect (NOT Z-scale/unit, but outcome scale/unit!) that the study is powered for: should we choose ourselves or compute from other numbers above ???
sdPower <- allsd[3]*sqrt(1-cor0j1^2)
n <- ceiling(2*2*((sdPower/deltaPower)^2)*(qnorm(1-beta)-qnorm(alpha))^2) #104 with Corine's data # should we choose ourselves or compute from the above numbers ???
                                       # inflate SS as required for interim

## adjust for expected withdrawal
if(missing){
    n <- n/(1-(Miss11+Miss21))
}

## * Seed
set.seed(140786598)
nsimAll <- n.iter_sim * nsim
allseeds <- sample.int(n = 1000000000, size = nsimAll, replace=FALSE) #x=1:(.Machine$integer.max) seems to be maximal possible

## * Planned boundaries
plannedB <- vector(mode = "list", length = 3)
plannedB[[1]] <- CalcBoundaries(kMax = kMax,  
                                alpha = alpha, 
                                beta = beta,  
                                InfoR.i = informationRates,  
                                InfoR.d = c(Id,1),  
                                rho_alpha = rho_alpha,  
                                rho_beta = rho_beta,  
                                method = 1,  
                                cNotBelowFixedc = FALSE,
                                bindingFutility = binding,
                                delta = deltaPower)
plannedB[[2]] <- CalcBoundaries(kMax = kMax,  
                                alpha = alpha, 
                                beta = beta,  
                                InfoR.i = informationRates,  
                                InfoR.d = c(Id,1),  
                                rho_alpha = rho_alpha,  
                                rho_beta = rho_beta,  
                                method = 1,  
                                cNotBelowFixedc = TRUE,
                                bindingFutility = binding,
                                delta = deltaPower)

nGSD <- 549
RES <- NULL

cat("Sample size: ",paste(nGSD, collapse = ", "),"\n",sep="")

## * Loop
myseedi <- 312274104## 402321297

## ** simulate
res <- GenData(n=max(nGSD), 
               N.fw=2,
               rand.block=block,
               allsd=allsd,
               mean0=mean0,
               delta=delta,
               ar=ar,
               cor.01.1=cor011,
               cor.ij.1=corij1,
               cor.0j.1=cor0j1,
               seed=myseedi,
               MissProb=MyMissProb,
               DigitsOutcome=2,
               TimeFactor=TimeFactor,
               DigitsTime=0
               )
d <- res$d

## ** prepare export
currentGSD.interim <- vector(mode = "list", length = 2)
currentGSD.decision <- vector(mode = "list", length = 2)
currentGSD.final <- vector(mode = "list", length = 2)
thets <- matrix(NA, nrow = length(method), ncol = kMax-1,
                dimnames = list(NULL, paste0("time.interim",1:(kMax-1))))
    
## ** interim 1
iMeth <- 1
thets[iMeth,1] <- d$t3[nGSD[iMeth]*PropForInterim[1]]
lmmI1 <- analyzeData(SelectData(d,t=thets[iMeth,1]), ddf = "nlme", data.decision = sum(d$t1 <= thets[iMeth,1] + theDelta.t*TimeFactor), getinfo = TRUE, trace = TRUE)
currentGSD.interim[[1]] <- update(plannedB[[1]], delta = lmmI1, trace = FALSE)
currentGSD.interim[[2]] <- update(plannedB[[2]], delta = lmmI1, trace = FALSE)

dDecision <- d[which(d$t1 <= thets[iMeth,1] + theDelta.t*TimeFactor),]
lmmD <- analyzeData(dDecision, ddf = "nlme", getinfo = TRUE, trace = TRUE)

GS <- update(currentGSD.interim[[1]], delta = lmmD, trace = FALSE)
summary(GS)
## estimate   lower  upper p.value
##  -0.1217 -0.6253 0.3818 0.68208

DelayedGSD.options(FCT.p_value = "FinalPvalue2", continuity.correction = 1)
test1 <- update(currentGSD.interim[[2]], delta = lmmD, trace = FALSE)
summary(test1)
## estimate  lower  upper p.value
##  -0.2085 -0.712 0.2951 0.79139

##      efficacy reversal continue
## [1,] 0.008379   0.5456   0.2374


DelayedGSD.options(FCT.p_value = "FinalPvalue", continuity.correction = 2)
test2 <- update(currentGSD.interim[[2]], delta = lmmD, trace = FALSE)
summary(test2)
## estimate   lower  upper p.value
##  -0.1217 -0.6253 0.3818 0.68208





##----------------------------------------------------------------------
### example-correction-2stage.R ends here
