### example-increasing-uk-3stage.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  9 2023 (11:06) 
## Version: 
## Last-Updated: okt 26 2023 (17:14) 
##           By: Brice Ozenne
##     update #: 40
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
iter_sim <- 95
n.iter_sim <- 100

missing <- TRUE 
binding <- TRUE 
cNotBelowFixedc <- FALSE 
ar.factor <- 5 
delta.factor <- 0.6 
n.method <- 3 


## * Settings
nsim <- 100 # number of simulations
method <- 1:3 # methods used to compute the boundaries
                                        #--- to plan the trial ----
kMax <- 3  #max number of analyses (including final)
alpha <- 0.025  #type I error (one sided)
beta <- 0.2  #type II error
informationRates <- c(0.40,0.65,1)  #planned  information rates
rho_alpha <- 2  # rho parameter for alpha error spending function
rho_beta <- 2  # rho parameter for beta error spending function
## deltaPower <- 0.75 # just to try another value when Id > Imax
Id <- c(0.50,0.75)  #(expected) information rate at each decision analysis
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

PropForInterim <- c(0.35,0.6) # Decide to have interim analysiz when PropForInterim % of all subjects have had the chance to have one follow-up measuement recorded in the data to be available for analysis.
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
plannedB <- vector(mode = "list", length = length(method))
for(iMeth in 1:length(method)){ ## iMeth <- 2
  plannedB[[iMeth]] <- CalcBoundaries(kMax = kMax,  
                                      alpha = alpha, 
                                      beta = beta,  
                                      InfoR.i = informationRates,  
                                      InfoR.d = c(Id,1),  
                                      rho_alpha = rho_alpha,  
                                      rho_beta = rho_beta,  
                                      method = method[iMeth],  
                                      cNotBelowFixedc = cNotBelowFixedc || (iMeth==3),
                                      bindingFutility = binding,
                                      delta = deltaPower)
  ## summary(plannedB[[1]])
  ## plot(plannedB[[1]])
  ## coef(plannedB[[iMeth]], type = "decision")
}
if(is.null(n.method)){
    inflationFactor <- unlist(lapply(plannedB,function(iP){iP$planned$InflationFactor}))
}else{
    inflationFactor <- rep(plannedB[[n.method]]$planned$InflationFactor, 3)
}
nGSD <- ceiling(n*inflationFactor)
RES <- NULL

cat("Sample size: ",paste(nGSD, collapse = ", "),"\n",sep="")

## * Loop
myseedi <- 942188850

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
currentGSD.interim <- setNames(lapply(1:(kMax-1), function(i){
    setNames(vector(mode = "list", length = length(method)), paste0("method ",method))
}), paste0("decision ", 1:(kMax-1)))
currentGSD.decision <- setNames(lapply(1:(kMax-1), function(i){
    setNames(vector(mode = "list", length = length(method)), paste0("method ",method))
}), paste0("interim ", 1:(kMax-1)))
currentGSD.final <- setNames(vector(mode = "list", length = length(method)), paste0("method ",method))
thets <- matrix(NA, nrow = length(method), ncol = kMax-1,
                dimnames = list(NULL, paste0("time.interim",1:(kMax-1))))
    
## ** interim 1
iMeth <- 1
thets[iMeth,1] <- d$t3[nGSD[iMeth]*PropForInterim[1]]
lmmI1 <- analyzeData(SelectData(d,t=thets[iMeth,1]), ddf = "nlme", data.decision = sum(d$t1 <= thets[iMeth,1] + theDelta.t*TimeFactor), getinfo = TRUE, trace = TRUE)
currentGSD.interim[[1]][[iMeth]] <- update(plannedB[[iMeth]], delta = lmmI1, trace = FALSE)

## ** interim 2
iMeth <- 1
thets[iMeth,2] <- d$t3[nGSD[iMeth]*PropForInterim[2]]
lmmI2 <- analyzeData(SelectData(d,t=thets[iMeth,2]), ddf = "nlme", data.decision = sum(d$t1 <= thets[iMeth,2] + theDelta.t*TimeFactor), getinfo = TRUE, trace = TRUE)
currentGSD.interim[[2]][[iMeth]] <- update(currentGSD.interim[[1]][[iMeth]], delta = lmmI2, trace = FALSE)

## ** debug
iLk <- currentGSD.interim[[1]][[iMeth]]$lk[1]
iUk <- currentGSD.interim[[1]][[iMeth]]$uk[1]
iInfo.max <- currentGSD.interim[[1]][[iMeth]]$planned$Info.max
iVinfo <- c(currentGSD.interim[[1]][[iMeth]]$Info.i[1],lmmI2$information["interim"])
iRho <- sqrt(iVinfo[1]/iVinfo[2])
iMinfo <- matrix(c(1,iRho,iRho,1),2,2)
iAlpha <- ErrorSpend(iVinfo[2], rho = 2, beta_or_alpha = 0.025, Info.max = iInfo.max) - ErrorSpend(iVinfo[1], rho = 2, beta_or_alpha = 0.025, Info.max = iInfo.max)

pmvnorm(lower = c(iLk,1.1*iUk),
        upper = c(iUk,Inf),
        mean = rep(0,2),
        sigma= iMinfo,
        abseps = 1e-6)


uniroot(function(x){pmvnorm(lower = c(iLk,x),
                            upper = c(iUk,Inf),
                            mean = rep(0,2),
                            sigma= iMinfo,
                            abseps = 1e-6) - iAlpha},
        lower = iLk,
        upper = 1.1*iUk,
        tol = 1e-6)$root





##----------------------------------------------------------------------
### example-increasing-uk-3stage.R ends here
