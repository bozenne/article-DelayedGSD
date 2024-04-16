rm(list=ls())

## * User interface
## cd /projects/biostat01/people/hpl802/DelayedGSD/
args <- commandArgs(TRUE) ## BATCH MODE

## arguments missing, binding, ... in BATCH model (e.g. when running on the server via slurm)
iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
if(length(args)>0){
    for (arg in args){
        eval(parse(text=arg))
    }
}
if(is.na(iter_sim)){ ## arguments for interactive R session (when not running on the server via slurm, iter_sim will be NA)
    iter_sim <- 19
    n.iter_sim <- 100

    if("missing" %in% ls() == FALSE){ missing <- TRUE }
    if("binding" %in% ls() == FALSE){ binding <- TRUE }
    if("cNotBelowFixedc" %in% ls() == FALSE){ cNotBelowFixedc <- TRUE }
    if("ar.factor" %in% ls() == FALSE){ ar.factor <- 5 }
    if("delta.factor" %in% ls() == FALSE){ delta.factor <- 0 }
    if("n.method" %in% ls() == FALSE){ n.method <- 3 }
}

name <- ""
if(missing>0){
    name <- paste(name,"missing",sep="_")
}else{
    name <- paste(name,"nomissing",sep="_")
}
if(cNotBelowFixedc>0){
    name <- paste(name,"fixC",sep="_")
}
if(binding>0){
    name <- paste(name,"binding",sep="_")
}else{
    name <- paste(name,"nonbinding",sep="_")
}
name <- paste0("2stage",name,"_ar",ar.factor)
if(delta.factor>0){
    name <- paste(name,"power",sep="_")
}else{
    name <- paste(name,"typeI",sep="_")
}

cat("BATCH ",name,": ",iter_sim," over ",n.iter_sim,"\n",sep="")
cat("Arguments:\n")
df.args <- data.frame(missing = missing,
                      binding = binding,
                      cNotBelowFixedc = cNotBelowFixedc,
                      ar.factor = ar.factor,
                      delta.factor = delta.factor)
if(!is.null(n.method)){df.args$n.method <- n.method}
print(df.args, row.names = FALSE)
cat("\n")

## * Settings
nsim <- 100 # number of simulations
method <- 1:3 # methods used to compute the boundaries
                                        #--- to plan the trial ----
kMax <- 2  #max number of analyses (including final)
alpha <- 0.025  #type I error (one sided)
beta <- 0.1  #type II error
informationRates <- c(0.53,1)  #planned  information rates
rho_alpha <- 2  # rho parameter for alpha error spending function
rho_beta <- 2  # rho parameter for beta error spending function
## deltaPower <- 0.75 # just to try another value when Id > Imax
Id <- 0.63  #(expected) information rate at each decision analysis
                                        #
                                        #---- to generate data -----------
                                        #
block <- c(1,1,0,0) 
allsd <- c(2.5,2.1,2.4) # sd, first from baseline measurement, then the two changes from baseline
mean0 <- c(10,0,0) # mean placebo group (again, first is absolute value, then change from baseline)
delta <- c(0,0.5,1)*delta.factor # treatment effect
ar <- 2*2.5*ar.factor # orginial accrual rate from data from Corine is 0.86 per week, hence we multiply by 2 for by 14 days. As too low, we further multiply by 2
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
PropForInterim <- 0.5 # Decide to have interim analysiz when PropForInterim % of all subjects have had the chance to have one follow-up measuement recorded in the data to be available for analysis.
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

## * Load dependencies
library(DelayedGSD) ## remotes::install_github("PauloWhite/DelayedGSD")
DelayedGSD.options(FCT.p_value = "FinalPvalue2", continuity.correction = TRUE)
source("FCT.R") ## exportGSD function

## * Planned boundaries
plannedB <- vector(mode = "list", length = 3)
for(iMeth in method){ ## iMeth <- 1
    plannedB[[iMeth]] <- CalcBoundaries(kMax=kMax,  
                                        alpha=alpha, 
                                        beta=beta,  
                                        InfoR.i=informationRates,  
                                        InfoR.d=c(Id,1),  
                                        rho_alpha=rho_alpha,  
                                        rho_beta=rho_beta,  
                                        method=iMeth,  
                                        cNotBelowFixedc=cNotBelowFixedc || (iMeth==3),
                                        bindingFutility=binding,
                                        delta=deltaPower)
    ## summary(plannedB[[1]])
    ## coef(plannedB[[iMeth]], type = "information")
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
allj <- seq(1+(iter_sim-1)*nsim, iter_sim*nsim, by = 1)
#allj <- 572:1000
for(j in allj){ ## j <- 1 ## 5
  startComp <- Sys.time()
  myseedi <- allseeds[j]
  #myseedi <- 402321297
  # {{{ TRACE info (e.g. to check the Rout)
  cat("seed ",myseedi," for ","j=",j," (index ",which(j==allj),") out of ",nsim,": ", sep="")
  # }}}
  
  # {{{ generate data
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
  ## head(d,n=20)
  # }}}
  # {{{ reformat data like those of Corine
  ## Make data long format
  ## dd <- FormatAsCase(d)
  ## head(dd)
  ## summary(dd)
  # }}}
  
  # {{{ make data available at interim
  # Here we stop inclusion data collection for the interim analysis as soon as
  # half of the participants have completed (or had the opportunity to complete) the follow-up 
  thets <- d$t3[ceiling(nGSD*PropForInterim)]
  cat("time = ",paste(thets, collapse = ", "),"\n",sep="")
  #thet <- d$t3[ceiling(n*PropForInterim)]
  
  ## ddi <- FormatAsCase(di) # needed ????
  ## head(d[d$id==52,])
  # }}}
  
  nX1.interim <- vector()
  nX2.interim <- vector()
  nX3.interim <- vector()
  currentGSD <- vector(mode = "list", length = 3)
  out.interim <- vector(mode = "list", length = 3)
  for(iMeth in method){ ## iMeth <- 3
    # {{{ make data available at interim
    di <- SelectData(d,t=thets[iMeth])
    
    nX1.interim[iMeth] <- sum(!is.na(di$X1))
    nX2.interim[iMeth] = sum(!is.na(di$X2))
    nX3.interim[iMeth] = sum(!is.na(di$X3))
    
    ## {{{ analyze data at at interim
    ## ** interim
    lmmI <- analyzeData(di, ddf = "nlme", data.decision = sum(d$t1 <= thets[iMeth] + theDelta.t*TimeFactor), getinfo = TRUE, trace = TRUE)
    ## lmmI <- analyzeData(di, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    
    currentGSD[[iMeth]] <- update(plannedB[[iMeth]], delta = lmmI, trace = FALSE)
    
    out.interim[[iMeth]] <- exportGSD(currentGSD[[iMeth]],
                                      export.statistic = TRUE,
                                      export.ML = TRUE,
                                      export.MUE = FALSE,
                                      export.info = TRUE,
                                      export.predinfo = TRUE,
                                      export.boundary = TRUE,
                                      export.decision = TRUE)
  }
  ## currentGSD[[3]]
  ## plot(currentGSD[[1]], legend.x = "bottomleft")
  
  
  out.decision <- vector(mode = "list", length = 3)
  for(iMeth in method){ ## iMeth <- 1
    
    ## ** decision
    dDecision <- d[which(d$t1 <= thets[iMeth] + theDelta.t*TimeFactor),]
    lmmD <- analyzeData(dDecision, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    
      ## Non binding: never stop for futility when simulating under the null and always stop for futility when simulating under the alternative
      ## (then the observed rejection rate should match the nominal type 1 or type 2 error)
      if(out.interim[[iMeth]]$decision == "stop" && (out.interim[[iMeth]]$reason!="futility" || binding == TRUE || delta.factor > 0)){

          currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmmD, trace = FALSE)
          ## plot(currentGSD[[iMeth]])
      
          out.decision[[iMeth]] <- exportGSD(currentGSD[[iMeth]],
                                             export.statistic = TRUE,
                                             export.ML = TRUE,
                                             export.MUE = TRUE,
                                             export.info = TRUE,
                                             export.predinfo = FALSE,
                                             export.boundary = TRUE,
                                             export.decision = TRUE)
      
      
      }else{
          if(out.interim[[iMeth]]$decision == "stop"){ ## overrule futility boundary
              currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], overrule.futility = TRUE)
          }
          ## update information
          currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmmD, k = 1, type.k = "decision", trace = FALSE)
      }
  }
                                        # }}}
                                        # {{{ Analyze data at decision
  
  ## ** finale
  out.final <- vector(mode = "list", length = 3)
  for(iMeth in method){ ## iMeth <- 3
    dFinal <- d[1:nGSD[iMeth],]
    lmmF <- analyzeData(dFinal, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    
    if(currentGSD[[iMeth]]$stage[,"type"]=="decision"){
      
      out.final[[iMeth]] <- cbind(method = iMeth, stage = 2, type = "final", exportGSD(NA))
      
    }else{
      currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmmF, trace = FALSE)
      
      out.final[[iMeth]] <- exportGSD(currentGSD[[iMeth]],
                                      export.statistic = TRUE,
                                      export.ML = TRUE,
                                      export.MUE = TRUE,
                                      export.info = TRUE,
                                      export.predinfo = FALSE,
                                      export.boundary = TRUE,
                                      export.decision = TRUE)
      
    }
  }
  # }}}
  
  stopComp <- Sys.time()
  # {{{ Save results
  outMerge <- do.call(rbind,lapply(method, function(iMeth){ ## iMeth <- 3
      cbind(rbind(out.interim[[iMeth]],
                  out.decision[[iMeth]],
                  out.final[[iMeth]]),
            time.interim = thets[iMeth], nX1.interim = nX1.interim[iMeth], nX2.interim = nX2.interim[iMeth], nX3.interim = nX3.interim[iMeth])
  }))
  
  ## outMerge[outMerge$method==3,]
  
  out <- cbind(
      ## results
      outMerge,
      ## simulation details
      seed=myseedi,             
      ## computation time
      computation.time=as.double(round(difftime(stopComp,startComp,units="secs"),3))
  )
  out$sample.size <- list(nGSD)
  ## names(out) <- myColNames
  RES <- rbind(RES,out)
  if(j %in% round(quantile(allj, probs = (1:10)/10))){
      saveRDS(RES,file=file.path("Results",name,paste0("sim-",name,"-",iter_sim,"(tempo)_",nsim,".rds")))
  }
  # }}}
}

## * Export
rownames(RES) <- NULL
saveRDS(RES,file=file.path("Results",name,paste0("sim-",name,"-",iter_sim,"_",nsim,".rds")))

sessionInfo()
