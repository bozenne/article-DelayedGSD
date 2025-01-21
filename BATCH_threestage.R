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
    iter_sim <- 132
    n.iter_sim <- 200

    if("missing" %in% ls() == FALSE){ missing <- TRUE }
    if("binding" %in% ls() == FALSE){ binding <- TRUE }
    if("cNotBelowFixedc" %in% ls() == FALSE){ cNotBelowFixedc <- TRUE }
    if("ar.factor" %in% ls() == FALSE){ ar.factor <- 1 }
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
name <- paste0("3stage",name,"_ar",ar.factor)
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
kMax <- 3  #max number of analyses (including final)
alpha <- 0.025  #type I error (one sided)
beta <- 0.1  #type II error
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
theDelta.t <- 3*7 # number of days to process the data and make them ready to analyze after collecting them
TimeFactor <- 6*7 ## number of days between two visits
ar <- (TimeFactor/7)*2.5*ar.factor ## number of subjects accrued between two visits (orginial accrual rate from data from Corine is 0.86 per week, this about three times larger)
                                        #
                                        #--- actually for both planing the trial  and generating data-----
                                        #
                                        #
deltaPower <- 1 # effect (NOT Z-scale/unit, but outcome scale/unit!) that the study is powered for: should we choose ourselves or compute from other numbers above ???
sdPower <- allsd[3]*sqrt(1-cor0j1^2)
n <- 2*2*((sdPower/deltaPower)^2)*(qnorm(1-beta)-qnorm(alpha))^2 #104 with Corine's data # should we choose ourselves or compute from the above numbers ???
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
DelayedGSD.options(FCT.p_value = "FinalPvalue2", continuity.correction = TRUE, tolerance = 2.5e-3) ## tolerance for CI/MUE = alpha/10
source("FCT.R") ## exportGSD function

## * Planned boundaries
plannedB <- vector(mode = "list", length = length(method))
for(iMeth in 1:length(method)){ ## iMeth <- 1
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
allj <- seq(1+(iter_sim-1)*nsim, iter_sim*nsim, by = 1)
#allj <- 572:1000
for(j in allj){ ## j <- 1 ## 5
  startComp <- Sys.time()
  myseedi <- allseeds[j]
  # {{{ TRACE info (e.g. to check the Rout)
  print(paste0("seed ",myseedi," for ","j=",j," (index ",which(j==allj),") out of ",nsim))
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


  ## ** prepare export
  currentGSD.interim <- setNames(lapply(1:(kMax-1), function(i){
      setNames(vector(mode = "list", length = length(method)), paste0("method ",method))
  }), paste0("decision ", 1:(kMax-1)))
  currentGSD.decision <- setNames(lapply(1:(kMax-1), function(i){
      setNames(vector(mode = "list", length = length(method)), paste0("method ",method))
  }), paste0("interim ", 1:(kMax-1)))
  currentGSD.final <- setNames(vector(mode = "list", length = length(method)), paste0("method ",method))
  
  outMerge <- NULL
 
  ## ** early stop
  stopGSD <- setNames(rep(FALSE,length(method)), paste0("method ",method))
  thets <- matrix(NA, nrow = length(method), ncol = kMax-1,
                  dimnames = list(NULL, paste0("time.interim",1:(kMax-1))))
  
  for(iStage in 1:(kMax-1)){ ## iStage <- 1
                                        
      for(iMeth in which(stopGSD==FALSE)){ ## iMeth <- 1
          thets[iMeth,iStage] <- d$t3[ceiling(nGSD[iMeth]*PropForInterim[iStage])]

                                        # {{{ make data available at interim
          di <- SelectData(d,t=thets[iMeth,iStage])

          ## *** interim
          lmmI <- analyzeData(di, ddf = "nlme", data.decision = sum(d$t1 <= thets[iMeth,iStage] + theDelta.t), getinfo = TRUE, trace = TRUE)
          if(iStage == 1){
              currentGSD.interim[[iStage]][[iMeth]] <- update(plannedB[[iMeth]], delta = lmmI, trace = FALSE)
          }else{
              currentGSD.interim[[iStage]][[iMeth]] <- update(currentGSD.interim[[iStage-1]][[iMeth]], delta = lmmI, trace = FALSE)
          }
          
          outMerge <- rbind(outMerge,
                            cbind(time = thets[iMeth,iStage],
                                  nX1 = sum(!is.na(di$X1)),
                                  nX2 = sum(!is.na(di$X2)),
                                  nX3 = sum(!is.na(di$X3)),
                                  exportGSD(currentGSD.interim[[iStage]][[iMeth]],
                                            export.statistic = TRUE,
                                            export.ML = TRUE,
                                            export.MUE = FALSE,
                                            export.info = TRUE,
                                            export.predinfo = TRUE,
                                            export.boundary = TRUE,
                                            export.decision = TRUE))
                            )
    
          ## *** decision
          dDecision <- d[which(d$t1 <= thets[iMeth,iStage] + theDelta.t),]
          lmmD <- analyzeData(dDecision, ddf = "nlme", getinfo = TRUE, trace = TRUE)

          ## Non binding: never stop for futility when simulating under the null and always stop for futility when simulating under the alternative
          ## (then the observed rejection rate should match the nominal type 1 or type 2 error)
          iConclusion <- currentGSD.interim[[iStage]][[iMeth]]$conclusion
          if(iConclusion["interim",iStage] == "stop" && (iConclusion["reason.interim",iStage]!="futility" || binding == TRUE || delta.factor > 0)){

              currentGSD.decision[[iStage]][[iMeth]] <- update(currentGSD.interim[[iStage]][[iMeth]], delta = lmmD, trace = FALSE)
              ## plot(currentGSD[[iMeth]])
      
              countNA <- rowSums(is.na(dDecision[,paste0("X",2:3)]))

              outMerge <- rbind(outMerge,
                                cbind(time = max(dDecision$t3),
                                      nX1 = sum(!is.na(dDecision$X1)),
                                      nX2 = sum(!is.na(dDecision$X2)),
                                      nX3 = sum(!is.na(dDecision$X3)),
                                      exportGSD(currentGSD.decision[[iStage]][[iMeth]],
                                                export.statistic = TRUE,
                                                export.ML = TRUE,
                                                export.MUE = TRUE,
                                                export.info = TRUE,
                                                export.predinfo = FALSE,
                                                export.boundary = TRUE,
                                                export.decision = TRUE))
                                )
              stopGSD[iMeth] <- TRUE              
      
          }else{
              if(iConclusion["interim",iStage] == "stop"){ ## overrule futility boundary                  
                  currentGSD.interim[[iStage]][[iMeth]] <- update(currentGSD.interim[[iStage]][[iMeth]], overrule.futility = TRUE)
              }
              ## update information
              currentGSD.interim[[iStage]][[iMeth]] <- update(currentGSD.interim[[iStage]][[iMeth]], delta = lmmD, k = iStage, type.k = "decision", trace = FALSE)
          }
      }
  }
                                        # }}}
                                        # {{{ Analyze data at decision
  ## ** finale
  for(iMeth in which(stopGSD==FALSE)){ ## iMeth <- 1
      dFinal <- d[1:nGSD[iMeth],]
      lmmF <- analyzeData(dFinal, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    
      currentGSD.final[[iMeth]] <- update(currentGSD.interim[[iStage]][[iMeth]], delta = lmmF, trace = FALSE)
      
      countNA <- rowSums(is.na(dFinal[,paste0("X",2:3)]))

      ## confint(currentGSD.final[[iMeth]], method = c("ML","MUE"))
      outMerge <- rbind(outMerge,
                        cbind(time = max(dFinal$t3),
                              nX1 = sum(!is.na(dFinal$X1)),
                              nX2 = sum(!is.na(dFinal$X2)),
                              nX3 = sum(!is.na(dFinal$X3)),
                              exportGSD(currentGSD.final[[iMeth]],
                                        export.statistic = TRUE,
                                        export.ML = TRUE,
                                        export.MUE = TRUE,
                                        export.info = TRUE,
                                        export.predinfo = FALSE,
                                        export.boundary = TRUE,
                                        export.decision = TRUE)
                              ))
      
  }
                                        # }}}
  stopComp <- Sys.time()
                                        # {{{ Save results
  rownames(outMerge) <- NULL
  out <- cbind(outMerge,
                     seed=myseedi,
                     computation.time=as.double(round(difftime(stopComp,startComp,units="secs"),3))
               )
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
