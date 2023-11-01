### FCT.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  1 2022 (15:45) 
## Version: 
## Last-Updated: okt 10 2023 (13:11) 
##           By: Brice Ozenne
##     Update #: 129
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * exportGSD 
##' @description Take a GSD object and extract relevant elements
exportGSD <- function(object,
                      export.statistic = TRUE,
                      export.ML = TRUE,
                      export.MUE = TRUE,
                      export.info = TRUE,
                      export.predinfo = TRUE,
                      export.boundary = TRUE,
                      export.decision = TRUE,
                      export.sigma = TRUE){

    ## ** initialize
    out <- data.frame(statistic = NA,
                      estimate_ML = NA,
                      se_ML = NA,
                      p.value_ML = NA,
                      lower_ML = NA,
                      upper_ML = NA,
                      estimate_MUE = NA,
                      p.value_MUE = NA,
                      lower_MUE = NA,
                      upper_MUE = NA,
                      info = NA,
                      infoPC = NA,
                      info.pred = NA,
                      infoPC.pred = NA,
                      uk = NA,
                      lk = NA,
                      ck = NA,
                      decision = NA,
                      reason = NA,
                      sigma = NA)

    ## ** check user input
    if(identical(object,NA)){
       return(out) 
    }else if(is.null(object)){
        stop("Argument \'object\' is NULL. \n")
    }else if(!inherits(object,"delayedGSD")){
        stop("Argument \'object\' must inherits from \"delayedGSD\". \n")
    }

    ## ** extract information
    stage <- object$stage
    if(stage$type %in% "interim"){
        object.confint <- confint(object)
    }else  if(stage$type %in% c("decision","final")){
        object.confint <- confint(object, method = c("ML","MUE"))
    }
    object.confint <- object.confint[object.confint$stage==stage$k,,drop=FALSE]
    
    object.info <- coef(object, type = "information")
    object.info <- object.info[object.info$stage==stage$k,,drop=FALSE]

    object.boundary <- coef(object, type = "boundary")
    object.boundary <- object.boundary[object.boundary$stage==stage$k,,drop=FALSE]

    object.decision <- coef(object, type = "decision")

    ## ** fill output
    if(export.statistic){
        out$statistic <- object.confint[1,"statistic"]
    }
    if(export.ML){
        if(stage$type %in% "interim"){
            out$estimate_ML <- object.confint[,"estimate"]
            out$se_ML <- object.confint[,"se"]
        }else if(stage$type %in% c("decision","final")){
            out$estimate_ML <- object.confint[object.confint$method == "ML","estimate"]
            out$se_ML <- object.confint[object.confint$method == "ML","se"]
            out$p.value_ML <- object.confint[object.confint$method == "ML","p.value"]
            out$lower_ML <- object.confint[object.confint$method == "ML","lower"]
            out$upper_ML <- object.confint[object.confint$method == "ML","upper"]
        }
    }

    if(export.MUE){
        if(stage$type %in% c("decision","final")){
            out$estimate_MUE <- object.confint[object.confint$method == "MUE","estimate"]
            out$p.value_MUE <- object.confint[object.confint$method == "MUE","p.value"]
            out$lower_MUE <- object.confint[object.confint$method == "MUE","lower"]
            out$upper_MUE <- object.confint[object.confint$method == "MUE","upper"]
        }
    }

    if(export.info){
        if(stage$type %in% c("interim","final")){
            out$info <- object.info[,"Interim"]
            out$infoPC <- object.info[,"Interim.pc"]
        }else if(stage$type == "decision"){
            out$info <- object.info[,"Decision"]
            out$infoPC <- object.info[,"Decision.pc"]
        }
    }

    if(export.predinfo){
        if(stage$type %in% "interim"){
            out$info.pred <- object.info[,"Decision"]
            out$infoPC.pred <- object.info[,"Decision.pc"]
        }
    }

    if(export.boundary){
        if(stage$type %in% "interim"){
            out$uk <- object.boundary[,"Ebound"]
            out$lk <- object.boundary[,"Fbound"]
        }else if(stage$type %in% c("decision","final")){
            out$ck <- object.boundary[,"Cbound"]
        }
    }

    if(export.decision){
        out$decision <- object.decision["decision",NCOL(object.decision)]
        out$reason <- object.decision["comment",NCOL(object.decision)]
    }
    if(export.sigma){
        index.lmm <- utils::tail(which(sapply(object$lmm,is.null)==FALSE),1)
        out$sigma <- sigma(object$lmm[[index.lmm]]$fit)
    }

    ## ** export
    return(cbind(method = object$method, stage = object$stage[1,"k"], type = object$stage[1,"type"],out))

}

## * FCT_simGSD (documentation)
##' @param method [integer, 1-14] method to be used (see method2num)
##' @param alpha [numeric, 0-1] type I error (one sided)
##' @param beta [numeric, 0-1] max number of analyses (including final)
##' @param kMax [integer] max number of analyses (including final)
##' @param InfoR.i [numeric vector of length kMax] planned information rates at interim and final.
##' @param InfoR.d [numeric vector of length kMax] planned information rates at decision and final.
##' @param rho_alpha [numeric, >0] rho parameter for alpha error spending function
##' @param rho_beta [numeric, >0] rho parameter for beta error spending function
##' @param ar [numeric, >0] accrual rate
##' @param deltaPlan [numeric vector of length Kmax+1] planned treatment effect (i.e. difference in mean) on primary outcome at each visit,
##' used to decide on the sample size and planned boundaries. 
##' First value from baseline measurement, following for from baseline.
##' @param deltaTrue [numeric vector of length Kmax+1] (true) treatment effect (i.e. difference in mean) on primary outcome at each visit,
##' used to generate the data. 
##' First value from baseline measurement, following for from baseline.
##' If \code{NULL}, taken to be the same as the planned treatment effect.
##' @param allsd [numeric vecotr of length Kmax+1] standard deviation.
##' First value from baseline measurement, following for from baseline.
##' @param n.sample [integer] sample size used in each simulation. If \code{NULL}, will be deduced from the treatment effect, standard deviation,  type I and type II error.
##' @param mean0 [numeric vecotr of length Kmax+1] mean placebo group.
##' First value from baseline measurement, following for from baseline.
##' @param PropForInterim [numeric vector of length Kmax-1] decide to have interim analysiz when PropForInterim % of all subjects have had the chance to have one follow-up measuement recorded in the data to be available for analysis.
##' @param theDelta.t [numeric] time lag to process the data and make them ready to analyze after collecting them (unit is time between two follow-up visits).
##' @param path [character] directory where to export the results.
##' If \code{NULL}, will be exported in the current directory.
##' If \code{NA}, will not be exported.
##' @param name [character] start of the filename.
##' @param export.tempo [logical] should results be exported at each iteration in a temporary file.
##' @param n.sim [integer,>0] number of simulations.
##' Can also be of length 2, in that case the second element indicates which replicate of the simulation is being run.
##' @param trace [logical] When \code{TRUE}, ouput is generate in the console as the execution of the function progresses.
##' @param seed,N.fw,rand.block,mean0,cor.01.1,cor.ij.1,cor.0j.1,TimeFactor,MissProb arguments passed to GenData
##' 
##'
##' @details Comments about the arguments that have been removed from the code: \itemize{
##' \item \bold{ar}: orginial accrual rate from data from Corine is 0.86 per week, hence we multiply by 2 for by 14 days. As to low, we further multiply by 2
##' \item \bold{cor011,corij1,cor0j1}: from data from Corine
##' \item \bold{MissProb}: miss both V1 and V2 (1,1), miss V1 and but not V2 (1,2), do not miss V1 and but miss V2 (2,1), miss none (2,2)
##' }

## * FCT_simGSD (code)
FCT_simGSD <- function(method, kMax, InfoR.i, InfoR.d, PropForInterim, deltaPlan, allsd, mean0, 
                       deltaTrue = NULL, n.sample = NULL, alpha = 0.025, beta = 0.2, rho_alpha = 2, rho_beta = 2, theDelta.t = 1.50001,
                       path = NULL, name = "simGSD", export.tempo = TRUE, n.sim = 250, seed = NA, trace = TRUE, 
                       rand.block = c(1,1,0,0), ar = (0.86*2)*2,
                       cor011 = -0.15, corij1 = 0.68, cor0j1 = -0.27, TimeFactor = 14,
                       MissProb = structure(c(5/104, 1/104, 6/104, 92/104), dim = c(2L, 2L), dimnames = list(c("V1 missing", "V1 not missing"), c("V2 missing", "V2 not missing")))
                       ){

    require(DelayedGSD)
                                        
    ## ** normalize user input
    if(length(InfoR.i)!=kMax){
        stop("Length of argument \'InfoR.i\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax i.e. ",kMax,".\n",sep="")
    }
    if(InfoR.i[kMax]!=1){
        stop("Argument \'InfoR.i\' at kMax should be 1. \n")
    }
    if(length(InfoR.d)!=kMax){
        stop("Length of argument \'InfoR.d\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax i.e. ",kMax,".\n",sep="")
    }
    if(InfoR.d[kMax]!=1){
        stop("Argument \'InfoR.d\' at kMax should be 1. \n")
    }
    if(length(PropForInterim)!=(kMax-1)){
        stop("Length of argument \'PropForInterim\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax-1 i.e. ",kMax-1,".\n",sep="")
    }
    if(kMax>2 && any(diff(PropForInterim)<0)){
        stop("Argument \'PropForInterim\' should be increasing. \n")
    }
    if(length(deltaPlan)!=(kMax+1)){
        stop("Length of argument \'deltaPlan\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax+1 i.e. ",kMax+1,".\n",sep="")
    }
    if(is.null(deltaTrue)){
        deltaTrue <- deltaPlan
    }else if(length(deltaTrue)!=(kMax+1)){
        stop("Length of argument \'deltaTrue\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax+1 i.e. ",kMax+1,".\n",sep="")
    }
    if(length(allsd)!=(kMax+1)){
        stop("Length of argument \'allsd\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax+1 i.e. ",kMax+1,".\n",sep="")
    }
    if(length(mean0)!=(kMax+1)){
        stop("Length of argument \'mean0\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax+1 i.e. ",kMax+1,".\n",sep="")
    }
    if(length(n.sim)==1){
        n.run <- 1
        multirun <- FALSE
    }else if(length(n.sim)==2){
        n.run <- n.sim[2]
        n.sim <- n.sim[1]
        multirun <- TRUE
    }else{
        stop("Argument \'n.sim\' must have length 1 or 2. \n")
    }
    vec.sim <- ((n.run-1)*n.sim + 1):(n.run*n.sim) ## indices of all iterations for this replicate
    
    n.method <- length(method)
    grid.method <- method2num[method,,drop=FALSE]

    ## ** Sample size
    deltaPower <- abs(deltaPlan[kMax+1]) # effect (NOT Z-scale/unit, but outcome scale/unit!) that the study is powered for: should we choose ourselves or compute from other numbers above ???
    if(is.null(n.sample)){ ## deltaPlan
        n.sample <- ceiling(2*2*((allsd[kMax+1]/deltaPower)^2)*(qnorm(1-beta)-qnorm(alpha))^2) #104 with Corine's data # should we choose ourselves or compute from the above numbers ???
    }
    
    ## ** Seed
    if(!is.na(seed)){
        if(n.sim*n.run>1e5){
            stop("Could not generate the seeds, argument \'n\' in function sample.int need to be changed internally. \n")
        }
        set.seed(seed)
        allseeds <- sample.int(n = 1e5, size = n.sim*n.run, replace=FALSE) #x=1:(.Machine$integer.max) seems to be maximal possible
    }

    ## ** Welcome message
    if(trace){
        if(multirun){
            cat("\t\tSimulation study for Group Sequential trial number ",n.run,". \n",
                "\t\t(",n.sim," iterations ",vec.sim[1],":",vec.sim[n.sim],", sample size: ",n.sample,", seed: ",seed,"). \n\n",sep="")
        }else{
            cat("\t\tSimulation study for Group Sequential trial. \n",
                "\t\t(",n.sim," iterations, sample size: ",n.sample,", seed: ",seed,"). \n\n",sep="")
        }
    }

    ## ** Compute inflation factor and sample size
    plannedB <- vector(mode = "list", length = n.method)
    for(iMeth in 1:n.method){ ## iMeth <- 1

        plannedB[[iMeth]] <- CalcBoundaries(kMax = kMax,  
                                            alpha = alpha, 
                                            beta = beta,  
                                            InfoR.i = InfoR.i,  
                                            InfoR.d = InfoR.d,  
                                            rho_alpha = rho_alpha,  
                                            rho_beta = rho_beta,  
                                            method = method2num[iMeth,"method"],  
                                            cNotBelowFixedc = method2num[method[iMeth],"fixC"],
                                            bindingFutility = method2num[method[iMeth],"binding"],
                                            PowerCorrection = method2num[method[iMeth],"correction"],
                                            delta = utils::tail(deltaPlan,1))
        ## summary(plannedB[[1]])
        ## coef(plannedB[[iMeth]], type = "information")
    }
    inflationFactor <- unlist(lapply(plannedB,function(iP){iP$planned$InflationFactor}))
    nGSD <- ceiling(n.sample*inflationFactor)

    ## ** Loop
    RES <- NULL # initialize results to save
    for(iSim in 1:n.sim){ ## iSim <- 1
        startComp <- Sys.time()
        myseedi <- allseeds[vec.sim[iSim]]

        if(trace){
            if(multirun){
                cat("simulation ",iSim,"(",vec.sim[iSim],")/",n.sim,": seed ",myseedi,"\n",sep="")
            }else{
                cat("simulation ",iSim,"/",n.sim,": seed ",myseedi,"\n",sep="")
            }
        }

        ## *** Generate data
        resData <- GenData(n=n.sample, 
                           N.fw=kMax,
                           rand.block=rand.block,
                           allsd=allsd,
                           mean0=mean0,
                           delta=deltaTrue, 
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
        data <- resData$d
   
        ## *** interim
        ## data
        ## Here we stop inclusion data collection for the interim analysis as soon as
        ## half of the participants have completed (or had the opportunity to complete) the follow-up
        t.interim <- data$t3[ceiling(n.sample*PropForInterim)]
        data.interim <- SelectData(data, t=t.interim)

        ## analysis
        lmm.interim <- analyzeData(data.interim, ddf = "nlme", data.decision = sum(data$t1 <= t.interim + theDelta.t*TimeFactor), getinfo = TRUE, trace = TRUE)

        ## update boundary
        currentGSD <- vector(mode = "list", length = n.method)
        out.interim <- vector(mode = "list", length = n.method)
        for(iMeth in 1:n.method){ ## iMeth <- 1

            currentGSD[[iMeth]] <- update(plannedB[[iMeth]], delta = lmm.interim, trace = FALSE)

            iConfint.interim <- confint(currentGSD[[iMeth]])
            iInfo.interim <- coef(currentGSD[[iMeth]], type = "information")
            iBoundary.interim <- coef(currentGSD[[iMeth]], type = "boundary")
            iDecision.interim <- coef(currentGSD[[iMeth]], type = "decision")

            out.interim[[iMeth]] <-  data.frame(statistic = iConfint.interim[1,"statistic"],
                                                estimate_ML = iConfint.interim[1,"estimate"],
                                                se_ML = iConfint.interim[1,"se"],
                                                n = utils::tail(stats::na.omit(currentGSD[[iMeth]]$n.obs),1),
                                                info = iInfo.interim[1,"Interim"],
                                                infoPC = iInfo.interim[1,"Interim.pc"],
                                                info.pred = iInfo.interim[1,"Decision"],
                                                infoPC.pred = iInfo.interim[1,"Decision.pc"],
                                                uk = iBoundary.interim[1,"Ebound"],
                                                lk = iBoundary.interim[1,"Fbound"],
                                                decision = iDecision.interim["decision","stage 1"],
                                                reason = iDecision.interim["reason.interim","stage 1"],
                                                I.decreasing = sum(iDecision.interim["reason.interim",]=="decreasing information", na.rm = TRUE),
                                                Imax.reached = sum(iDecision.interim["reason.interim",]=="Imax reached", na.rm = TRUE))
        }
        ## currentGSD[[1]]
        ## plot(currentGSD[[1]])

        ## *** decision
        data.decision <- data[which(data$t1 <= t.interim + theDelta.t*TimeFactor),]
        lmm.decision <- analyzeData(data.decision, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    
        out.decision <- vector(mode = "list", length = n.method)
        for(iMeth in 1:n.method){ ## iMeth <- 1
          
            if(out.interim[[iMeth]]$decision == "stop"){
                currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmm.decision, trace = FALSE)
                ## plot(currentGSD[[iMeth]])

                iConfint.decision <- confint(currentGSD[[iMeth]], method = c("ML","MUE"))
                iInfo.decision <- coef(currentGSD[[iMeth]], type = "information")
                iBoundary.decision <- coef(currentGSD[[iMeth]], type = "boundary")
                iDecision.decision  <- coef(currentGSD[[iMeth]], type = "decision")

                out.decision[[iMeth]] <- data.frame(statistic = iConfint.decision[1,"statistic"],
                                                    p.value_ML = iConfint.decision[iConfint.decision$method == "ML","p.value"],
                                                    lower_ML = iConfint.decision[iConfint.decision$method == "ML","lower"],
                                                    upper_ML = iConfint.decision[iConfint.decision$method == "ML","upper"],
                                                    estimate_ML = iConfint.decision[iConfint.decision$method == "ML","estimate"],
                                                    p.value_MUE = iConfint.decision[iConfint.decision$method == "MUE","p.value"],
                                                    lower_MUE = iConfint.decision[iConfint.decision$method == "MUE","lower"],
                                                    upper_MUE = iConfint.decision[iConfint.decision$method == "MUE","upper"],
                                                    estimate_MUE = iConfint.decision[iConfint.decision$method == "MUE","estimate"],
                                                    n = utils::tail(stats::na.omit(currentGSD[[iMeth]]$n.obs),1),
                                                    info = iInfo.decision[1,"Interim"],
                                                    infoPC = iInfo.decision[1,"Interim.pc"],
                                                    ck = iBoundary.decision[1,"Cbound"],
                                                    decision = unname(iDecision.decision["decision","stage 2"]),
                                                    I.decreasing = sum(currentGSD[[iMeth]]$conclusion=="decreasing information", na.rm = TRUE),
                                                    Imax.reached = sum(currentGSD[[iMeth]]$conclusion=="Imax reached", na.rm = TRUE)
                                                    )

            }else{
                ## update information
                currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmm.decision, k = 1, type.k = "decision", trace = FALSE)

                iInfo.decision <- coef(currentGSD[[iMeth]], type = "information")
                iBoundary.decision <- coef(currentGSD[[iMeth]], type = "boundary")

                out.decision[[iMeth]] <- data.frame(statistic = NA,
                                                    p.value_ML = NA,
                                                    lower_ML = NA,
                                                    upper_ML = NA,
                                                    estimate_ML = NA,
                                                    p.value_MUE = NA,
                                                    lower_MUE = NA,
                                                    upper_MUE = NA,
                                                    estimate_MUE = NA,
                                                    n = NA,
                                                    info = iInfo.decision[1,"Decision"],
                                                    infoPC = iInfo.decision[1,"Decision.pc"],
                                                    ck = iBoundary.decision[1,"Cbound"],
                                                    decision = NA,
                                                    I.decreasing = NA,
                                                    Imax.reached = NA)
            }
        }
                                        # }}}
                                        # {{{ Analyze data at decision

        ## *** finale
        data.final <- data
        lmm.final <- analyzeData(data.final, ddf = "nlme", getinfo = TRUE, trace = TRUE)

        out.final <- vector(mode = "list", length = n.method)
        for(iMeth in 1:n.method){ ## iMeth <- 1
            if(out.interim[[iMeth]]$decision == "stop"){
                out.final[[iMeth]] <- data.frame(statistic = NA,
                                                 p.value_ML = NA,
                                                 lower_ML = NA,
                                                 upper_ML = NA,
                                                 estimate_ML = NA,
                                                 p.value_MUE = NA,
                                                 lower_MUE = NA,
                                                 upper_MUE = NA,
                                                 estimate_MUE = NA,
                                                 n = NA,
                                                 info = NA,
                                                 infoPC = NA,
                                                 ck = NA,
                                                 decision = NA,
                                                 I.decreasing = NA,
                                                 Imax.reached = NA)

            }else{
                currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmm.final, trace = FALSE)

                ## plot(test)
                ## summary(test)
                iConfint.final <- confint(currentGSD[[iMeth]], method = c("ML","MUE"))
                iInfo.final <- coef(currentGSD[[iMeth]], type = "information")
                iBoundary.final <- coef(currentGSD[[iMeth]], type = "boundary")
                iDecision.final  <- coef(currentGSD[[iMeth]], type = "decision")

                out.final[[iMeth]] <- data.frame(statistic = iConfint.final[1,"statistic"],
                                                 p.value_ML = iConfint.final[iConfint.final$method == "ML","p.value"],
                                                 lower_ML = iConfint.final[iConfint.final$method == "ML","lower"],
                                                 upper_ML = iConfint.final[iConfint.final$method == "ML","upper"],
                                                 estimate_ML = iConfint.final[iConfint.final$method == "ML","estimate"],
                                                 p.value_MUE = iConfint.final[iConfint.final$method == "MUE","p.value"],
                                                 lower_MUE = iConfint.final[iConfint.final$method == "MUE","lower"],
                                                 upper_MUE = iConfint.final[iConfint.final$method == "MUE","upper"],
                                                 estimate_MUE = iConfint.final[iConfint.final$method == "MUE","estimate"],
                                                 n = utils::tail(stats::na.omit(currentGSD[[iMeth]]$n.obs),1),
                                                 info = iInfo.final[1,"Interim"],
                                                 infoPC = iInfo.final[1,"Interim.pc"],
                                                 ck = iBoundary.final[1,"Cbound"],
                                                 decision = unname(coef(currentGSD[[iMeth]], type = "decision")["decision","stage 2"]),
                                                 I.decreasing = sum(currentGSD[[iMeth]]$conclusion=="decreasing information", na.rm = TRUE),
                                                 Imax.reached = sum(currentGSD[[iMeth]]$conclusion=="Imax reached", na.rm = TRUE)
                                                 )
            }
        }
                                        # }}}

        stopComp <- Sys.time()
                                        # {{{ Save results
        outMerge <- do.call(rbind,lapply(1:n.method, function(iMeth){ ## iMeth <- 
            iNames <- unique(c(names(out.interim[[iMeth]]),names(out.decision[[iMeth]]),names(out.final[[iMeth]])))
            iMerge <- data.frame(matrix(NA, ncol = length(iNames)+3, nrow = 3, dimnames = list(NULL, c("method", "stage", "type", iNames))))
            iMerge[1,c("method","stage","type",names(out.interim[[iMeth]]))] <- data.frame(method = iMeth, stage = 1, type = "interim", out.interim[[iMeth]]) 
            iMerge[2,c("method","stage","type",names(out.decision[[iMeth]]))] <- data.frame(method = iMeth, stage = 1, type = "decision", out.decision[[iMeth]]) 
            iMerge[3,c("method","stage","type",names(out.final[[iMeth]]))] <- data.frame(method = iMeth, stage = 2, type = "final", out.final[[iMeth]])
            return(iMerge)
        }))

    ## outMerge[outMerge$method==3,]

        out <- cbind(
            ## results
            outMerge,
            ## simulation details
            time.interim = t.interim,
            seed=myseedi,             
            nX1.interim = sum(!is.na(data.interim$X1)),
            nX2.interim = sum(!is.na(data.interim$X2)),
            nX3.interim = sum(!is.na(data.interim$X3)),
            ## computation time
            computation.time=as.double(round(difftime(stopComp,startComp,units="secs"),3))
        )
        ## names(out) <- myColNames
        RES <- rbind(RES,out)
        save(RES,file=file.path(path.res,paste0(name,"(tempo)-",iter_sim,".rda")))
                                        # }}}
    }

    ## ** final export
    rownames(RES) <- NULL
    save(RES,file=file.path(path,paste0(name,"-",iter_sim,".rda")))

    return(invisible(RES))
}

## * method2num
method2num <- rbind(expand.grid(method = 1,
                                binding = c(TRUE,FALSE),
                                correction = c(TRUE,FALSE),
                                fixC = c(TRUE,FALSE)),
                    expand.grid(method = 2,
                                binding = c(TRUE,FALSE),
                                correction = FALSE,
                                fixC = c(TRUE,FALSE)),
                    expand.grid(method = 3,
                                binding = c(TRUE,FALSE),
                                correction = FALSE,
                                fixC = TRUE))
method2num <- data.frame(index = 1:NROW(method2num), method2num)
##    index method binding correction  fixC
## 1      1      1    TRUE       TRUE  TRUE
## 2      2      1   FALSE       TRUE  TRUE
## 3      3      1    TRUE      FALSE  TRUE
## 4      4      1   FALSE      FALSE  TRUE
## 5      5      1    TRUE       TRUE FALSE
## 6      6      1   FALSE       TRUE FALSE
## 7      7      1    TRUE      FALSE FALSE
## 8      8      1   FALSE      FALSE FALSE
## 9      9      2    TRUE      FALSE  TRUE
## 10    10      2   FALSE      FALSE  TRUE
## 11    11      2    TRUE      FALSE FALSE
## 12    12      2   FALSE      FALSE FALSE
## 13    13      3    TRUE      FALSE  TRUE
## 14    14      3   FALSE      FALSE  TRUE

##----------------------------------------------------------------------
### FCT.r ends here
