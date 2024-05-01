### FCT.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  1 2022 (15:45) 
## Version: 
## Last-Updated: apr 29 2024 (09:40) 
##           By: Brice Ozenne
##     Update #: 131
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
