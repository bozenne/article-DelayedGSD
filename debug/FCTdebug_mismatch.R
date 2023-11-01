### FCTdebug_mismatch.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  5 2023 (13:23) 
## Version: 
## Last-Updated: okt  5 2023 (18:02) 
##           By: Brice Ozenne
##     Update #: 102
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * calcP_new
calcP_new <- function(delta, object){

    requireNamespace("mvtnorm")

    ## ** extract 
    kMax <- object$kMax
    stage <- object$stage$k
    lk <- object$lk
    uk <- object$uk
    ck.unrestricted <- object$ck.unrestricted
    ck <- object$ck
    Info.i <- object$Info.i
    Info.d <- object$Info.d
    estimate <- object$delta[(object$delta$type!="interim")*(object$delta$method=="ML")*(object$delta$stage==stage)>0,"estimate"]

    bindingFutility <- object$bindingFutility
    conclusion <- object$conclusion
    method <- object$method
    reason.interim <- conclusion["reason.interim",]

    ## ** prepare

    ## statistic
    statistic <- estimate * sqrt(Info.d[stage])
    Fstatistic <- statistic - max(0, statistic - ck.unrestricted[stage]) + max(0, statistic - ck[stage])
    
    ## information
    Info.vec <- c(Info.i,Info.d)
    Info.type <- c(rep("interim",length(Info.i)),
                   rep("decision",length(Info.d)-1),
                   rep("final",1))
    Info.indexInterim <- which(Info.type=="interim")
    Info.indexDecision <- which(Info.type=="decision")
    Info.indexFinal <- which(Info.type=="final")
        
    Info.matrix <- diag(1, length(Info.vec))
    Info.matrix[lower.tri(Info.matrix)] <- sqrt((1/Info.vec) %*% t(Info.vec))[lower.tri(Info.matrix)]
    Info.matrix[upper.tri(Info.matrix)] <- t(Info.matrix)[upper.tri(Info.matrix)]

    ## futility bounds
    if(bindingFutility){
        lk.continue <- lk
    }else{
        lk.continue <- rep(-Inf,length(lk))
    }

    ## ** evaluate p-value at each stage
    ls.pvalue <- lapply(1:stage, function(iStage){

        if(reason.interim[iStage]=="decreasing information"){ ## SPECIAL CASE: decreasing information from on interim analysis to another interim analysis

            iOut <- 0
            attr(iOut,"reason") <- "decreasing information"

        }else{

            ## previous interim analyses
            if(iStage==1){
                iSeq_interimM1 <- NULL
            }else{ ## index of previous interim that are not skipped
                iSeq_interimM1 <- intersect(1:(iStage-1),which(reason.interim[1:(iStage-1)]!="decreasing information"))
            }
            
            if(iStage==kMax || reason.interim[iStage]=="Imax reached"){ ## SPECIAL CASE: no interim (either final analysis or decision after skipping interim)

                iIndex <- c(Info.indexInterim[iSeq_interimM1], c(Info.indexDecision,Info.indexFinal)[iStage])

                    iOut <- mvtnorm::pmvnorm(lower = c(lk.continue[iSeq_interimM1], Fstatistic),  
                                             upper = c(uk[iSeq_interimM1],                 Inf),
                                             mean = delta * sqrt(Info.vec[iIndex]),
                                             sigma = Info.matrix[iIndex,iIndex,drop=FALSE])
        
            }else{ ## NORMAL CASE
                ## (and special case decreasing information between interim and decision analysis?)
                iIndex_interim <- c(Info.indexInterim[iSeq_interimM1], Info.indexInterim[iStage])
                iIndex <- c(iIndex_interim, Info.indexDecision[iStage])

                if(iStage == stage){
                    iStatistic <- Fstatistic

                    ## 1- prob to continue (and therefore have a more extreme result than stopping now for futility)
                    if((iStatistic < ck[stage]) || (method == 3 && reason.interim[stage]=="futility")){
                        iTerm1 <- mvtnorm::pmvnorm(lower = lk.continue[c(iSeq_interimM1,stage)],  
                                                   upper = uk[c(iSeq_interimM1,stage)],
                                                   mean = delta * sqrt(Info.vec[iIndex_interim]),
                                                   sigma = Info.matrix[iIndex_interim,iIndex_interim,drop=FALSE])                        
                    }else{
                        iTerm1 <- 0
                    }
                    
                }else{
                    iStatistic <- ck.unrestricted[iStage]
                    iTerm1 <- 0
                }
                
                ## 2- probability to stop for efficacy and conclude more extreme efficacy
                if(method == 3 && iStage == stage && reason.interim[stage]=="futility"){
                    iTerm2 <- mvtnorm::pmvnorm(lower = c(lk.continue[iSeq_interimM1], uk[iStage], min(statistic, ck.unrestricted[stage])),  
                                               upper = c(uk[iSeq_interimM1],                 Inf,        Inf),
                                               mean = delta * sqrt(Info.vec[iIndex]),
                                               sigma = Info.matrix[iIndex,iIndex,drop=FALSE])
                }else{
                    iTerm2 <- mvtnorm::pmvnorm(lower = c(lk.continue[iSeq_interimM1], uk[iStage], iStatistic),  
                                               upper = c(uk[iSeq_interimM1],                 Inf,        Inf),
                                               mean = delta * sqrt(Info.vec[iIndex]),
                                               sigma = Info.matrix[iIndex,iIndex,drop=FALSE])
                }

                ## 3- probability to stop for futility and conclude more extreme efficacy
                if(method %in% c(1,2) || (iStage == stage && iStatistic < ck[stage])){
                    iTerm3 <- mvtnorm::pmvnorm(lower = c(lk.continue[c(iSeq_interimM1)],       -Inf, iStatistic),   
                                               upper = c(uk[iSeq_interimM1],             lk[iStage],        Inf),
                                               mean = delta * sqrt(Info.vec[iIndex]),
                                               sigma = Info.matrix[iIndex,iIndex,drop=FALSE])
                }else if(method == 3 && iStage == stage && reason.interim[stage]=="futility"){
                    iTerm3 <- mvtnorm::pmvnorm(lower = c(lk.continue[c(iSeq_interimM1)],       -Inf, statistic),  ## no correction (i.e. -ck+ck.unrestricted) 
                                               upper = c(uk[iSeq_interimM1],             lk[iStage],        Inf),
                                               mean = delta * sqrt(Info.vec[iIndex]),
                                               sigma = Info.matrix[iIndex,iIndex,drop=FALSE])
                }else{
                    iTerm3 <- 0
                }

                iOut <- unname(iTerm1 + iTerm2 + iTerm3)
                attr(iOut, "terms") <- unname(c(iTerm1,iTerm2,iTerm3))
            }
        }
        return(iOut)

    })

    ## ** export
    out <- do.call("+",ls.pvalue)
    
    return(out)

}


## calcP_futility <- function(delta, object){

##     ## truncated test statistic     
##     if(object$delta$statistic[2]<object$ck[1] && object$delta$statistic[2]>object$ck.unrestricted[1]){
##         statistic <- object$ck.unrestricted[1]
##     }else{
##         statistic <- object$delta$statistic[2]
##     }

##     if(object$delta$statistic[2]<object$ck[1] || (object$method == 3 && object$conclusion["reason.interim",1]=="futility")){
##         term1 <- mvtnorm::pmvnorm(lower = -Inf,  
##                                   upper = object$uk[1],
##                                   mean = delta * sqrt(object$Info.i),
##                                   sigma = matrix(1,1,1))
##     }else{
##         term1 <- 0
##     }

##     if(object$method == 3 && object$conclusion["reason.interim",1]=="futility"){

##         term2 <- mvtnorm::pmvnorm(lower = c(-Inf        , statistic),
##                                   upper = c(object$lk[1], Inf),
##                                   mean = delta * c(sqrt(object$Info.i), sqrt(object$Info.d[1])),
##                                   sigma = cbind(c(1, sqrt(object$Info.i/object$Info.d[1])), 
##                                                 c(sqrt(object$Info.i/object$Info.d[1]), 1)
##                                                 )
##                                   )
##     }else if(object$delta$statistic[2]>=object$ck[1]){
##         ## continuity correction
##         term2 <- mvtnorm::pmvnorm(lower = c(object$uk[1],object$ck.unrestricted[1]),  
##                                   upper = c(Inf,object$ck[1]),
##                                   mean = delta * c(sqrt(object$Info.i), sqrt(object$Info.d[1])),
##                                   sigma = cbind(c(1, sqrt(object$Info.i/object$Info.d[1])), 
##                                                 c(sqrt(object$Info.i/object$Info.d[1]), 1)
##                                                 )
##                                   )

##         if(object$method %in% 1:2){
##             term2 <- term2 + mvtnorm::pmvnorm(lower = c(-Inf,object$ck.unrestricted[1]),  
##                                               upper = c(object$lk[1],object$ck[1]),
##                                               mean = delta * c(sqrt(object$Info.i), sqrt(object$Info.d[1])),
##                                               sigma = cbind(c(1, sqrt(object$Info.i/object$Info.d[1])), 
##                                                             c(sqrt(object$Info.i/object$Info.d[1]), 1)
##                                                             )
##                                               )
##         }

##     }else{
##         ## (this should only be added in case of futility for method 3, as for method 3 flips from fut to eff
##         ##  are not allowed and a futile result should not be considered more extreme than a positive result)
##         term2 <- mvtnorm::pmvnorm(lower = c(-Inf, statistic),
##                                   upper = c(object$lk[1],Inf),
##                                   mean = delta * c(sqrt(object$Info.i), sqrt(object$Info.d[1])),
##                                   sigma = cbind(c(1, sqrt(object$Info.i/object$Info.d[1])), 
##                                                 c(sqrt(object$Info.i/object$Info.d[1]), 1)
##                                                 )
##                                   )
##     }

##     term3 <- mvtnorm::pmvnorm(lower = c(object$uk[1], statistic),  
##                               upper = c(Inf, Inf),
##                               mean = delta * c(sqrt(object$Info.i), sqrt(object$Info.d[1])),
##                               sigma = cbind(c(1, sqrt(object$Info.i/object$Info.d[1])), 
##                                             c(sqrt(object$Info.i/object$Info.d[1]), 1)
##                                             )
##                               )

##     if(object$method %in% 1:2){
##         term3 <- term3 + mvtnorm::pmvnorm(lower = c(-Inf, statistic),  
##                                           upper = c(object$lk[1], Inf),
##                                           mean = delta * c(sqrt(object$Info.i), sqrt(object$Info.d[1])),
##                                           sigma = cbind(c(1, sqrt(object$Info.i/object$Info.d[1])), 
##                                                         c(sqrt(object$Info.i/object$Info.d[1]), 1)
##                                                         )
##                                           )
##     }

##     ## p-value
##     p.value <- unname(term1 + term2 + term3)
    
##     ## export
##     attr(p.value, "terms") <- unname(c(term1,term2,term3))
##     return(p.value)
## }

##----------------------------------------------------------------------
### FCTdebug_mismatch.R ends here
