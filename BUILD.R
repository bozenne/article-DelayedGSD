### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2022 (16:40) 
## Version: 
## Last-Updated: okt 25 2023 (19:33) 
##           By: Brice Ozenne
##     Update #: 81
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

## cd /projects/biostat01/people/hpl802/DelayedGSD/
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "x:/DelayedGSD"
}else if(system("whoami",intern=TRUE)=="hpl802"){  
    path <- "."
}
path.results <- file.path(path,"Results")
path.Bresults <- file.path(path,"Results-built")
export <- TRUE

## * function used to collect the results from different files
loadRes <- function(path, tempo.file = FALSE, type = NULL,
                    export.attribute = NULL, trace = 2, space = "     "){
    all.files <- list.files(path)
    file.tempo <- grep("(tempo)",all.files,value = TRUE)
    file.final <- setdiff(all.files, file.tempo)

    if(tempo.file){
        file.read <- file.tempo
    }else{
        file.read <- file.final
    }
    if(!is.null(type)){
        file.read <- grep(pattern=type,x=file.read,value=TRUE)
    }

    n.file <- length(file.read)
    if(n.file==0){
        if(trace>1){
            cat(space,"No file found in ",path,". \n\n",sep="")
        }
        return(NULL)
    }
    myApply <- switch(as.character(as.logical(trace)),
                      "TRUE" = pbapply::pblapply,
                      "FALSE" = lapply)

    ls.out <- do.call(myApply, args = list(X = 1:n.file, FUN = function(iFile){
        if(grepl("\\.rds$",file.read[iFile])){
            iRead <- try(readRDS(file = file.path(path,file.read[iFile])))
        }else if(grepl("\\.rda$", file.read[iFile])){
            iRead <- try(load(file = file.path(path,file.read[iFile])))
            if(!inherits(iRead,"try-error")){
                iRead <- eval(parse(text = iRead))
            }
        }else{
            return(NULL)
        }
        if(inherits(iRead,"try-error")){
            return(NULL)
        }else{
            iOut <- cbind(data.table::as.data.table(iRead),
                          file = file.read[iFile])
            return(iOut)
        }
    }))

    out <- do.call(rbind, ls.out)

    Urow <- unique(sapply(ls.out, NROW))
    if(length(Urow)==1){
        cat(space,length(ls.out)," files with ",Urow," lines \n\n",sep="")
    }else{
        cat(space,length(ls.out)," files, with from ",min(Urow)," to ",max(Urow)," lines \n\n",sep="")
    }
    return(out)
}

## * Process results (2 stages)
## ** Aggregate files and export 
dir.2stage <- grep("2stage",list.dirs(path = path.results), value = TRUE)
name.2stage <- stats::setNames(paste0("res2stage_",sapply(strsplit(dir.2stage, split = "2stage_"),"[",2)), dir.2stage)

for(iDir in dir.2stage){ ## iDir <- dir.2stage[1]
    iName <- name.2stage[iDir]
    cat(which(iDir == dir.2stage),") read \'",iDir,"\' \n",
        "     save in \'",iName,"\' \n", sep = "")
    assign(x = iName, value = loadRes(iDir))
    if(export){
        saveRDS(eval(parse(text=iName)), file = file.path(path.results,paste0(iName,".rds") ))
    }
}

## dt <- loadRes("x:/DelayedGSD/Results/2stage_missing_binding_ar10_power", tempo.file = TRUE)
## dt[, .N, by = "file"]
## length(unique(dt$file))
## dt[file=="sim-2stage_missing_binding_ar10_power-1(tempo)_100.rds"]


## length(unique(res2stage_missing_binding_ar10_power$file)) ## 99
## length(unique(res2stage_missing_binding_ar10_typeI$file))        ## 99
## length(unique(res2stage_missing_binding_ar5_power$file))         
## length(unique(res2stage_missing_binding_ar5_typeI$file))      
## length(unique(res2stage_missing_fixC_binding_ar10_power$file)) ## 98
## length(unique(res2stage_missing_fixC_binding_ar10_typeI$file))   ## 98
## length(unique(res2stage_missing_fixC_binding_ar5_power$file))    
## length(unique(res2stage_missing_fixC_binding_ar5_typeI$file))    
## length(unique(res2stage_missing_fixC_nonbinding_ar10_power$file)) ## 99
## length(unique(res2stage_missing_fixC_nonbinding_ar10_typeI$file)) ## 99
## length(unique(res2stage_missing_fixC_nonbinding_ar5_power$file)) 
## length(unique(res2stage_missing_fixC_nonbinding_ar5_typeI$file)) 
## length(unique(res2stage_missing_nonbinding_ar10_power$file))     
## length(unique(res2stage_missing_nonbinding_ar10_typeI$file))     
## length(unique(res2stage_missing_nonbinding_ar5_power$file))      
## length(unique(res2stage_missing_nonbinding_ar5_typeI$file))      
## length(unique(res2stage_nomissing_binding_ar5_power$file))       
## length(unique(res2stage_nomissing_binding_ar5_typeI$file))       

## ** Aggregate scenario and export 
legend.2stage <- data.frame(name = name.2stage,
                            scenario = 1:length(name.2stage),
                            missing = grepl("nomissing", name.2stage)== FALSE,
                            binding = grepl("nonbinding", name.2stage)== FALSE,
                            fixC = grepl("fixC", name.2stage),
                            ar = gsub("ar","",sapply(strsplit(name.2stage, split = "_"), function(iVec){tail(iVec,2)[1]})),
                            hypo = sapply(strsplit(name.2stage, split = "_"), function(iVec){tail(iVec,1)})
                            )

res2stage <- do.call(rbind,lapply(name.2stage, function(iName){ ## iName <- name.Bresults[1]
    data.table(legend.2stage[name.2stage==iName,-1], eval(parse(text = iName)))
}))
res2stage$computation.time <- NULL
##res2stage$file <- NULL
res2stage$sigma <- NULL
res2stage$ar <- as.numeric(res2stage$ar)

if(export){
    saveRDS(res2stage, file = file.path(path.Bresults,"res2stage.rds"))
}

## res2stage <- readRDS(file = file.path(path.Bresults,"res2stage.rds"))
## unique(res2stage[,.N, by = c("scenario","method","type")]$N)
## [1] 10000

## * Process results (3 stages)
## ** Aggregate files and export 
dir.3stage <- grep("3stage",list.dirs(path = path.results), value = TRUE)
name.3stage <- stats::setNames(paste0("res3stage_",sapply(strsplit(dir.3stage, split = "3stage_"),"[",2)), dir.3stage)

for(iDir in dir.3stage){ ## iDir <- dir.3stage[1]
    iName <- name.3stage[iDir]
    cat(which(iDir == dir.3stage),") read \'",iDir,"\' \n",
        "     save in \'",iName,"\' \n", sep = "")
    assign(x = iName, value = loadRes(iDir, tempo.file = TRUE))
    if(export){
        saveRDS(eval(parse(text=iName)), file = file.path(path.results,paste0(iName,".rds") ))
    }
}

## ** Aggregate scenario and export 
legend.3stage <- data.frame(name = name.3stage,
                            scenario = 1:length(name.3stage),
                            missing = grepl("nomissing", name.3stage)== FALSE,
                            binding = grepl("nonbinding", name.3stage)== FALSE,
                            fixC = grepl("fixC", name.3stage),
                            ar = gsub("ar","",sapply(strsplit(name.3stage, split = "_"), function(iVec){tail(iVec,2)[1]})),
                            hypo = sapply(strsplit(name.3stage, split = "_"), function(iVec){tail(iVec,1)})
                            )

res3stage <- do.call(rbind,lapply(name.3stage, function(iName){ ## iName <- name.Bresults[2]
    iValue <- eval(parse(text = iName))
    if(is.null(iValue)){
        return(NULL)
    }else{
        return(data.table(legend.3stage[name.3stage==iName,-1], iValue))
    }
}))
res3stage$computation.time <- NULL
##res3stage$file <- NULL
res3stage$sigma <- NULL
res3stage$ar <- as.numeric(res3stage$ar)

if(export){
    saveRDS(res3stage, file = file.path(path.Bresults,"res3stage.rds"))
}

## res3stage <- readRDS(file = file.path(path.Bresults,"res3stage.rds"))
## unique(res3stage[,.N, by = c("scenario","method","type")]$N)
## [1] 10000



##----------------------------------------------------------------------
### BUILD.R ends here
