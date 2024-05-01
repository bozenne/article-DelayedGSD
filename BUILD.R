### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2022 (16:40) 
## Version: 
## Last-Updated: apr 30 2024 (18:05) 
##           By: Brice Ozenne
##     Update #: 94
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

readSampleSize <- function(n.stage){
    out.dir <- grep(paste0(n.stage,"stage"),list.files("output"), value = TRUE)
    out.file <- sapply(out.dir, function(iDir){list.files(file.path("output",iDir),full.names = TRUE)[1]})
    out.nchar <- sapply(out.file, function(iFile){gsub("Sample size:| ","",grep(readLines(iFile), pattern = "^Sample size", value = TRUE))})
    out.n <-  apply(do.call(rbind,strsplit(out.nchar,split =",",fixed=TRUE)),2,as.numeric)
    rownames(out.n) <- names(out.nchar)
    return(out.n)
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


## quantile(readSampleSize(2)[,1])
##  0%  25%  50%  75% 100% 
## 656  734  734  741  741 


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
                            ar = gsub("^.*_ar([0-9]+)_.*$","\\1",name.2stage),
                            hypo = c("power","typeI")[grepl("_power", name.2stage)+2*grepl("_typeI", name.2stage)],
                            infoBias = -grepl("missinfo0x5",name.2stage)+grepl("missinfo1x5",name.2stage)
                            )

res2stage <- do.call(rbind,lapply(name.2stage, function(iName){ ## iName <- name.Bresults[1]
    iValue <- eval(parse(text = iName))
    if(!is.null(iValue)){
        return(data.table(legend.2stage[name.2stage==iName,-1], iValue))
    }else{
        return(NULL)
    }
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

## ** find missing seed
## test.scenario <- res2stage[, .(seed= length(unique(seed)),dir=file[1]), by = "scenario"]
## test.scenario[seed!=max(seed),]
## 
## Useed <- res2stage[, unique(seed)]
## NAseed <- res2stage[scenario == 9, setdiff(Useed,unique(seed))]
## NAseed[1]
## [1] 506346813
## pbfile <- unique(sapply(strsplit(res2stage[seed %in% NAseed & scenario == 1, file],split = "-"), "[",3))
## [1] "33_100.rds"

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

## quantile(readSampleSize(3)[,1])
 ##  0%  25%  50%  75% 100% 
 ## 671  750  750  762  762 

## ** Aggregate scenario and export 
legend.3stage <- data.frame(name = name.3stage,
                            scenario = 1:length(name.3stage),
                            missing = grepl("nomissing", name.3stage)== FALSE,
                            binding = grepl("nonbinding", name.3stage)== FALSE,
                            fixC = grepl("fixC", name.3stage),
                            ar = gsub("^.*_ar([0-9]+)_.*$","\\1",name.3stage),
                            hypo = c("power","typeI")[grepl("_power", name.3stage)+2*grepl("_typeI", name.3stage)]
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


## * debug
library(DelayedGSD)
## Sample size: 270, 270, 270
## "seed 529169576 for j=9903 (index 3) out of 100"
##  missing binding cNotBelowFixedc ar.factor delta.factor n.method
##     TRUE    TRUE           FALSE         1            1        3


df.method <- data.frame(method = 1:3, binding = TRUE, fixC = c(FALSE,FALSE,TRUE))
 
MyMissProb <- matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538),  nrow = 2, ncol = 2,
                     dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")))
 
args.GenData <- list(rand.block = c(1, 1, 0, 0),
                     allsd = c(2.5, 2.1, 2.4),
                     mean0 = c(10, 0, 0),
                     delta = c(0, 0.5, 1)*1,
                     ar = 15,
                     cor.01.1 = -0.15,
                     cor.ij.1 = 0.68,
                     cor.0j.1 = -0.27,
                     MissProb = MyMissProb,
                     DigitsOutcome = 2,
                     TimeFactor = 42,
                     DigitsTime = 0)

res3stage <- operatingDelayedGSD(n.sim = 1, 
                                 method = df.method,
                                 args.GenData = args.GenData,
                                 kMax = 3, InfoR.i = c(0.40,0.65,1), InfoR.d = c(0.50,0.75, 1), delta = 1,
                                 PropForInterim = c(0.35,0.6), lag = 21,
                                 seed = 529169576)

##----------------------------------------------------------------------
### BUILD.R ends here
