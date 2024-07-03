#!/bin/bash

#SBATCH --job-name=2stageGST  # Job name
#SBATCH --array=1-200     # Number of jobs
#SBATCH --partition=long    # Queue. Use long for run time >8h and standard otherwise
#SBATCH --mem=1G # Memory limit, e.g. reserve 1 GB memory 
#SBATCH --output=output/slurm/slurm-%x-%a-%J.out # Direct output to subdirectory
#SBATCH --error=output/slurm/slurm-%x-%a-%J.out # Same file for output and error

####---- no missing value
R CMD BATCH --vanilla '--args missing=FALSE binding=TRUE cNotBelowFixedc=FALSE ar.factor=1 delta.factor=1 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_nomissing_binding_ar1_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=FALSE binding=TRUE cNotBelowFixedc=FALSE ar.factor=1 delta.factor=0 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_nomissing_binding_ar1_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
####---- missing values:  (non-binding,binding) ar(5,10) (power,type1)
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=FALSE ar.factor=1 delta.factor=1 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_binding_ar1_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=FALSE ar.factor=1 delta.factor=0 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_binding_ar1_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=FALSE ar.factor=2 delta.factor=1 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_binding_ar2_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=FALSE ar.factor=2 delta.factor=0 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_binding_ar2_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=1 delta.factor=1 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_nonbinding_ar1_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=1 delta.factor=0 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_nonbinding_ar1_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=2 delta.factor=1 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_nonbinding_ar2_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=2 delta.factor=0 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_nonbinding_ar2_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
####---- missing values: fix C (non-binding,binding) ar(5,10) (power,type1)
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=TRUE ar.factor=1 delta.factor=1 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_fixC_binding_ar1_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=TRUE ar.factor=1 delta.factor=0 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_fixC_binding_ar1_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=TRUE ar.factor=2 delta.factor=1 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_fixC_binding_ar2_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=TRUE ar.factor=2 delta.factor=0 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_fixC_binding_ar2_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=TRUE ar.factor=1 delta.factor=1 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_fixC_nonbinding_ar1_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=TRUE ar.factor=1 delta.factor=0 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_fixC_nonbinding_ar1_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=TRUE ar.factor=2 delta.factor=1 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_fixC_nonbinding_ar2_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=TRUE ar.factor=2 delta.factor=0 n.method=3 incorrect.information=0' BATCH_twostage.R output/2stage_missing_fixC_nonbinding_ar2_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
####---- incorrect information:  non-binding ar 1 power fixC
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=1 delta.factor=1 n.method=3 incorrect.information=0.5' BATCH_twostage.R output/2stage_missing_nonbinding_ar1_power_missinfo0x5/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=1 delta.factor=1 n.method=3 incorrect.information=1.5' BATCH_twostage.R output/2stage_missing_nonbinding_ar1_power_missinfo1x5/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout

## go to directory    ## cd c:/Users/hpl802/Documents/Github/DelayedGSD/Simulations/brice-biostat-cluster/
                      ## cd /projects/biostat01/people/hpl802/DelayedGSD/

## clean outputs      ## rm -r output/slurm/*;
                      ## rm -r output/2stage_nomissing_binding_ar1_power/*;rm -r output/2stage_nomissing_binding_ar1_typeI/*;
                      ## rm -r output/2stage_missing_binding_ar1_power/*;rm -r output/2stage_missing_binding_ar1_typeI/*;
                      ## rm -r output/2stage_missing_binding_ar2_power/*;rm -r output/2stage_missing_binding_ar2_typeI/*;
                      ## rm -r output/2stage_missing_nonbinding_ar1_power/*;rm -r output/2stage_missing_nonbinding_ar1_typeI/*;
                      ## rm -r output/2stage_missing_nonbinding_ar2_power/*;rm -r output/2stage_missing_nonbinding_ar2_typeI/*;
                      ## rm -r output/2stage_missing_fixC_binding_ar1_power/*;rm -r output/2stage_missing_fixC_binding_ar1_typeI/*;
                      ## rm -r output/2stage_missing_fixC_binding_ar2_power/*;rm -r output/2stage_missing_fixC_binding_ar2_typeI/*;
                      ## rm -r output/2stage_missing_fixC_nonbinding_ar1_power/*;rm -r output/2stage_missing_fixC_nonbinding_ar1_typeI/*;
                      ## rm -r output/2stage_missing_fixC_nonbinding_ar2_power/*;rm -r output/2stage_missing_fixC_nonbinding_ar2_typeI/*;
                      ## rm -r output/2stage_missing_nonbinding_ar1_power_missinfo0x5/*;rm -r output/2stage_missing_nonbinding_ar1_power_missinfo1x5/*;


## clean results      ## 
                      ## rm -r Results/2stage_nomissing_binding_ar1_power/*;rm -r Results/2stage_nomissing_binding_ar1_typeI/*;
                      ## rm -r Results/2stage_missing_binding_ar1_power/*;rm -r Results/2stage_missing_binding_ar1_typeI/*;
                      ## rm -r Results/2stage_missing_binding_ar2_power/*;rm -r Results/2stage_missing_binding_ar2_typeI/*;
                      ## rm -r Results/2stage_missing_nonbinding_ar1_power/*;rm -r Results/2stage_missing_nonbinding_ar1_typeI/*;
                      ## rm -r Results/2stage_missing_nonbinding_ar2_power/*;rm -r Results/2stage_missing_nonbinding_ar2_typeI/*;
                      ## rm -r Results/2stage_missing_fixC_binding_ar1_power/*;rm -r Results/2stage_missing_fixC_binding_ar1_typeI/*;
                      ## rm -r Results/2stage_missing_fixC_binding_ar2_power/*;rm -r Results/2stage_missing_fixC_binding_ar2_typeI/*;
                      ## rm -r Results/2stage_missing_fixC_nonbinding_ar1_power/*;rm -r Results/2stage_missing_fixC_nonbinding_ar1_typeI/*;
                      ## rm -r Results/2stage_missing_fixC_nonbinding_ar2_power/*;rm -r Results/2stage_missing_fixC_nonbinding_ar2_typeI/*;
                      ## rm -r Results/2stage_missing_nonbinding_ar1_power_missinfo0x5/*;rm -r Results/2stage_missing_nonbinding_ar1_power_missinfo1x5/*;

## submission command ## sbatch SUBM_twostage.slurm
## submission output  ## Submitted batch job 6378441 (time using C-c t: onsdag 24-06-05 at 15:35) 

## init               ## module module load gcc/11.2.0 R/4.1.2
                      ## cat .bash_profile (in the root directory on the server to visualize the current configuration)
                      ## vi .bash_profile (in the root directory on the server to modify the current configuration)

## commands           ## squeue           : view current jobs 
                      ## squeue -u id     : view current jobs for user id
                      ## scancel          : delete the job with job id 1034 from the queue type
                      ## sinfo            : view state of the queues
                      ## sshare           : check own usage of the server
                        
## documentation      ## from SLURM: https://slurm.schedmd.com/pdfs/summary.pdf: key functions
                      ##           : https://slurm.schedmd.com/sbatch.html (filename pattern)
                      ## from KU   : https://hpc.ku.dk/documentation/slurm.html
                      ##           : https://github.com/aejensen/computing-biostat-ucph/blob/main/README.md
