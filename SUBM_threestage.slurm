#!/bin/bash

#SBATCH --job-name=3stageGST  # Job name
#SBATCH --array=1-100     # Number of jobs
#SBATCH --partition=long    # Queue. Use long for run time >8h and standard otherwise
#SBATCH --mem=1G # Memory limit, e.g. reserve 1 GB memory 
#SBATCH --output=output/slurm/slurm-%x-%a-%J.out # Direct output to subdirectory
#SBATCH --error=output/slurm/slurm-%x-%a-%J.out # Same file for output and error

## no missing value
R CMD BATCH --vanilla '--args missing=FALSE binding=TRUE cNotBelowFixedc=FALSE ar.factor=5 delta.factor=0.6 n.method=3' BATCH_threestage.R output/3stage_nomissing_binding_ar5_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=FALSE binding=TRUE cNotBelowFixedc=FALSE ar.factor=5 delta.factor=0 n.method=3' BATCH_threestage.R output/3stage_nomissing_binding_ar5_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
## missing values:  (non-binding,binding) ar(5,10) (power,type1)
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=FALSE ar.factor=5 delta.factor=0.6 n.method=3' BATCH_threestage.R output/3stage_missing_binding_ar5_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=FALSE ar.factor=5 delta.factor=0 n.method=3' BATCH_threestage.R output/3stage_missing_binding_ar5_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=FALSE ar.factor=10 delta.factor=0.6 n.method=3' BATCH_threestage.R output/3stage_missing_binding_ar10_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=FALSE ar.factor=10 delta.factor=0 n.method=3' BATCH_threestage.R output/3stage_missing_binding_ar10_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=5 delta.factor=0.6 n.method=3' BATCH_threestage.R output/3stage_missing_nonbinding_ar5_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=5 delta.factor=0 n.method=3' BATCH_threestage.R output/3stage_missing_nonbinding_ar5_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=10 delta.factor=0.6 n.method=3' BATCH_threestage.R output/3stage_missing_nonbinding_ar10_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=FALSE ar.factor=10 delta.factor=0 n.method=3' BATCH_threestage.R output/3stage_missing_nonbinding_ar10_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
## missing values: fix C (non-binding,binding) ar(5,10) (power,type1)
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=TRUE ar.factor=5 delta.factor=0.6 n.method=3' BATCH_threestage.R output/3stage_missing_fixC_binding_ar5_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=TRUE ar.factor=5 delta.factor=0 n.method=3' BATCH_threestage.R output/3stage_missing_fixC_binding_ar5_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=TRUE ar.factor=10 delta.factor=0.6 n.method=3' BATCH_threestage.R output/3stage_missing_fixC_binding_ar10_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=TRUE cNotBelowFixedc=TRUE ar.factor=10 delta.factor=0 n.method=3' BATCH_threestage.R output/3stage_missing_fixC_binding_ar10_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=TRUE ar.factor=5 delta.factor=0.6 n.method=3' BATCH_threestage.R output/3stage_missing_fixC_nonbinding_ar5_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=TRUE ar.factor=5 delta.factor=0 n.method=3' BATCH_threestage.R output/3stage_missing_fixC_nonbinding_ar5_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=TRUE ar.factor=10 delta.factor=0.6 n.method=3' BATCH_threestage.R output/3stage_missing_fixC_nonbinding_ar10_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla '--args missing=TRUE binding=FALSE cNotBelowFixedc=TRUE ar.factor=10 delta.factor=0 n.method=3' BATCH_threestage.R output/3stage_missing_fixC_nonbinding_ar10_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout

## go to directory    ## cd c:/Users/hpl802/Documents/Github/DelayedGSD/Simulations/brice-biostat-cluster/
                      ## cd /projects/biostat01/people/hpl802/DelayedGSD/

## clean outputs      ## rm -r output/slurm/*;
                      ## rm -r output/3stage_nomissing_binding_ar5_power/*;rm -r output/3stage_nomissing_binding_ar5_typeI/*;
                      ## rm -r output/3stage_missing_binding_ar5_power/*;rm -r output/3stage_missing_binding_ar5_typeI/*;
                      ## rm -r output/3stage_missing_binding_ar10_power/*;rm -r output/3stage_missing_binding_ar10_typeI/*;
                      ## rm -r output/3stage_missing_nonbinding_ar5_power/*;rm -r output/3stage_missing_nonbinding_ar5_typeI/*;
                      ## rm -r output/3stage_missing_nonbinding_ar10_power/*;rm -r output/3stage_missing_nonbinding_ar10_typeI/*;
                      ## rm -r output/3stage_missing_fixC_binding_ar5_power/*;rm -r output/3stage_missing_fixC_binding_ar5_typeI/*;
                      ## rm -r output/3stage_missing_fixC_binding_ar10_power/*;rm -r output/3stage_missing_fixC_binding_ar10_typeI/*;
                      ## rm -r output/3stage_missing_fixC_nonbinding_ar5_power/*;rm -r output/3stage_missing_fixC_nonbinding_ar5_typeI/*;
                      ## rm -r output/3stage_missing_fixC_nonbinding_ar10_power/*;rm -r output/3stage_missing_fixC_nonbinding_ar10_typeI/*;


## clean results      ## 
                      ## rm -r Results/3stage_nomissing_binding_ar5_power/*;rm -r Results/3stage_nomissing_binding_ar5_typeI/*;
                      ## rm -r Results/3stage_missing_binding_ar5_power/*;rm -r Results/3stage_missing_binding_ar5_typeI/*;
                      ## rm -r Results/3stage_missing_binding_ar10_power/*;rm -r Results/3stage_missing_binding_ar10_typeI/*;
                      ## rm -r Results/3stage_missing_nonbinding_ar5_power/*;rm -r Results/3stage_missing_nonbinding_ar5_typeI/*;
                      ## rm -r Results/3stage_missing_nonbinding_ar10_power/*;rm -r Results/3stage_missing_nonbinding_ar10_typeI/*;
                      ## rm -r Results/3stage_missing_fixC_binding_ar5_power/*;rm -r Results/3stage_missing_fixC_binding_ar5_typeI/*;
                      ## rm -r Results/3stage_missing_fixC_binding_ar10_power/*;rm -r Results/3stage_missing_fixC_binding_ar10_typeI/*;
                      ## rm -r Results/3stage_missing_fixC_nonbinding_ar5_power/*;rm -r Results/3stage_missing_fixC_nonbinding_ar5_typeI/*;
                      ## rm -r Results/3stage_missing_fixC_nonbinding_ar10_power/*;rm -r Results/3stage_missing_fixC_nonbinding_ar10_typeI/*;

## submission command ## sbatch SUBM_threestage.slurm
## submission output  ## Submitted batch job 6056333 (time using C-c t: tirsdag 23-10-31 at 17:22) 

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
