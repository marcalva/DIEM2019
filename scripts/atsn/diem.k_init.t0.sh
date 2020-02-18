#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 4
#$ -l h_data=8000M,h_rt=4:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -t 1-11
#$ -V
#$ -o diem.k_init.t0.sh.log.$TASK_ID

. /u/local/Modules/default/init/modules.sh
module load R/3.5.1
module load hdf5

export HOME=/u/project/pajukant/malvarez/
export R_LIBS_USER=/u/project/pajukant/malvarez/lib/R_%V

Rscript diem.k_init.t0.R $SGE_TASK_ID

