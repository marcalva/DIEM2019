#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4000M,h_rt=16:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -V
#$ -o diem.integrated.fltr.1.sh.log

. /u/local/Modules/default/init/modules.sh
module load R/3.5.1
module load hdf5

export HOME=/u/project/pajukant/malvarez/
export R_LIBS_USER=/u/project/pajukant/malvarez/lib/R_%V

Rscript diem.integrated.fltr.1.R
