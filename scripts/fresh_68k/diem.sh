#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4000M,h_rt=10:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -V
#$ -o diem.sh.log

export OMP_NUM_THREADS=4


. /u/local/Modules/default/init/modules.sh
module load R/3.5.1
module load hdf5
module load python/3.6.1_shared

export HOME=/u/project/pajukant/malvarez/

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/local/compilers/gcc/4.7.2/lib:/u/local/compilers/gcc/4.7.2/lib64:/u/project/pajukant/malvarez/lib/
export PATH=$PATH:/u/project/pajukant/malvarez/bin:/u/project/pajukant/malvarez/.local/bin/

export R_LIBS_USER=/u/project/pajukant/malvarez/lib/R_%V

# For custom headers and libraries
export CPATH=$CPATH:/u/project/pajukant/malvarez/include/
export LIBRARY_PATH=$LIBRARY_PATH:/u/project/pajukant/malvarez/lib/
export PYTHONPATH=/u/project/pajukant/malvarez/.local/lib/python3.6/site-packages/:${PYTHONPATH}

Rscript diem.R

