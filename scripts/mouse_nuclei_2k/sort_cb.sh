#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4000M,h_rt=6:00:00,highp
#$ -M malvarez@mail
#$ -v QQAPP=openmp
#$ -V
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o sort_cb.sh.log

. /u/local/Modules/default/init/modules.sh
module load samtools
. /u/home/m/malvarez/.bashrc

export OMP_NUM_THREADS=8

cd ../../

bamdir="data/raw/mouse_nuclei_2k/bams/"

bamfile="${bamdir}/possorted_genome_bam.bam"
bamout="${bamdir}/cellsorted_genome_bam.bam"

# remove tmp files if any
rm -rf ${bamout}.tmp*

samtools sort \
    -l 7 \
    -m 3G \
    -t CB \
    -O BAM \
    -@ 8 \
    -o $bamout \
    $bamfile

