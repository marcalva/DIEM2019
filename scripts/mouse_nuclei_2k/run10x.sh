#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 16
#$ -l h_data=4000M,h_rt=10:00:00,highp
#$ -M malvarez@mail
#$ -v QQAPP=openmp
#$ -V
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o run10x.sh.log

. /u/local/Modules/default/init/modules.sh
module load python/3.6.1_shared
module load samtools

. /u/home/m/malvarez/.bashrc

export OMP_NUM_THREADS=16

echo $PYTHONPATH

velocyto="/u/project/pajukant/malvarez/.local/bin/velocyto"

cd ../../

out_dir="data/processed/mouse_nuclei_2k/velocyto/"
mkdir -p $out_dir
rm -rf ${out_dir}/velocyto

subject="mouse_nuclei_2k"
# bamfile="data/raw/mouse_nuclei_2k/bams/cellsorted_genome_bam.bam"
bamfile="data/raw/mouse_nuclei_2k/bams/possorted_genome_bam.bam"
bcfile="data/processed/mouse_nuclei_2k/velocyto/top10K_ids.txt"
bcfile_1="data/processed/mouse_nuclei_2k/velocyto/top10K_ids.tmp.txt"
awk '{print $1"-1"}' $bcfile > $bcfile_1

gtf="data/raw/mouse_nuclei_2k/ref/Mus_musculus.GRCm38.93.pli.gtf"

$velocyto run \
    -@ 16 \
    --sampleid $subject \
    -b $bcfile_1 \
    -o $out_dir \
    $bamfile \
    --samtools-threads 16 \
    --samtools-memory 3600 \
    $gtf

rm $bcfile_1

