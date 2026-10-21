#!/bin/sh

#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --time=90:00:00
#SBATCH --mail-user=vvelasqz@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=QC.%J.err
#SBATCH --output=QC.%J.out
#SBATCH --job-name="QC"
#SBATCH --mem=50G


source micromamba/etc/profile.d/micromamba.sh

#start with fastqc and multiqc

RUN_PATH=NGS/annotation/ANVA32/RNASeq/
OUT_FLDR=NGS/annotation/ANVA32/RNASeq_QC/

mkdir ${OUT_FLDR}
mkdir ${OUT_FLDR}fastqc
mkdir ${OUT_FLDR}multiqc

module load fastqc
for file in $(ls ${RUN_PATH}*fastq)
do
    SAMPLE=`basename $file`
    fastqc -t 30 $RUN_PATH/${SAMPLE} -o ${OUT_FLDR}fastqc
done

micromamba activate seqQC
cd ${OUT_FLDR}multiqc
multiqc ${OUT_FLDR}fastqc -o ${OUT_FLDR}multiqc


#CANO

RUN_PATH=NGS/annotation/CANO08/RNASeq/
OUT_FLDR=NGS/annotation/CANO08/RNASeq_QC/

mkdir ${OUT_FLDR}
mkdir ${OUT_FLDR}fastqc
mkdir ${OUT_FLDR}multiqc

module load fastqc
for file in $(ls ${RUN_PATH}*fastq)
do
    SAMPLE=`basename $file`
    fastqc -t 30 $RUN_PATH/${SAMPLE} -o ${OUT_FLDR}fastqc
done


#module load multiqc
cd ${OUT_FLDR}multiqc
multiqc ${OUT_FLDR}fastqc -o ${OUT_FLDR}multiqc

