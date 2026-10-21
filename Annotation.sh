#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --time=7-90:00:00
#SBATCH --mail-user=vvelasqz@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=annotation.%J.err
#SBATCH --output=annotation.%J.out
#SBATCH --job-name="annotation"
#SBATCH --mem=50G



source micromamba/etc/profile.d/micromamba.sh


#Assemblies
#NGS/annotation/CANO08/ragtag.scaffold.simplified.nodup2.fasta

#Masking
micromamba activate braker-env

cd NGS/annotation/ANVA32/

BuildDatabase -name persea_db ragtag.scaffold.simplified.nodup2.fasta

RepeatModeler -database persea_db -threads 12

RepeatMasker -pa 30 -lib persea_db-families.fa ragtag.scaffold.simplified.nodup2.fasta

RepeatMasker -pa 30 -lib persea_db-families.fa -xsmall -dir masked_output ragtag.scaffold.simplified.nodup2.fasta

cd NGS/annotation/ANVA32/masked_output

#ProcessRepeats -gff ragtag.scaffold.simplified.nodup2.fasta.out

awk 'NR>3 {
  seq=$5; start=$6; end=$7; strand=$9;
  rep=$10; class=$11; id=$15;
  print seq "\tRepeatMasker\trepeat\t" start "\t" end "\t.\t" strand "\t.\tID=" id ";Name=" rep ";Class=" class
}' ragtag.scaffold.simplified.nodup2.fasta.out > repeats.gff

micromamba deactivate

cd NGS/annotation/ANVA32/braker

module load braker3
module load genemark-et
module load augustus
module load samtools
module load bamtools
module load ncbi-rmblastn
module load hisat2  

export AUGUSTUS_CONFIG_PATH=augustus_config/config
export AUGUSTUS_BIN_PATH=augustus_config/bin
export AUGUSTUS_SCRIPTS_PATH=augustus_config/scripts
export PATH=$AUGUSTUS_BIN_PATH:$AUGUSTUS_SCRIPTS_PATH:$PATH
export GENEMARK_PATH=/opt/rit/el9/20230413/app/linux-rhel9-x86_64_v3/gcc-11.2.1/genemark-et-4.69-r34gcbqxoc5e4asjxjweodetom27vw57/bin/gmes_petap.pl
export BAMTOOLS_PATH=/opt/rit/el9/20230413/app/linux-rhel9-x86_64_v3/gcc-11.2.1/bamtools-2.5.2-j4iaiww5rf2aongra4ucq4apln6yl7aw/bin/
export BLAST_PATH=/opt/rit/el9/20230413/app/linux-rhel9-x86_64_v3/gcc-11.2.1/ncbi-rmblastn-2.14.0-vqnew3zgcsdjebs4cozc5uldnjq54epw/bin/

echo $AUGUSTUS_CONFIG_PATH

braker.pl --genome=NGS/annotation/ANVA32/masked_output/ragtag.scaffold.simplified.nodup2.fasta.masked --prot_seq=NGS/annotation/Viridiplantae.fa \
    --rnaseq_sets_ids=ANVA_Hoja_1,ANVA_Hoja_2,ANVA_Raiz_1,ANVA_Raiz_2,ANVA_Tallo_1,ANVA_Tallo_2 \
    --rnaseq_sets_dirs=NGS/annotation/ANVA32/RNASeq \
    --softmasking --threads=30 --species=ANVA32_inhouse \
    --AUGUSTUS_CONFIG_PATH=augustus_config/config

micromamba activate seqtools

mkdir NGS/annotation/ANVA32/liftoff_softmasking

cd NGS/annotation/ANVA32/liftoff_softmasking

micromamba activate liftoff_env

liftoff \
  -g NGS/assembly/Avocado_Hass_HiFi_Gr80_gene_annotation.gff \
  -p 30 \
  -o ANVA32_liftoff_output.gff3 \
  -copies \
  NGS/annotation/ANVA32/masked_output/ragtag.scaffold.simplified.nodup2.fasta.masked \
  NGS/assembly/Avocado_Hass_genome_HiFiAsm_Gr80.fasta

micromamba deactivate

micromamba activate seqtools

gffread ANVA32_liftoff_output.gff3 \
  -g NGS/annotation/ANVA32/masked_output/ragtag.scaffold.simplified.nodup2.fasta.masked \
  -w ANVA32_liftoff_transcripts.fa

micromamba deactivate


#busco embryophyta

cd NGS/annotation/ANVA32/braker/braker_inhouse/

mkdir busco_embryophyta_odb12_V3

micromamba activate busco_env

busco -i NGS/annotation/ANVA32/braker/braker_inhouse/braker/braker.aa -f \
 -l embryophyta_odb12 -o busco_embryophyta_odb12_V3 \
 -m proteins --download_path busco_downloads/ -c 14


#trans and rrnas
cd NGS/annotation/ANVA32/
# mkdir trna_predictions
# mkdir rrna_predictions
# Install latest tRNAscan-SE with GFF support
micromamba activate rna_annotation

# Run prediction with GFF output
cd trna_predictions
#tRNAscan-SE --gff NGS/annotation/ANVA32/ragtag.scaffold.simplified.nodup2.fasta > ANVA32_trna_predictions.gff3
# Eukaryotic genome, output results and structures
tRNAscan-SE -E -o ANVA32_trna_predictions.txt -f ANVA32_trna_structures.txt -m ANVA32_trna_summary.txt NGS/annotation/ANVA32/ragtag.scaffold.simplified.nodup2.fasta


cd ../rrna_predictions
barrnap NGS/annotation/ANVA32/ragtag.scaffold.simplified.nodup2.fasta > ANVA32_rrna_predictions.gff



