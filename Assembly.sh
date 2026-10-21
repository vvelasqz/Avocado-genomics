#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --time=90:00:00
#SBATCH --mail-user=vvelasqz@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=assembly.%J.err
#SBATCH --output=assembly.%J.out
#SBATCH --job-name="Assembly"
#SBATCH --mem=350G


source micromamba/etc/profile.d/micromamba.sh

#cp NGS/reads/ANVA32-2_CAN08-2_Pooled/HiFi/PR-28_Avocado_ANVA32-2_CAN08-2_Pooled_HiFi_300pM-Cell4_\(all_samples\)/6950_import-collection/hifi_reads/ANVA32.fastq.gz NGS/assembly/.
#cp NGS/reads/ANVA32-2_CAN08-2_Pooled/HiFi/PR-28_Avocado_ANVA32-2_CAN08-2_Pooled_HiFi_300pM-Cell4_\(all_samples\)/6950_import-collection/hifi_reads/CANO08.fastq.gz NGS/assembly/.


micromamba activate nextflow
module load singularity

cd NGS/assembly/

nextflow run nf-core/genomeassembler \
   -profile hifi_hifiasm,singularity \
   --input NGS/assembly/samplesheet.csv \
   --outdir NGS/assembly/ref_asembly \
   -params-file NGS/assembly/params.yaml


nextflow run nf-core/genomeassembler \
   -profile hifi_hifiasm,singularity -resume\
   --input NGS/assembly/samplesheet2.csv \
   --outdir NGS/assembly/ref_asembly2 \
   --skip_assembly true --skip_alignments true
   -params-file NGS/assembly/params.yaml


run busco on the side with the fist try assemblies

module load busco

micromamba activate busco_env

gunzip NGS/assembly/ref_asembly/ANVA32/assembly/hifiasm/fasta/ANVA32.asm.bp.p_ctg.fa.gz

busco -i NGS/assembly/ref_asembly/ANVA32/assembly/hifiasm/fasta/ANVA32.asm.bp.p_ctg.fa \
 -f -l eudicotyledons_odb12 -o NGS/assembly/ref_asembly/ANVA32/QC/busco_eudicotyledons -m genome -c 14

gunzip NGS/assembly/ref_asembly/CANO08/assembly/hifiasm/fasta/CANO08.asm.bp.p_ctg.fa.gz
busco -i NGS/assembly/ref_asembly/CANO08/assembly/hifiasm/fasta/CANO08.asm.bp.p_ctg.fa \
 -f -l eudicotyledons_odb12 -o NGS/assembly/ref_asembly/CANO08/QC/busco_eudicotyledons -m genome -c 14

quast
micromamba activate quast-env

quast NGS/assembly/ref_asembly/ANVA32/assembly/hifiasm/fasta/ANVA32.asm.bp.p_ctg.fa.gz \
-o NGS/assembly/ref_asembly/ANVA32/QC/quast_output

quast NGS/assembly/ref_asembly/CANO08/assembly/hifiasm/fasta/CANO08.asm.bp.p_ctg.fa.gz \
-o NGS/assembly/ref_asembly/CANO08/QC/quast_output

#move busco downloads 

#mv ~/busco_downloads/ .

#run mabs to optimize assemblies
micromamba activate mabs-env

python Mabs/mabs-hifiasm.py --long_reads NGS/assembly/ANVA32.fastq.gz \
   --local_busco_dataset busco_downloads/lineages/eudicotyledons_odb12.tar.gz --threads 30 \
   --output_folder NGS/assembly/ref_asembly/ANVA32/mabs_output_ANVA32

python Mabs/mabs-hifiasm.py --long_reads NGS/assembly/CANO08.fastq.gz \
   --local_busco_dataset busco_downloads/lineages/eudicotyledons_odb12.tar.gz --threads 30 \
   --output_folder NGS/assembly/ref_asembly/CANO08/mabs_output_CANO08

download busco lineage files
python Mabs/mabs-hifiasm.py --long_reads NGS/assembly/ANVA32.fastq.gz \
   --download_busco_dataset eudicots_odb10.2024-01-08.tar.gz --threads 30 \
   --output_folder NGS/assembly/ref_asembly/ANVA32/mabs_output_ANVA32_2

python Mabs/mabs-hifiasm.py --long_reads NGS/assembly/CANO08.fastq.gz \
   --download_busco_dataset eudicots_odb10.2024-01-08.tar.gz --threads 30 \
   --output_folder NGS/assembly/ref_asembly/CANO08/mabs_output_CANO08_2



micromamba activate scaffolding

cd NGS/assembly/ref_asembly/ANVA32/longstitch_output

#longstich soft
micromamba activate scaffolding

cd NGS/assembly/ref_asembly/ANVA32/longstitch_output

longstitch ntLink \
  draft=ANVA32_assembly \
  reads=ANVA32 \
  k_ntLink=28 \
  w=80 \
  z=300 \
  max_gap=100000 \
  span=20 \
  threads=30


cd NGS/assembly/ref_asembly/CANO08/longstitch_output

longstitch ntLink \
  draft=CANO08_assembly \
  reads=CANO08 \
  k_ntLink=28 \
  w=80 \
  z=300 \
  max_gap=100000 \
  span=20 \
  threads=30

#ragtag
micromamba activate scaffolding
ragtag.py scaffold \
  -o NGS/assembly/ref_asembly/ANVA32/ragtag_output \
  -t 30 -u \
  --aligner nucmer \
  --nucmer-params="--maxmatch -l 100 -c 500" \
  --remove-small \
  NGS/assembly/Avocado_Hass_genome_HiFiAsm_Gr80.fasta \
  NGS/assembly/ref_asembly/ANVA32/longstitch_output/ANVA32_assembly.fa.k28.w80.z300.ntLink.scaffolds.fa


ragtag.py scaffold \
  -o NGS/assembly/ref_asembly/CANO08/ragtag_output \
  -t 30 -u \
  --aligner nucmer \
  --nucmer-params="--maxmatch -l 100 -c 500" \
  --remove-small \
  NGS/assembly/Avocado_Hass_genome_HiFiAsm_Gr80.fasta \
  NGS/assembly/ref_asembly/CANO08/longstitch_output/CANO08_assembly.fa.k28.w80.z300.ntLink.scaffolds.fa


micromamba activate LTR_retriever
cd NGS/assembly/ref_asembly/CANO08/QC/ragtag_ltr
ltr_finder NGS/assembly/ref_asembly/CANO08/ragtag_output/ragtag.scaffold.simplified.nodup2.fasta > CANO08_ragtag_finder.scn
LTR_retriever -genome NGS/assembly/ref_asembly/CANO08/ragtag_output/ragtag.scaffold.simplified.nodup2.fasta -infinder CANO08_ragtag_finder.scn

#sed '/^>/ s/_.*//' NGS/assembly/ref_asembly/ANVA32/ragtag_output/ragtag.scaffold.fasta > NGS/assembly/ref_asembly/ANVA32/ragtag_output/ragtag.scaffold.simplified.fasta


#merqury
micromamba activate LTR_retriever
#mkdir NGS/assembly/ref_asembly/ANVA32/QC/merqury
cd NGS/assembly/ref_asembly/ANVA32/QC/merqury

#$MERQURY/merqury.sh F1.k18.meryl NGS/assembly/ref_asembly/ANVA32/ragtag_output/ragtag.scaffold.simplified.fasta ANVA32

meryl count k=21 output read_k21.meryl NGS/assembly/ANVA32.fastq.gz

$MERQURY/merqury.sh read_k21.meryl NGS/assembly/ref_asembly/ANVA32/ragtag_output/ragtag.scaffold.simplified.fasta ANVA32


#mkdir NGS/assembly/ref_asembly/CANO08/QC/merqury
cd NGS/assembly/ref_asembly/CANO08/QC/merqury

meryl count k=21 output read_k21.meryl NGS/assembly/CANO08.fastq.gz

$MERQURY/merqury.sh read_k21.meryl NGS/assembly/ref_asembly/CANO08/ragtag_output/ragtag.scaffold.simplified.fasta CANO08

$MERQURY/merqury.sh F1.k18.meryl NGS/assembly/ref_asembly/CANO08/ragtag_output/ragtag.scaffold.simplified.fasta CANO08

