#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --time=7-90:00:00
#SBATCH --mail-user=vvelasqz@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=phylo.%J.err
#SBATCH --output=phylo.%J.out
#SBATCH --job-name="phylo"
#SBATCH --mem=100G


source micromamba/etc/profile.d/micromamba.sh

# Load modules (adapt to your HPC environment)
module load orthofinder/2.5.4-py310-4cstbgx
module load iqtree2/2.2.2-py310-crc5oqm
module load diamond/2.0.15-agvctiy
module load mafft
module load bcftools
module load samtools
module load psmc
module load paml
module load beast


cd NGS/phylogenetics

micromamba activate seqtools
cd genomes

sed 's/^>lcl|/>/' Persea_americana_var._drymifolia.faa > Persea_americana_var._drymifolia_headers.faa

gffread Persea_americana_var._drymifolia_annos1-cds1-id_typename-nu1-upa1-add_chr0.gid68395.gff \
  -g Persea_americana_var._drymifolia_headers.faa \
  -y Persea_americana_var._drymifolia_proteins.fa

  gffread Persea_americana_var._drymifolia_annos1-cds0-id_typenum-nu1-upa1-add_chr0.gid68395.gff \
  -g Persea_americana_var._drymifolia_headers.faa \
  -y Persea_americana_var._drymifolia_proteins.fa

wget https://zenodo.org/record/10089692/files/treemutation_save.zip?download=1 -O treemutation_save.zip
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/037/013/555/GCA_037013555.1_HS1_Sru_omap1_HYBRID_SCAFFOLD/GCA_037013555.1_HS1_Sru_omap1_HYBRID_SCAFFOLD_genomic.fna.gz  

gffread sixto.braker.gtf \
  -g NGS/phylogenetics/genomes/treemutation_save/01_genomes/05_te/HS1_Sru_omap1_hap1_HYBRID_SCAFFOLD.fa.masked \
  -y sixto_braker_proteins.fa


########################################
# 1. Phylogenomics
########################################

orthofinder -f NGS/phylogenetics/genomes/ \
            -t 16 -a 16 \
            -o NGS/phylogenetics/phylogenomics/orthofinder_results/



module load mafft
micromamba activate trimal_env 
module load python  # for AMAS if needed

Input folder from OrthoFinder
OG_DIR="NGS/phylogenetics/phylogenomics/orthofinder_results/Results_Apr07/Single_Copy_Orthologue_Sequences/"
OUTDIR="NGS/phylogenetics/phylogenomics/alignments/"
mkdir -p $OUTDIR
mkdir -p $OUTDIR/filtered_alignments/

for og in $OG_DIR/*.fa; do
    base=$(basename $og .fa)

    echo "Processing $base"

    # 1. Align with MAFFT
    mafft --maxiterate 1000 --localpair $og > $OUTDIR/${base}.aln.fa

    # 2. Trim with trimAl
    trimal -automated1 \
        -in $OUTDIR/${base}.aln.fa \
        -out $OUTDIR/${base}.trim.fa
done

echo "All alignments and trimming complete."

3. Concatenate using AMAS (recommended)
echo "Concatenating alignments..."

AMAS.py concat \
    -i $OUTDIR/*.trim.fa \
    -f fasta \
    -d aa \
    -p partitions.txt \
    -t supermatrix.fa

fix concat with renamed alin files 
ALIGNDIR=NGS/phylogenetics/phylogenomics/alignments/

for f in "$ALIGNDIR"/*.trim.fa; do
    base=$(basename "$f" .trim.fa)

    sed \
        -e 's/^>file.*/>Sixto/' \
        -e 's/^>.*Persea.*/>PaDr/' \
        -e 's/^>XXG.*/>PaWI/' \
        -e 's/^>.*PaAn.*/>PaAn/' \
        -e 's/^>.*PaCa.*/>PaCa/' \
        -e 's/^>.*PaHa.*/>PaHa/' \
        -e 's/^>.*PaGw.*/>PaGw/' \
        "$f" > "$OUTDIR/filtered_alignments/${base}.trim.renamed.fa"
done


AMAS.py concat \
    -i $OUTDIR/filtered_alignments/*trim.renamed.fa \
    -f fasta \
    -d aa \
    -p partitions_filtered_renamed.txt \
    -t supermatrix_filtered_renamed.fa

echo "Concatenation complete."
echo "Output:"
echo "  supermatrix.fa"
echo "  partitions.txt"

sed -E 's/^(p[0-9]+_OG[0-9]+ =)/AA, \1/' partitions.txt > partitions_fixed.txt

sed -E 's/^(p[0-9]+_OG[0-9]+ =)/AA, \1/' partitions_filtered.txt > partitions_fixed_filtered.txt

sed -E 's/^(p[0-9]+_OG[0-9]+ =)/AA, \1/' partitions_filtered_renamed.txt > partitions_filtered_renamed_fixed.txt

module load iqtree2/2.2.2-py310-crc5oqm
1. Create an output directory
mkdir -p iqtree_out

# 2. Run IQ-TREE with your supermatrix and partitions
iqtree2 -s supermatrix_filtered.fa \
        -p partitions_fixed_filtered.txt \
        -m MFP \
        -bb 1000 \
        -nt AUTO \
        -pre iqtree_out/supermatrix \
        #-alrt 100 \


mkdir -p iqtree_out_renamed

# 2. Run IQ-TREE with your supermatrix and partitions
iqtree2 \
  -s supermatrix_filtered_renamed.fa \
  -p partitions_filtered_renamed_fixed.txt \
  -m MFP \
  -bb 1000 \
  -nt 28 \
  -pre iqtree_out_renamed/supermatrix

iqtree2 -s supermatrix_filtered.fa -p partitions_fixed_filtered.txt -m TESTONLY

iqtree2 -s supermatrix.fa \
        -p partitions_fixed.txt \
        -m MFP \
        -bb 1000 \
        -nt AUTO \
        -pre iqtree_out/supermatrix


#make cds alignment and tree 
micromamba activate seqtools

ALIGNDIR=NGS/phylogenetics/phylogenomics/alignments/
# Directory containing CDS FASTAs for all species
CDSDIR=NGS/phylogenetics/genomes/cds
OUTDIR="NGS/phylogenetics/phylogenomics/cds_alignments"
mkdir -p $OUTDIR
# mkdir NGS/phylogenetics/phylogenomics/alignments/cds_tree

# Loop through each orthogroup alignment
for aln in $ALIGNDIR/*.aln.fa; do
    base=$(basename $aln .aln.fa)
    echo "Processing $base..."

    #if we include outgroup
    grep -E "^>(PaAn|PaCa|PaHa|PaGw|file)" $aln \
        | sed 's/>//' \
        | awk '{print $1}' > $OUTDIR/${base}.ids.default

    grep "^>Persea americana var drymifolia" "$aln" \
        | sed 's/>//' \
        | awk -F'\\|\\|' '{print "drymifolia_" $5}' > "$OUTDIR/${base}.ids.drymifolia"
    
    grep "^>XXG" "$aln" \
        | sed 's/>//' \
        | awk '{
           match($0,/Pa[0-9]+g[0-9]+(t[0-9]+)?/)
           if(RSTART>0){print substr($0,RSTART,RLENGTH)}
        }' > "$OUTDIR/${base}.ids.orthogroup"


    # Combine all cases into one .ids file
    cat $OUTDIR/${base}.ids.default \
        $OUTDIR/${base}.ids.drymifolia \
        $OUTDIR/${base}.ids.orthogroup \
        > $OUTDIR/${base}.ids
    
    # 2. Collect CDS sequences for those IDs across species
    seqkit grep -f $OUTDIR/${base}.ids $CDSDIR/*_clean.fa > $OUTDIR/${base}.cds.fa
done

module load mafft
micromamba activate trimal_env 

cp NGS/phylogenetics/dnds_results_branchmodel/*.cds.fa $OUTDIR/.

for og in $OUTDIR/*.cds.fa; do
    base=$(basename $og .cds.fa)

    echo "Processing $base and $og"

    # 1. Align with MAFFT
    mafft --maxiterate 1000 --localpair $og > $OUTDIR/${base}.cds.aln.fa

    # 2. Trim with trimAl
    trimal -automated1 \
        -in $OUTDIR/${base}.cds.aln.fa \
        -out $OUTDIR/${base}.cds.trim.fa
    
    sed \
        -e 's/^>file.*/>Sixto/' \
        -e 's/^>.*drymifolia.*/>PaDr/' \
        -e 's/^>Pa[0-9][0-9]g[0-9][0-9]*/>PaWI/' \
        -e 's/^>.*PaAn.*/>PaAn/' \
        -e 's/^>.*PaCa.*/>PaCa/' \
        -e 's/^>.*PaHa.*/>PaHa/' \
        -e 's/^>.*PaGw.*/>PaGw/' \
        "$OUTDIR/${base}.cds.trim.fa" > "$OUTDIR/${base}.cds.trim.renamed.fa"

done

module load iqtree2/2.2.2-py310-crc5oqm
# mkdir -p iqtree_out_cds_renamed
cd iqtree_out_cds_renamed

AMAS.py concat \
    -i $OUTDIR/*trim.renamed.fa \
    -f fasta \
    -d dna \
    -p partitions_cds_filtered_renamed.txt \
    -t supermatrix_cds_filtered_renamed.fa


sed -E 's/^(p[0-9]+_OG[0-9]+ =)/DNA, \1/' partitions_cds_filtered_renamed.txt > partitions_cds_filtered_renamed_fixed.txt


iqtree2 \
  -s supermatrix_cds_filtered_renamed.fa \
  -p partitions_cds_filtered_renamed_fixed.txt \
  -m MFP \
  -bb 1000 \
  -nt 28 \
  -pre supermatrix



# ============================
# 2. TREETIME DATING
# ============================


treetime --tree "$OUTDIR/species_tree_rooted.nwk" \
         --aln "$ALIGN" \
         --clock-rate 0 \
         --dates NGS/phylogenetics/iqtree_out_cds_renamed_dating_results/dates.csv \
         --name-column name \
         --date-column date \
         --outdir "$OUTDIR/treetime_relative"


if [ -n "$RATE" ]; then

    treetime --tree "$OUTDIR/species_tree_rooted.nwk" \
             --aln "$ALIGN" \
             --clock-rate "$RATE" \
             --dates NGS/phylogenetics/iqtree_out_cds_renamed_dating_results/dates.csv \
             --name-column name \
             --date-column date \
             --outdir "$OUTDIR/treetime_rate_${RATE}"

fi


