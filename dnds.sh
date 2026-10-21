#!/bin/sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --time=10-90:00:00
#SBATCH --mail-user=vvelasqz@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=dnds.%J.err
#SBATCH --output=dnds.%J.out
#SBATCH --job-name="dnds"
#SBATCH --mem=100G



source micromamba/etc/profile.d/micromamba.sh


cd NGS/phylogenetics


module load mafft
micromamba activate trimal_env 

selection analysis

#Load required modules (adjust to your HPC environment)
module load mafft 
module load iqtree2 
module load paml 
module load blast-plus/2.13.0-py310-irl7uxo
micromamba activate seqtools

# Directory containing orthogroup alignments
ALIGNDIR=NGS/phylogenetics/phylogenomics/alignments/
# Directory containing CDS FASTAs for all species
CDSDIR=NGS/phylogenetics/genomes/cds
# Output directory
OUTDIR=NGS/phylogenetics/dnds_results_branchmodel

mkdir -p $OUTDIR
#sed -E 's/^(>Pa[0-9]+g[0-9]+)[0-9]+/\1/' phylogenetics/genomes/cds/WI_pame_cleanV0.fa > phylogenetics/genomes/cds/WI_pame_clean.fa
#grep "^>" phylogenetics/genomes/cds/WI_pame_clean.fa | sort | uniq -d

# sed -E 's/^>.*(Super-Scaffold_[0-9]+).*/>\1/' NGS/phylogenetics/genomes/sixto_treemutation_save/01_genomes/06_genes/GCA_037013555.1_HS1_Sru_omap1_HYBRID_SCAFFOLD_genomic.fna > NGS/phylogenetics/genomes/sixto_treemutation_save/01_genomes/06_genes/GCA_037013555.1_HS1_Sru_omap1_HYBRID_SCAFFOLD_genomic_cleaned.fa

# gffread NGS/phylogenetics/genomes/sixto_treemutation_save/01_genomes/06_genes/sixto.braker.gtf \
#     -g NGS/phylogenetics/genomes/sixto_treemutation_save/01_genomes/06_genes/GCA_037013555.1_HS1_Sru_omap1_HYBRID_SCAFFOLD_genomic_cleaned.fa \
#     -x NGS/phylogenetics/genomes/cds/sixto.braker.cds.fa

# sed -E 's/\.t([0-9]+)/t\1/' NGS/phylogenetics/genomes/cds/sixto.braker.cds.fa > NGS/phylogenetics/genomes/cds/sixto.braker.cds_clean.fa


# Loop through each orthogroup alignment
for aln in $ALIGNDIR/*.aln.fa; do
    base=$(basename $aln .aln.fa)
    echo "Processing $base..."

    # 1. Extract IDs from protein alignment
    #grep ">" $aln | sed 's/>//' > $OUTDIR/${base}.ids

    # Case A: PaAn / PaCa / PaHa / PaGw → first word only
    # grep -E "^>(PaAn|PaCa|PaHa|PaGw)" $aln \
    #     | sed 's/>//' \
    #     | awk '{print $1}' > $OUTDIR/${base}.ids.default

    #if we include outgroup
    grep -E "^>(PaAn|PaCa|PaHa|PaGw|file)" $aln \
        | sed 's/>//' \
        | awk '{print $1}' > $OUTDIR/${base}.ids.default

    # Case B: Drymifolia RaGOO → drymifolia_file_1_file_1_gXXXX
    grep "^>Persea americana var drymifolia" "$aln" \
        | sed 's/>//' \
        | awk -F'\\|\\|' '{print "drymifolia_" $5}' > "$OUTDIR/${base}.ids.drymifolia"


    # Case C: Orthogroup headers → extract PaXXgYYYY
    
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
    #seqkit grep -r -f $OUTDIR/${base}.ids $CDSDIR/*_clean.fa > $OUTDIR/${base}.cds.fa

    
    # Create a cleaned alignment file with normalized IDs

    # Copy original alignment to intermediate

    #"$ALIGNDIR/${base}.trim.fa" > "$OUTDIR/${base}.trim.id.fa"
    cp "$ALIGNDIR/${base}.trim.fa" "$OUTDIR/${base}.trim.intermediate.fa"

    # Replace "Persea" header with the single Drymifolia ID
    dry_id=$(cat "$OUTDIR/${base}.ids.drymifolia")
    sed -i "s/^>Persea/>${dry_id}/" "$OUTDIR/${base}.trim.intermediate.fa"

    # Replace "XXG..." header with the single Outgroup ID
    out_id=$(cat "$OUTDIR/${base}.ids.orthogroup")
    sed -i "s/^>XXG[0-9]\+/>${out_id}/" "$OUTDIR/${base}.trim.intermediate.fa"

    # Extract only sequences with valid headers
    seqkit grep -f "$OUTDIR/${base}.ids" "$OUTDIR/${base}.trim.intermediate.fa" > "$OUTDIR/${base}.trim.id.fa"
    

    # 3. Back-translate to codon alignment
    #pal2nal.pl $OUTDIR/OG0013809.trim.id.fa $OUTDIR/OG0013809.cds.fa -output paml > OG0013809.codon.phy
    #translatorx_vLocal.pl -a $OUTDIR/OG0013809.trim.id.fa -i $OUTDIR/OG0013809.cds.fa -o OG0013809.codon
    #translatorx_vLocal.pl -i $OUTDIR/OG0013809.cds.fa -o OG0013809.codon

    #pal2nal.pl $OUTDIR/${base}.trim.id.fa $OUTDIR/${base}.cds.fa -output paml > $OUTDIR/${base}.codon.phy
    translatorx_vLocal.pl -i $OUTDIR/${base}.cds.fa -o $OUTDIR/${base}.codon

    # 4. Build gene tree
    iqtree2 -s $OUTDIR/${base}.codon.nt_ali.fasta -m MFP -nt AUTO -pre $OUTDIR/${base}

    # 5. Create codeml control file, pairwise model = 2 
    cat <<EOF > $OUTDIR/${base}.ctl
seqfile = $OUTDIR/${base}.codon.nt_ali.fasta
treefile = $OUTDIR/${base}.labeled.treefile
outfile = $OUTDIR/${base}.out

noisy = 9
verbose = 1
seqtype = 1        * 1 = codons (this is critical!)
runmode = -2

model = 0        * branch model (per species/branch ω)
NSsites = 0
icode = 0
fix_kappa = 0
fix_omega = 0
omega = 1
EOF

    # 6. Run codeml
    codeml $OUTDIR/${base}.ctl
done

echo "Pipeline finished. Results in $OUTDIR"





#summarize results

# Extract ML dN/dS vs outgroup for each OG
OUTDIR=NGS/phylogenetics/dnds_results_branchmodel
OUTGROUP_PATTERN="file"   # outgroup names start with "file"
# mkdir NGS/phylogenetics/dnds_results_branchmodel_summary/

RESULTS="NGS/phylogenetics/dnds_results_branchmodel_summary/dnds_vs_outgroup.tsv"

echo -e "OG\tAccession\tOmega\tDN\tDS\tOutgroup" > $RESULTS

for out in $OUTDIR/OG*.out; do
    OG=$(basename "$out" .out)

    # Detect the outgroup name (first sequence starting with "file")
    OUTGROUP=$(grep -o "([Ff]ile[^)]*)" "$out" \
               | head -n 1 \
               | sed -E 's/[()]//g')

    awk -v og="$OUTGROUP" -v OG="$OG" -v OUTGROUP="$OUTGROUP" '
        # Detect pairwise block header lines
        /\.\.\./ {
            # Extract left and right names explicitly
            # Example: 7 (file_1_file_1_g21007t1) ... 1 (drymifolia_file_1_file_1_g11638)
            match($0, /\(([^\)]+)\)/, L)
            match($0, /\.\.\. [0-9]+ \(([^\)]+)\)/, R)

            left = L[1]
            right = R[1]

            # Determine which one is the outgroup
            if (left == og) {
                acc = right
                inblock = 1
            } else if (right == og) {
                acc = left
                inblock = 1
            } else {
                inblock = 0
            }
        }

        # Inside block: find the line containing dN/dS
        inblock && /dN\/dS=/ {
            match($0, /dN\/dS= *([0-9.]+)/, o)
            match($0, /dN = *([0-9.]+)/, dn)
            match($0, /dS = *([0-9.]+)/, ds)

            print OG "\t" acc "\t" o[1] "\t" dn[1] "\t" ds[1] "\t" OUTGROUP
            inblock = 0
        }
    ' "$out" >> $RESULTS
done




echo -e "Orthogroup\tBranch\tOmega" > $OUTDIR/all_results.tsv

# Loop through all codeml outputs
for out in $OUTDIR/*.out; do
    base=$(basename $out .out)

    # Extract omega lines and append to master table
    grep -A1 "omega" "$out" \
        | awk -v b="$base" '{if(NR%2==1){branch=$1} else {print b "\t" branch "\t" $1}}' \
        >> $OUTDIR/all_dnds_results.tsv
done


