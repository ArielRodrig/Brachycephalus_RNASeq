# Estimate gene expression with callisto

# prepare reference transcriptome (all tissues)
WORKDIR=$actaeus"/kallisto_quantification"
mkdir $WORKDIR

cp $actaeus/ensemble_transcriptome/PerLocality/entap_outfiles/frame_selection/TransDecoder/processed/ORF-containing_transcripts_no_contam.fasta $WORKDIR
cd $WORKDIR

# Pseudo-quantification with kallisto, build the index first
$Software/kallisto_linux-v0.44.0/kallisto index -i evigene_corset.idx  ORF-containing_transcripts_no_contam.fasta

# This will list all the samples

cd $actaeus/clean_reads
samples=$(ls *_R1_001.P.fq | sed "s/_R1_001.P.fq//g" | uniq )

# Make directories for each sample:
cd $WORKDIR
mkdir kallisto_quants
cd kallisto_quants
for i in $samples; do mkdir $i; done
cd ..

# Perform the actual pseudo-quantification for all of the samples 
readsdir=$actaeus/clean_reads
parallel -j 16 /data/biolinux/Software/kallisto_linux-v0.44.0/kallisto quant -i evigene_corset.idx -o kallisto_quants/{} -b 100 \
$readsdir{}_R1_001.P.fq $readsdir{}_R2_001.P.fq ::: $samples
done
#######


