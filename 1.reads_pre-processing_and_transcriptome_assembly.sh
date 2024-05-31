# --------------- Preprocessing --------------#
# concatenating reads per locality

mkdir $actaeus/ensemble_transcriptome/PerLocality/reads
cat $actaeus/raw_reads/ARO01_RNA1021_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1022_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1023_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1024_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1030_R1_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/Morro_R1_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1021_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1022_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1023_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1024_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1030_R2_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/Morro_R2_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1032_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1033_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1034_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1035_R1_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/Cacho_R1_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1032_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1033_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1034_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1035_R2_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/Cacho_R2_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1037_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1038_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1039_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1040_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1041_R1_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/Palha_R1_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1037_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1038_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1039_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1040_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1041_R2_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/Palha_R2_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1043_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1045_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1046_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1047_R1_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/BracoGreen_R1_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1043_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1045_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1046_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1047_R2_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/BracoGreen_R2_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1048_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1050_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1051_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1052_R1_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/BracoBrown_R1_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1048_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1050_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1051_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1052_R2_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/BracoBrown_R2_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1059_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1060_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1061_R1_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1062_R1_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/Casarao_R1_001.P.fq

cat $actaeus/raw_reads/ARO01_RNA1059_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1060_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1061_R2_001.P.fq \
$actaeus/raw_reads/ARO01_RNA1062_R2_001.P.fq \
> $actaeus/ensemble_transcriptome/PerLocality/reads/Casarao_R2_001.P.fq

cd $actaeus/ensemble_transcriptome/PerLocality/reads

# Normalizing read pools to 50X
for f1 in *_R1_001.P.fq; do
    f2=${f1%%_R1_001.P.fq}"_R2_001.P.fq"
$Software/bbmap/bbnorm.sh -Xmx100g in1=$f1 in2=$f2 out1=${f1%%.fq}".norm.fq" out2=${f2%%.fq}".norm.fq" target=50 mindepth=2;
for file in *.norm.fq; do
gzip $file; done
rm $f1
rm $f2
done

# rRNA de-contamination with SortMeRNA
WORKDIR="$actaeus/ensemble_transcriptome/PerLocality/sortmerna"
mkdir $WORKDIR
for f1 in *_R1_001.P.norm.fq.gz; do
    f2=${f1%%_R1_001.P.norm.fq.gz}"_R2_001.P.norm.fq.gz"
$Software/sortmerna/bin/sortmerna --ref $Software/sortmerna/database/smr_v4.3_fast_db.fasta \
--reads $f1 \
--reads $f2 \
--threads 32 -workdir $WORKDIR -blast 1 \
-num_alignments 1 -v 1 -fastx 1 -paired_in 1 --zip-out --out2 1 \
--aligned ${f1%%_R1_001.P.norm.fq.gz}".rRNAcontaminated" \
--other ${f1%%_R1_001.P.norm.fq.gz}".rRNAfiltered"
rm -r "$WORKDIR/"readb
rm -r "$WORKDIR/"kvdb
rm -r "$WORKDIR/"idx
done

########################
#       ASSEMBLIES     #
########################


localities="Cacho Morro Palha BracoBrown BracoGreen Casarao"
for loc in $localities; do
READSDIR=$actaeus"/ensemble_transcriptome/PerLocality/reads"
WORKDIR=$actaeus"/ensemble_transcriptome/PerLocality/assemblies"
READ1=$READSDIR/$loc".rRNAfiltered_fwd.fq.gz"
READ2=$READSDIR/$loc".rRNAfiltered_rev.fq.gz"
TRINITYDIR=$actaeus"/ensemble_transcriptome/PerLocality/assemblies/trinity_temp"
mkdir $TRINITYDIR
cd $TRINITYDIR
$Software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left $READ1  --right $READ2 --min_kmer_cov 2 --SS_lib_type RF --CPU 30 --max_memory 100G --bflyCalculateCPU --full_cleanup --no_normalize_reads --output "trinity.$loc"
cp "trinity.$loc.Trinity.fasta" $WORKDIR
rm -r $TRINITYDIR
done

# rnaSPADES 
localities="Cacho Morro Palha BracoBrown BracoGreen Casarao"
for loc in $localities; do
READSDIR=$actaeus"/ensemble_transcriptome/PerLocality/reads"
WORKDIR=$actaeus"/ensemble_transcriptome/PerLocality/assemblies"
GZIP_READ1=$READSDIR/$loc".rRNAfiltered_fwd.fq.gz"
GZIP_READ2=$READSDIR/$loc".rRNAfiltered_rev.fq.gz"
gunzip -k $GZIP_READ1
gunzip -k $GZIP_READ2
READ1=${GZIP_READ1%%.fq.gz}".fq"
READ2=${GZIP_READ2%%.fq.gz}".fq"
SPADESDIR=$actaeus"/ensemble_transcriptome/PerLocality/assemblies/spades_temp"
mkdir $SPADESDIR
cd $WORKDIR
rnaspades.py -k auto -o $SPADESDIR --ss rf -1 $READ1 -2 $READ2 -m 124 -t 24 --disable-gzip-output  --checkpoints "last"
cp $SPADESDIR/"transcripts.fasta" $WORKDIR"/spades_"$loc".transcripts.fasta"
rm -r $SPADESDIR
rm $READ1
rm $READ2
done

# OASES

localities="Cacho Morro Palha BracoBrown BracoGreen Casarao"
for loc in $localities; do
READSDIR=$actaeus"/ensemble_transcriptome/PerLocality/reads"
WORKDIR=$actaeus"/ensemble_transcriptome/PerLocality/assemblies"
GZIP_READ1=$READSDIR/$loc".rRNAfiltered_fwd.fq.gz"
GZIP_READ2=$READSDIR/$loc".rRNAfiltered_rev.fq.gz"
gunzip -k $GZIP_READ1
gunzip -k $GZIP_READ2
READ1=${GZIP_READ1%%.fq.gz}".fq"
READ2=${GZIP_READ2%%.fq.gz}".fq"
TEMPDIR=$actaeus"/ensemble_transcriptome/PerLocality/assemblies/oases_temp"
mkdir $TEMPDIR
kmerValues="19 37 55 73"
for kmer in $kmerValues; do
OASESDIR=$TEMPDIR"/oases_k"$kmer
mkdir $OASESDIR
cd $OASESDIR
echo VELVETH_START; velveth $OASESDIR $kmer -shortPaired -fastq -separate -strand_specific $READ1 $READ2; echo VELVETG_START; velvetg $OASESDIR -read_trkg yes -min_contig_lgth 100 -cov_cutoff 4 -ins_length 400 -clean yes; echo OASES_START; oases $OASESDIR -cov_cutoff 4 -min_trans_lgth 200
cp $OASESDIR/"transcripts.fa" $WORKDIR"/oases_k"$k"_"$loc".transcrips.fa"
rm -r $OASESDIR
done
rm $READ1
rm $READ2
done


# SOAPdenovo-trans

localities="Cacho Morro Palha BracoBrown BracoGreen Casarao"
for loc in $localities; do
READSDIR=$actaeus"/ensemble_transcriptome/PerLocality/reads"
WORKDIR=$actaeus"/ensemble_transcriptome/PerLocality/assemblies"
GZIP_READ1=$READSDIR/$loc".rRNAfiltered_fwd.fq.gz"
GZIP_READ2=$READSDIR/$loc".rRNAfiltered_rev.fq.gz"
gunzip -k $GZIP_READ1
gunzip -k $GZIP_READ2
READ1=${GZIP_READ1%%.fq.gz}".fq"
READ2=${GZIP_READ2%%.fq.gz}".fq"
TEMPDIR=$actaeus"/ensemble_transcriptome/PerLocality/assemblies/soap_temp"
mkdir $TEMPDIR
kmerValues="19 37 55 73"
for kmer in $kmerValues; do
SOAPDIR=$TEMPDIR"/soap_k"$kmer
mkdir $SOAPDIR
cd $SOAPDIR
echo -e "#maximal read length\nmax_rd_len=150\n[LIB]\n#maximal read length in this lib\nrd_len_cutof=150\n#average insert size\navg_ins=400\n#if sequence needs to be reversed \nreverse_seq=0\n#in which part(s) the reads are used\nasm_flags=3\n#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\nmap_len=32\n#fastq file for read 1\nq1=$READ1\n#fastq file for read 2\nq2=$READ2\n" > SOAP_Conf.txt

SOAPdenovo-Trans-127mer all -s SOAP_Conf.txt -K $kmer -o SOAP_k$kmer.$loc -F -p 30
cp *.scafSeq $WORKDIR/SOAP_k$kmer.$loc.transcrips.fa
rm -r $SOAPDIR
done
rm $READ1
rm $READ2
done


