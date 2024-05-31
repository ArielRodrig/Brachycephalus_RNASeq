# Script for calling variants using the RNAseq reads and taking the SuperTranscriptome as reference
# adapted from: https://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail

genomeDir=$actaeus/ensemble_transcriptome/PerLocality/Lace
workDir=$actaeus/GATK_ST
mkdir $workDir

# Align reads to SuperTranscriptome with STAR
cd $genomeDir
$Software/hisat2-2.1.0/hisat2-build $genomeDir/SuperDuper.fasta SuperTranscriptome -p 24
sort -k1,1V -k4,4n -k5,5n SuperDuper.gff > SuperDuper.sorted.gff
sort -k1,1V -k4,4n -k5,5n SuperDuperTrans.gff > SuperDuperTrans.sorted.gff

# Build SuperTranscriptome index
# The sjdbOverhang parameter is helpful in setting some internal options, and is recommended to be set as read_lenght - 1
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genomeDir \
--genomeFastaFiles SuperDuper.fasta --genomeSAindexNbases 12 \
--sjdbGTFfile $genomeDir/SuperDuperTrans.sorted.gff \
--sjdbOverhang 149 --limitGenomeGenerateRAM 100000000000


cd $workDir
for f1 in $actaeus/raw_reads/*_R1_001.P.fq; do
    f2=${f1%%_R1_001.P.fq}"_R2_001.P.fq"
    sample_name=$(basename -s _R1_001.P.fq $f1)
    echo -e "["$(date)"]\tAligning.."
    STAR --outFileNamePrefix $workDir/$sample_name --outSAMtype BAM Unsorted --outSAMstrandField intronMotif \
	--outSAMattrRGline ID:$sample_name CN:ZoolInst_TiHo LB:PairedEnd PL:Illumina PU:Unknown SM:$sample_name \
	--genomeDir $genomeDir --runThreadN 24 --readFilesIn $f1 $f2 --twopassMode Basic \
	--outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5
samtools sort -@ 12 -o $workDir/$sample_name"_sorted.bam" $workDir/$sample_name"Aligned.out.bam"
rm $workDir/$sample_name"Aligned.out.bam"
done

# Remove duplicates
cd $workDir
for srt in *_sorted.bam ; do
java -jar $Software/gatk-4.1.2.0/picard.jar MarkDuplicates I=$srt O=${srt%%.bam}".dedupped.bam" TMP_DIR=$workDir CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$srt".picard-output.metrics" 2>$workDir/$srt.MarkDuplicates.log
done
# re-run for the biggest sample which failed with defaults: add TMP_DIR=$workDir for future uses
java -jar $Software/gatk-4.1.2.0/picard.jar MarkDuplicates I=ARO01_RNA1038_sorted.bam O=ARO01_RNA1038_sorted.dedupped.bam  TMP_DIR=$workDir CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=ARO01_RNA1038_sorted.bam.picard-output.metrics 2>ARO01_RNA1038_sorted.bam.MarkDuplicates.log

# Continue with all samples now:
#SplitNCigarReads options not compatible with GATK 4 removed: # -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 
java -jar  $Software/gatk-4.1.2.0/picard.jar CreateSequenceDictionary -R $genomeDir/SuperDuper.fasta -O $genomeDir/SuperDuper.dict
cd $workDir
for ddp in *.dedupped.bam ; do
echo -e "["$(date)"]\tSpliting reads.."
$Software/gatk-4.1.2.0/gatk SplitNCigarReads -R $genomeDir/SuperDuper.fasta -I $ddp -O ${ddp%%.bam}".split.bam" 2>$ddp.SplitNCigarReads.log
done

# 6. Variant calling

cd $workDir
samtools faidx $genomeDir/SuperDuper.fasta 
java -jar  $Software/gatk-4.1.2.0/picard.jar CreateSequenceDictionary R=$genomeDir/SuperDuper.fasta O=$genomeDir/SuperDuper.dict

cd $workDir
$Software/gatk-4.1.2.0/gatk HaplotypeCaller --java-options '-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=28' --native-pair-hmm-threads 24 -R $genomeDir/SuperDuper.fasta \
-I $actaeus/GATK_ST/ARO01_RNA1021_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1022_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1023_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1024_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1030_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1032_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1033_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1034_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1035_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1037_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1038_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1039_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1040_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1041_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1043_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1045_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1046_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1047_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1048_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1050_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1051_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1052_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1055_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1059_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1060_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1061_sorted.dedupped.split.bam \
-I $actaeus/GATK_ST/ARO01_RNA1062_sorted.dedupped.split.bam \
--dont-use-soft-clipped-bases false --standard-min-confidence-threshold-for-calling 20.0 -O actaeus_outgroup.raw.vcf

# The avove comand completed.
# ProgressMeter - Traversal complete. Processed 297363 total regions in 3919.0 minutes.

# Variant Filtration

cd $workDir
java -jar $Software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration -R $genomeDir/SuperDuper.fasta -V actaeus_outgroup.raw.vcf -window 35 -cluster 3 -filter "FS > 30.0 || QD < 2.0" --filter-name FSQD -O actaeus_outgroup.filtered.vcf 2>actaeus_outgroup.filtered.VariantFilter.log
java -jar $Software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants -R $genomeDir/SuperDuper.fasta -V actaeus_outgroup.filtered.vcf -O actaeus_outgroup.filtered.biallelic.vcf --restrict-alleles-to BIALLELIC --set-filtered-gt-to-nocall true
java -Xmx60g -Xms60g -XX:+UseParallelGC -jar $Software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
-R $genomeDir/SuperDuper.fasta -V actaeus_outgroup.filtered.biallelic.vcf --select-type-to-include SNP -O actaeus_outgroup.filtered.biallelic.SNPs.vcf


# Divide into SNPs and Indels
java -Xmx60g -Xms60g -XX:+UseParallelGC -jar $Software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
-R $genomeDir/SuperDuper.fasta -V actaeus_outgroup.filtered.biallelic.vcf --select-type-to-include SNP -O actaeus_outgroup.filtered.biallelic.SNPs.vcf 

java -Xmx60g -Xms60g -XX:+UseParallelGC -jar $Software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
-R $genomeDir/SuperDuper.fasta -V actaeus_outgroup.filtered.biallelic.vcf --select-type-to-include INDEL -O actaeus_outgroup.filtered.biallelic.INDELS.vcf 

java -Xmx60g -Xms60g -XX:+UseParallelGC -jar $Software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
-R $genomeDir/SuperDuper.fasta -V actaeus_outgroup.raw.vcf --select-type-to-include SNP -O actaeus_outgroup.raw.SNPs.vcf 

java -Xmx60g -Xms60g -XX:+UseParallelGC -jar $Software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
-R $genomeDir/SuperDuper.fasta -V actaeus_outgroup.raw.vcf --select-type-to-include INDEL -O actaeus_outgroup.raw.INDELS.vcf 

$Software/gatk-4.1.2.0/gatk VariantsToTable -V actaeus_outgroup.raw.SNPs.vcf \
             -F CHROM -F POS -F TYPE -GF AD \
             -O actaeus_outgroup.raw.SNPs.tsv 2> actaeus_outgroup.raw.SNPs.tsv.log
$Software/gatk-4.1.2.0/gatk VariantsToTable -V actaeus_outgroup.raw.INDELS.vcf \
             -F CHROM -F POS -F TYPE -GF AD \
             -O actaeus_outgroup.raw.INDELS.tsv 2> actaeus_outgroup.raw.INDELS.tsv.log

# These outputs contain the outgroup samples. These were removed later for population genetic analyses