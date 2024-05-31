# Input file preparation
mkdir $actaeus/population_genetics
cp  $actaeus/GATK_ST/actaeus_outgroup.filtered.biallelic.SNPs.vcf.idx $actaeus/population_genetics/ 
cp $actaeus/GATK_ST/actaeus_outgroup.filtered.biallelic.SNPs.vcf $actaeus/population_genetics/
cd $actaeus/population_genetics/

# Remove the prefix introduced by LACE, leaving only numbers for the SuperTranscripts names
sed -i 's/Cluster\-//' actaeus_outgroup.filtered.biallelic.SNPs.vcf 
# Sort VCF by SuperTranscript name and coordinates
java -jar $Software/gatk-4.1.2.0/picard.jar SortVcf \
      I=actaeus_outgroup.filtered.biallelic.SNPs.vcf \
      O=actaeus_outgroup.filtered.biallelic.SNPs.sorted.vcf
# Redo the sorting by "chromosomes"
cat actaeus_outgroup.filtered.biallelic.SNPs.EDIT.vcf | vcf-sort --chromosomal-order > actaeus_outgroup.filtered.biallelic.SNPs.sortedOK.vcf

# Estimate LD from the correlation of genotypes and export unlinked SNPs
# The method was taken from: "the variants were linkage-disequilibrium-pruned to obtain a set of variants in approximate linkage equilibrium (unlinked sites) using the --indep-pairwise 50 5 0.2 option in PLINK v1.0.7 " *https://www.nature.com/articles/s41559-018-0717-x

# This comand remove outgroup and filter for minor allele frequencies (exclude reare alleles and invariant loci) and then searches for correlation among loci as signals of linkage dissequilibrium
# maf threshold = 0.038 calculated as 2 copies of the allele in 52 haploid genomes (as the N is 26 indivs.)
plink2 --vcf actaeus_outgroup.filtered.biallelic.SNPs.sortedOK.vcf \
--allow-extra-chr --set-all-var-ids @:#\$r,\$a --remove outgroup.txt --maf 0.04 --indep-pairwise 50 5 0.2 --bad-ld 

plink2 --vcf actaeus_outgroup.filtered.biallelic.SNPs.sortedOK.vcf \
--allow-extra-chr --set-all-var-ids @:#\$r,\$a --remove outgroup.txt --maf 0.04 --exclude plink2.prune.out --make-bed --max-alleles 2 --geno --out actaeus.filtered.biallelic.SNPs.sorted.ingroup.PLINK

#Output VCF with no missing data
plink2 --vcf actaeus_outgroup.filtered.biallelic.SNPs.sortedOK.vcf \
--allow-extra-chr --set-all-var-ids @:#\$r,\$a --remove outgroup.txt --maf 0.04 --exclude plink2.prune.out --max-alleles 2 --geno --recode vcf --out actaeus.filtered.biallelic.SNPs.sorted.ingroup
# import VCF and output VCF with 25% missing data
plink2 --vcf actaeus_outgroup.filtered.biallelic.SNPs.sortedOK.vcf \
--allow-extra-chr --set-all-var-ids @:#\$r,\$a --remove outgroup.txt --maf 0.04 --exclude plink2.prune.out --max-alleles 2 --geno 0.25 --recode vcf --out actaeus.filtered.biallelic.SNPs.sorted.ingroup.miss25


##########################################################################
#
#		Admixture calculations
#
###########################################################################

mkdir $actaeus/population_genetics/admixture
#copy all PLINK files to working dir.
cp $actaeus/population_genetics/*PLINK.* $actaeus/population_genetics/admixture
cd $actaeus/population_genetics/admixture

$Software/plink-1.07-x86_64/plink --noweb --bfile actaeus.filtered.biallelic.SNPs.sorted.ingroup.PLINK --recode12 --out actaeus.biallelic.SNPsLDprune

prefix=actaeus.biallelic.SNPsLDprune
mkdir admixture_output
awk '{print $2}' ${prefix}.nosex > individuals.txt
for r in {1..20}; do for K in {1..10}; do
$Software/dist/admixture_linux-1.3.0/admixture -s ${RANDOM} --cv ${prefix}.ped $K | tee admixture_output/log"_"K$K"_run"$r.out
mv ${prefix}.$K.Q admixture_output/${prefix}"_"K$K"_run"$r.Q
mv ${prefix}.$K.P admixture_output/${prefix}"_"K$K"_run"$r.P
done
done
grep CV admixture_output/log*.out > cross_validation_errors.txt

