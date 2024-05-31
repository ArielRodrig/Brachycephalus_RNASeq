# EvidentialGene consensus transcriptome creation
evigene=$actaeus"/ensemble_transcriptome/PerLocality/evidentialgene"
mkdir $evigene
cd $evigene
myspecies="actaeus"
assemblies=$actaeus"/ensemble_transcriptome/PerLocality/assemblies"

# Step1:	Preprocessing of assemblies: remove duplicates and rename transcripts

# SOAP-denovo multi-kmer outputs to transcripts set
$Software/evigene_new/evigene/scripts/rnaseq/trformat.pl -pre actaeus -out $evigene/SOAP_actaeus.transcripts.fa -log -in $actaeus"/ensemble_transcriptome/PerLocality/assemblies/SOAP_k"*".transcrips.fa"
# Trinity
$Software/evigene_new/evigene/scripts/rnaseq/trformat.pl -pre actaeus -out $evigene/Trinity_actaeus.transcripts.fa -log -in $actaeus"/ensemble_transcriptome/PerLocality/assemblies/trinity"*".fasta"
# SPADES
$Software/evigene_new/evigene/scripts/rnaseq/trformat.pl -pre actaeus -out $evigene/SPADES_actaeus.transcripts.fa -log -in $actaeus"/ensemble_transcriptome/PerLocality/assemblies/spades"*"transcripts.fasta"

ls $evigene/*.fa
cat $evigene/*.fa > $evigene/transcript_collection.fa
cd $evigene

# Step2:	Run evidential gene pipeline
$Software/evigene_new/evigene/scripts/prot/tr2aacds4.pl -tidy -NCPU 24 -MAXMEM 122880 -MINAA 30 -reorient -pHeterozygosity=2 -species Brachycephalus_actaeus -log -cdna transcript_collection.fa

# Completness assesment of resulting transcriptome using BUSCO
cd $evigene/okayset
busco -m prot -i transcript_collection.okay.aa.pep -o transcript_collection.okay.aa.pep.BUSCO -l vertebrata_odb10 -c 30 -f



