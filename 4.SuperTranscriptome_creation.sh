#corset
cd $actaeus/ensemble_transcriptome/PerLocality/entap_outfiles/frame_selection/TransDecoder/processed
cat *_genes.fnn > ORF-containing_transcripts.fa

awk 'BEGIN{while((getline<$actaeus"/ensemble_transcriptome/PerLocality/entap_outfiles/final_results/contaminants_ids.list")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' ORF-containing_transcripts.fa > ORF-containing_transcripts_no_contam.fasta

mkdir $actaeus/ensemble_transcriptome/PerLocality/corset
cd $actaeus/ensemble_transcriptome/PerLocality/corset
salmon index --index entap_coding_nocontam --threads 24 --transcripts $actaeus/ensemble_transcriptome/PerLocality/entap_outfiles/frame_selection/TransDecoder/processed/ORF-containing_transcripts_no_contam.fasta

#sleep 4h
FILES=`ls $actaeus/raw_reads/*_R1_001.P.fq | sed 's/_R1_001.P.fq//g'`
for F in $FILES ; do
        R1=${F}_R1_001.P.fq
        R2=${F}_R2_001.P.fq
        salmon quant --index entap_coding_nocontam --libType ISR --dumpEq --hardFilter --skipQuant -1 $R1 -2 $R2 --output ${F}.out
done

gunzip -k -f $actaeus/raw_reads/ARO01_RNA*.out/aux_info/eq_classes.txt.gz

cd $actaeus/ensemble_transcriptome/PerLocality/corset
# Run corset without the log-likelihood ratio test to allow for subsequent differential transcript expression analysis (setting -D very high)
$Software/corset-1.09/corset -D 99999999999 -g 1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,2,2,2,2,4,4,4,4 -n JPCM1021,JPCM1022,JPCM1023,JPCM1024,JPCM1030,JPCM1032,JPCM1033,JPCM1034,JPCM1035,JPCM1037,JPCM1038,JPCM1039,JPCM1040,JPCM1041,JPCM1043,JPCM1045,JPCM1046,JPCM1047,JPCM1048,JPCM1050,JPCM1051,JPCM1052,JPCM1059,JPCM1060,JPCM1061,JPCM1062 -f true -i salmon_eq_classes $actaeus/raw_reads/ARO01_RNA*.out/aux_info/eq_classes.txt

# Run Lace 1.4 on the Corset results
#	select transcripts correctly classified by Corset
seqtk subseq $actaeus/ensemble_transcriptome/PerLocality/entap_outfiles/frame_selection/TransDecoder/processed/ORF-containing_transcripts_no_contam.fasta $actaeus/ensemble_transcriptome/PerLocality/corset/cluster.lace.list > $actaeus/ensemble_transcriptome/PerLocality/Lace/ORF-containing_transcripts_no_contam.ok.fasta
cd $actaeus/ensemble_transcriptome/PerLocality/Lace
$Software/necklace-1.01/tools/bin/lace ORF-containing_transcripts_no_contam.ok.fasta $actaeus/ensemble_transcriptome/PerLocality/corset/clusters.spaces.txt --outputDir $actaeus/ensemble_transcriptome/PerLocality/Lace --cores 24 --alternate --maxTran 100



