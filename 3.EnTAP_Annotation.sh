# Annotate the transcriptome with EnTAP

# Index the databases with diamond
diamond makedb --in $Software/EnTAP/databases/refseq_vert_2022.fasta -d $Software/EnTAP/databases/refseq_vert.dmnd
diamond makedb --in $Software/EnTAP/databases/uniprot_sprot.fasta -d $Software/EnTAP/databases/uniprot_sprot.dmnd

# configure EnTAP
$Software/EnTAP/EnTAP --config -d $Software/EnTAP/databases/uniprot_sprot.dmnd -d $Software/EnTAP/databases/refseq_vert.dmnd

# run EnTAP
cd $Software/EnTAP
EnTAP --runP -i $actaeus/ensemble_transcriptome/PerLocality/evidentialgene/okayset/transcript_collection.okay.mrna \
-d /data/biolinux/GenBank/uniprot/uniprot_sprot.dmnd -d /data/biolinux/GenBank/refseq/refseq_vert.dmnd \
-t 24


