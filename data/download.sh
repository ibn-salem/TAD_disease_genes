
#=======================================================================
# this scirpt is supposed to document data download and transforamtion
#=======================================================================


# tss coordintas from Enrique (mail from 2017-03-13)
# saved as file geneTSSFromMGT_v1.0.txt.bz2

# unzip file:
bzip2 -dk geneTSSFromMGT_v1.0.txt.bz2

# coordintas from Enrique (mail from 2017-04-28)
# saved as file geneTSSmiddleEndFromMGT_v2.0.txt.gz

gunzip geneTSSmiddleEndFromMGT_v2.0.txt.gz


#===============================================================================
# TADs by number of genes
#===============================================================================
# manuall copy file: "4Jonas_writeTable_dataFrameTADmodified_Rao_GM12878_v4.0.txt"
# from Enriques mail from 23.06.17

#===============================================================================
# Download ClinVar variants
#===============================================================================

mkdir -p ClinVar
wget -P ClinVar ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip ClinVar/variant_summary.txt.gz
