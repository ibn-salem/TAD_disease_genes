#!/usr/bin/Rscript
#=======================================================================
#
#   Quanitfication of all TSS contacts in Hi-C data
#
#=======================================================================

require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(BiocParallel)   # for parallel computing
require(chromint) 	# devtools::install_github("ibn-salem/chromint") to parse Hi-C from Rao et al. see: https://github.com/ibn-salem/chromint
require(InteractionSet)
require(tidyverse)


# TSS_FILE <- "data/geneTSSFromMGT_v1.0.txt"
TSS_FILE <- "data/geneTSSmiddleEndFromMGT_v2.0.txt"

#-----------------------------------------------------------------------
# Rao et al. 2014 Hi-C sample
#-----------------------------------------------------------------------
#~ CELL <- "K562"
CELL <- "GM12878_combined"
HIC_RESOLUTION <- 5*10^3 # 5kb
#HIC_RESOLUTION <- 5*10^4 # 50kb
#HIC_RESOLUTION <- 1*10^3 # 1kb
USE_LOCAL_HIC <- TRUE


#~ HIC_DATA_DIR_MURO="/project/jgu-cbdm/andradeLab/download/flat/databases/uncompressed/muro/ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GM12878_combined"
#~ HIC_DATA_DIR_MURO="/project/jgu-cbdm/andradeLab/download/flat/databases/uncompressed/muro/ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl"

HIC_DATA_DIR <- "/project/jgu-cbdm/ibnsalem/translocations/data/Rao2014"

RANDOM_SEED=13521


VERSION = "v03"

outPrefix = paste0("results/", VERSION, "/tss_contacts.", VERSION)

# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)


#-----------------------------------------------------------------------
# Options for parallel computation
#-----------------------------------------------------------------------

# use all available cores but generate random number streams on each worker
multicorParam <- MulticoreParam(RNGseed = RANDOM_SEED)
# set options
register(multicorParam)
# bpparam() # to print current options


txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqInfo <- seqinfo(txdb_hg19)

#-----------------------------------------------------------------------
# Load Hi-C data from Rao et al. 2014
#-----------------------------------------------------------------------
if ( !USE_LOCAL_HIC) {

	# parse normalized Hi-C map from Rao et al.
#~ 	gi = parseRaoHiCtoGI(CELL, HIC_RESOLUTION, HIC_DATA_DIR_MURO, seqInfo)
	# gi = parseRaoHiCtoGI(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo)
	gi = parseRaoHiCtoGI(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo, interChromosomal = FALSE)

	# annotate GI with genomic distance:
	gi$dist <- pairdist(gi)
	gi$cis <- intrachr(gi)


	# save data for faster loading next time
	save(gi, file = paste0(outPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".gi.RData"))

}else{
	load(paste0(outPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".gi.RData"))
}


#-----------------------------------------------------------------------
# Parse TSS coordintates
#-----------------------------------------------------------------------

regData <- read_tsv(
	TSS_FILE,
	col_types = cols(
	geneId = col_character(),
	finalLocus = col_character(),
	chr = col_character(),
	TSS = col_integer(),
	middle = col_integer(),
	end3 = col_integer()
	)
)

coordCols = c("TSS", "middle", "end3")

for (coordCol in coordCols) {
  
  outPrefixCoord = paste0(outPrefix, ".", coordCol)
  regionGR <- GRanges(
    regData$chr,
    IRanges(regData[[coordCol]], regData[[coordCol]]),
    seqinfo = seqInfo
  )
  
  #-----------------------------------------------------------------------
  # get contacts
  #-----------------------------------------------------------------------
  
  giRegHits <- findOverlaps(gi, regionGR)

  giRegHitsDF <- as_tibble(as.data.frame(giRegHits)) %>% 
    transmute(
      id = queryHits,
      reg_id = subjectHits
    )
  
  giDF <- as_tibble(as.data.frame(mcols(gi))) %>%
  	mutate(id = 1:n())
  
  regContacts <- bind_cols(
    giRegHitsDF,
    select(giDF[giRegHitsDF$id, ], -id)
  )
  
  #-----------------------------------------------------------------------
  # save reg contacts
  #-----------------------------------------------------------------------
  save(regContacts, file = paste0(outPrefixCoord, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".regContacts.RData"))
  
  
  #-----------------------------------------------------------------------
  # calculate contacts by reg
  #-----------------------------------------------------------------------
  regAllContacts <- regContacts %>%
  	group_by(reg_id) %>%
  	summarize(
  		sum_raw = sum(raw),
  		mean_raw = mean(raw, na.rm=TRUE),
  		n = n()
  		)
  
  #-----------------------------------------------------------------------
  # save reg contacts
  #-----------------------------------------------------------------------
  save(regAllContacts, file = paste0(outPrefixCoord, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".regAllContacts.RData"))
  

  #-----------------------------------------------------------------------
  # Annotate regData with counts
  #-----------------------------------------------------------------------
  regData <- regData %>%
    mutate(reg_id = seq(nrow(regData))) %>%
    left_join(regAllContacts, by = "reg_id")
  
  #-----------------------------------------------------------------------
  # write regData to output file
  #-----------------------------------------------------------------------
  write_tsv(regData, path = paste0(outPrefixCoord, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".reg_contacts.tsv"))
  
  
  #-----------------------------------------------------------------------
  # Some basic analysis
  #-----------------------------------------------------------------------
  # percent without any contacts
  
  noContacts <- regData %>%
    summarise(
      total = n(),
      na = sum(is.na(n)),
      na_percent = na / total * 100,
      zero = sum(n == 0, na.rm = TRUE)
    )
  
  write_tsv(
    noContacts,
    path = paste0(outPrefixCoord,
                  ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".contact_stats.tsv"))
  
  # distribution of contacts per region
  p <- regData %>%
    ggplot(aes(x = n)) +
    geom_histogram() +
    theme_bw() +
    labs(x = "Interacting regions per query region")
  ggsave(paste0(outPrefixCoord,
         ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".interacting_regions_per_query_region.pdf"),
         w = 7, h = 7)
  
  p <- regData %>%
    ggplot(aes(x = sum_raw)) +
    geom_histogram() +
    theme_bw() +
    labs(x = "Contact freq per query region")
  ggsave(paste0(outPrefixCoord,
         ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".contact_freq_per_query_region.pdf"),
         w = 7, h = 7)
  
  p <- regData %>%
    ggplot(aes(x = mean_raw)) +
    geom_histogram() +
    theme_bw() +
    labs(x = "Counts per interaction")
  ggsave(paste0(outPrefixCoord,
                   ".Hi-C.", CELL, ".", HIC_RESOLUTION,
                   ".counts_per_interaction.pdf"),
         w = 7, h = 7)

}

