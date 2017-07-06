#'==============================================================================
#'
#' Analysis of structural variants from ClinVar at TAD boundaries by the number
#'  of genes in each TAD
#'
#'==============================================================================

require(tidyverse)
require(GenomicRanges)


#-------------------------------------------------------------------------------
# some parameters
#-------------------------------------------------------------------------------

CLINVAR_FILE <- "data/ClinVar/variant_summary.txt"
TAD_FILE <- "data/4Jonas_writeTable_dataFrameTADmodified_Rao_GM12878_v4.0.txt"

outPrefix <- "results/SVs_at_TADs"

# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# parse TAD dataset
#-------------------------------------------------------------------------------

tadDF <- read_tsv(TAD_FILE, guess_max = 10000, col_types = 
                    cols(
                      TADlabel.chr.TADlabel = col_character(),
                      chr = col_character(),
                      s = col_integer(),
                      e = col_integer(),
                      TADid = col_character(),
                      clusterLabel.chr.clLabel.numTADsInCl = col_character(),
                      Clchr = col_character(),
                      ClS = col_integer(),
                      ClE = col_integer(),
                      nG = col_integer(),
                      nGD = col_integer(),
                      nGND = col_integer()
                    ))

tadGR <- GRanges(
  paste0("chr", tadDF$chr),
  IRanges(tadDF$s, tadDF$e),
  strand = "*",
  select(tadDF, -(2:4))
)


boundaryDF <- tadDF[rep(1:nrow(tadDF), 2), ] %>% 
  mutate(boundaryID = 1:n())

boundaryGR <- GRanges(
  paste0("chr", boundaryDF$chr),
  IRanges(
    c(tadDF$s, tadDF$e), 
    c(tadDF$s, tadDF$e)),
  strand = "*",
  select(boundaryDF, -(2:4))
)

#-------------------------------------------------------------------------------
# parse ClinVar
#-------------------------------------------------------------------------------

clinvar_raw <- read_tsv(CLINVAR_FILE, col_types = cols(
  .default = col_character(),
  `#AlleleID` = col_integer(),
  GeneID = col_integer(),
  ClinSigSimple = col_integer(),
  `RS# (dbSNP)` = col_integer(),
  Start = col_integer(),
  Stop = col_integer(),
  NumberSubmitters = col_integer(),
  SubmitterCategories = col_integer()
))

clinvar <- clinvar_raw %>% 
  dplyr::rename(AllelID = `#AlleleID`) %>% 
  filter(Assembly == "GRCh37") %>% 
  mutate(Size = Stop - Start + 1)

# sv <- clinvar %>% 
#   filter(Type %in% c("deletion")) %>% 
#   filter(ClinicalSignificance %in% c("Pathogenic", "Benign")) %>% 
#   mutate(svID = 1:n()) 
# 

# iterate over size thresholods
SIZE_TH = c(0, 5000, 10000, 10^5, 10^6, 10^7)

for (SIZE in SIZE_TH) {
  
  outPrefixSize <- paste0(outPrefix, ".", SIZE)  
  
  sv <- clinvar %>% 
    filter(Type %in% c("deletion")) %>% 
    filter(ClinicalSignificance %in% c("Pathogenic", "Benign")) %>% 
    filter(Size >= SIZE) %>% 
    mutate(svID = 1:n()) 

  svGR <- GRanges(
    paste0("chr", sv$Chromosome),
    IRanges(sv$Start, sv$Stop),
    strand = "*",
    select(sv, Type, ClinicalSignificance, Size, AllelID, svID)
  )
  
  
  #===============================================================================
  # SVs at boundaries
  #===============================================================================
  
  boudnaryToSV <- as_tibble(as.data.frame(findOverlaps(boundaryGR, svGR))) %>% 
    dplyr::rename(boundaryID = queryHits, svID = subjectHits)
  
  boundaryHits <- boundaryDF %>% 
    left_join(boudnaryToSV, by = "boundaryID") %>% 
    left_join(sv, by = "svID") 
    
  boundaryCounts <- boundaryHits %>% 
    group_by(boundaryID, nG, ClinicalSignificance) %>%
    summarize(n_hits = sum(!is.na(svID)))
  
  # print(boundaryCounts %>% ungroup() %>% count(n_hits))
  
  nGeneCounts <- boundaryCounts %>% 
    group_by(nG, ClinicalSignificance) %>% 
    summarize(n = n()) %>% 
    mutate(
      N = sum(n),
      freq = n / N,
      percent = freq * 100)
  
  valueDF <- nGeneCounts %>% 
    filter(ClinicalSignificance == "Pathogenic")
  
  #-------------------------------------------------------------------------------
  # percent of boundaries with pathogenic deletion
  #-------------------------------------------------------------------------------
  
  p <- ggplot(valueDF, 
              aes(x = nG , y = percent)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label = paste0(n, "/", N, " (",signif(percent, 3), "%)" )),
              angle = 90, hjust = "left") +
    ylim(0, 1.3 * max(valueDF$percent)) + 
    theme_bw() + 
    labs(y = "TAD boundaries with pathogenic deletion [%]",
         x = "Number of genes in TAD")
  
  ggsave(p, file = paste0(outPrefixSize, ".pathogenic_deletions_at_TADboundaries.barplot.pdf"), w = 5, h = 5)
  
  #-------------------------------------------------------------------------------
  # Number of tad boundaries with n genes
  #-------------------------------------------------------------------------------
  
  countNG <- boundaryDF %>% 
    count(nG)
  
  p <- ggplot(countNG, aes(x = nG , y = n)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label = n), angle = 90, hjust = "left") +
    theme_bw() + 
    labs(y = "Boundaries",
         x = "Number of genes in TAD")
  
  ggsave(p, file = paste0(outPrefixSize, ".number_of_Bounaries_by_nG.barplot.pdf"), w = 6, h = 6)
  
}

#===============================================================================
# Analyse ClinVar
#===============================================================================

#-------------------------------------------------------------------------------
# Number of variant Types
#-------------------------------------------------------------------------------

p <- ggplot(count(clinvar, Type), aes(x = Type, y = n, fill = Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), angle = 90, vjust = "right") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p, file = paste0(outPrefix, ".clinvar.Type.barplot.pdf"))

#-------------------------------------------------------------------------------
# Number of ClinicalSignificance
#-------------------------------------------------------------------------------
p <- ggplot(count(clinvar, ClinicalSignificance), aes(x = ClinicalSignificance, y = n, fill = ClinicalSignificance)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), hjust = -0.1) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  coord_flip()

ggsave(p, file = paste0(outPrefix, ".clinvar.ClinicalSignificance.barplot.pdf"), h = 12, w = 6)

#-------------------------------------------------------------------------------
# SV number by type and significance
#-------------------------------------------------------------------------------

countDF <- sv %>% 
  count(Type, ClinicalSignificance)

write_tsv(countDF, paste0(outPrefix, ".sv_n_by_clinsig.tsv"))

#-------------------------------------------------------------------------------
# SV size by type
#-------------------------------------------------------------------------------

p <- ggplot(sv, aes(x = Size, fill = ClinicalSignificance, color = ClinicalSignificance)) +
  geom_histogram(alpha = 0.4, , position="identity") +
  facet_grid(ClinicalSignificance ~ . , scales = "free_y") +
  scale_x_log10() +
  theme_bw()
ggsave(p, file = paste0(outPrefix, ".sv_deletions.size_by_ClinicalSignificance.histogram.pdf"), h = 3, w = 6)



