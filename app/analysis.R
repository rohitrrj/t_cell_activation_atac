# ### analysis.R ###

library(trackViewer)
library(ggplot2)
library(vroom)

File_location <- "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit"

# vstNormalizedCounts_t_cell_activation <- read.table(paste0(File_location,"/WebInputFiles/t_cell_activation/","vstNormalizedCounts_t_cell_activation_filt.txt"),header = T)
specification <- cols(
  SampleID = col_character(),
  Sample = col_character(),
  Hours = col_character(),
  Concentration = col_character(),
  Age = col_character(),
  Group = col_character(),
  GroupAge = col_character(),
  Peaks = col_character(),
  vstNormalizedCounts = col_double(),
  .delim = "\t"
)
vstNormalizedCounts_t_cell_activation <- vroom(paste0(File_location,"/WebInputFiles/t_cell_activation/","vstNormalizedCounts_t_cell_activation_filt.txt"),
                                               col_types = specification)
UCSC.hg19.genes<- read.table(paste0(File_location,"/WebInputFiles/commonFiles/","Gene_Symbols.txt"),header = F)

`%notin%` <- Negate(`%in%`)