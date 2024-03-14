library(ggplot2)
#library(hrbrthemes)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(seqinr)
require(ggseqlogo)
library(tidyverse)
library(grid)
library(readr)

################################################################################
###Rename 
################################################################################

barcode_counts <- read_csv("trimmed_overall_binding_Kds.csv")

barcode_counts <- barcode_counts[c(13,5,6)]

barcode_counts$avgcountxlog10Ka <- barcode_counts$avgcount * barcode_counts$log10Ka

#####################################################################
#Combine duplicated mutations (same mutation with multiple barcodes)
#####################################################################
barcode_combined <- barcode_counts %>% 
  group_by(aa_substitutions) %>% 
  summarise_all(funs(sum))

barcode_combined$logKa <- barcode_combined$avgcountxlog10Ka / barcode_combined$avgcount

mutation_to_Ka <- barcode_combined[c(1,5)]
fun_site <- function(x){
  substr(x, 2, nchar(x)-1)
}
site <- data.frame(sapply(mutation_to_Ka[1], fun_site))
site$aa_substitutions <- as.numeric(site$aa_substitutions)
site$aa_substitutions <- site$aa_substitutions+215
mutation_to_Ka$site <- site$aa_substitutions

write.csv(mutation_to_Ka, file = "mutation_to_Ka.csv", row.names = FALSE)
