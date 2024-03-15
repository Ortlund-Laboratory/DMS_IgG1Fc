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
antibody <- "FcgR2b"
ab_name <- antibody

inputdata <- read_csv("mutation_to_Ka.csv")

fun_first <- function(x) {
  substring(x, 1,1)
}

fun_last <- function(x){
  substr(x, nchar(x), nchar(x))
}

fun_site <- function(x){
  substr(x, 2, nchar(x)-1)
}

wildtype <- data.frame(sapply(inputdata[1], fun_first))
#colnames(wildtype) <- c("wildtype")
inputdata$wildtype <- wildtype[[1]]

mutation <- data.frame(sapply(inputdata[1], fun_last))
colnames(mutation) <- c("mutation")
inputdata$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(inputdata[1], fun_site))
colnames(site) <- c("site")
inputdata$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
inputdata <- inputdata[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
inputdata <- inputdata[
  with(inputdata, order(site, mutation)),
]

inputdata <- inputdata[which(inputdata$mutation != "*"),]

###############################################################################
###Add missing residues
###############################################################################
Fc <- read.fasta("Fc_prot.fasta")

###First reference data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_range <- min(inputdata$site):max(inputdata$site)
missing_aa <- seq_range[!seq_range %in% unique(inputdata$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           aa_substitutions = paste(toupper(Fc[["Fc_prot"]][missing_aa]), missing_aa, toupper(Fc[["Fc_prot"]][missing_aa])),
                           mutation = toupper(Fc[["Fc_prot"]][missing_aa]), 
                           wildtype = toupper(Fc[["Fc_prot"]][missing_aa]))
complete_input <- rbind.fill(inputdata, missing_data)

#Expand the dataset to include NA values for synonymous amino acids
all <- complete_input %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
inputdata = complete_input %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)

inputdata <- inputdata[
  with(inputdata, order(site, mutation)),
]

################################################################################
### Generate a matrix for plotting sequence logos
### This matrix is also used to calculate average escape scores
################################################################################
#remove the extra rwos that only contain the site number and nothing else 
inputdatafull <- inputdata
inputdata <- inputdata[which(!is.na(inputdata$aa_substitutions)),]

logo_matrix <- matrix(ncol = nrow(inputdata)/20,nrow=20)
row.names(logo_matrix) <- inputdata$mutation[1:20]
colnames(logo_matrix) <- seq(216,447, 1)


for(i in 0:ncol(logo_matrix)-1){
  logo_matrix[,i+1] <- complete_input$logKa[seq(from = 20*i+1, to=20*i+20, by = 1)]
}

################################################################################
##Plot heat maps
################################################################################

#Generate a dataframe containing the wild type amino acids 
#This will be used to mark wild type in the tiled heat map
tmp <- data.frame(sapply(inputdata[1], fun_first), 
                  sapply(inputdata[1], fun_site), 
                  sapply(inputdata[1], fun_first))
colnames(tmp) <- c("wildtype","site", "mutation")
frames <- distinct(tmp)
frames$site <- as.numeric(frames$site)
inputdata$site <- inputdata$site+215
frames$site <- frames$site+215

#Change data type to integer for "site"
#Required for proper heatmap plotting
frames$site <- as.integer(frames$site)
inputdata$site <- as.integer(inputdata$site)

write.table(inputdata,file="differential_matrix.txt")

#Set the order for amino acids in the heatmap (by aa property)
polar <- c("H", "C", "S", "T", "N", "Q")
nonpolar <- c("G", "A", "V", "L", "I", "M", "P")
aromatic <- c("F", "Y", "W")
positive <- c("K", "R")
negative <- c("D", "E")

aa_order <- c(negative,positive, polar, nonpolar, aromatic)
###ADK44 MUST CHANGE LIMITS 
#HeatMap
start = 0
end = 9000
mut_range <- subset(inputdata, site>start & site<end)
mut_range$site <- as.factor(mut_range$site)
frames_range <- subset(frames, site>start & site<end)
frames_range$site <- as.factor(frames_range$site)

ggsave(filename = "LogKa_Heatmap_no_limits.png", 
       ggplot(mut_range, aes(site, mutation)) + 
         geom_tile(color = "white",
                   fill = '#FCF0F0', #azure1 for blueish grey background
                   lwd = 0.1,
                   linetype = 1,
                   alpha = 1) +
         #scale_y_discrete(limits=rev) +
         ylim(rev(aa_order)) +
         geom_tile(data=frames_range, 
                   size=0,
                   height = 1,
                   fill='white', 
                   colour="white") +
         geom_point(aes(colour = logKa),alpha=1) +
         scale_colour_distiller(palette = "RdBu", direction = 1,
                                na.value = '#FCF0F0',
                            #    limits=c(min(inputdata$logKa, na.rm=T),max(inputdata$logKa, na.rm=T))) +
                                limits=c(4.77,6.93)) +
         scale_size(range = c(-0.2, 2)) +
         #The next section adds lightgrey tiles and black points for 
         #the wild type residues
         geom_point(inherit.aes = FALSE, 
                    data = frames_range,
                    aes(site, mutation),
                    shape=16,
                    size = 1.25,
                    colour="black") +
         xlab("Fc Site") + ylab("Mutation") +
         scale_x_discrete(breaks = levels(mut_range$site)[c(rep(F,4),T)]) +
         theme(# Hide panel borders and remove grid lines
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           # Change axis line
           #axis.line = element_line(colour = "black"),
           legend.key.size = unit(0.15, 'inch')
         ),
       width = 16, height = 3, dpi = 300, units = "in", device='png')

#Want a heatmap that pops. Therefore, normalize to the 10th lowest value and use
#that as bounds instead. In this case, the 10th lowest has a logKa of 5.41.
#Therefore, lower bound is 5.41. Calculate upper bound: 5.85-5.41=0.44, therefore:
#5.85+0.44=6.29

mut_range$logKa_edit <- mut_range$logKa
mut_range$logKa_edit <- ifelse(mut_range$logKa_edit>6.29,6.29,mut_range$logKa_edit)
mut_range$logKa_edit <- ifelse(mut_range$logKa_edit<5.41,5.41,mut_range$logKa_edit)

ggsave(filename = "LogKa_Heatmap_final_rep.png", 
       ggplot(mut_range, aes(site, mutation)) + 
         geom_tile(color = "white",
                   fill = '#FCF0F0', #azure1 for blueish grey background
                   lwd = 0.1,
                   linetype = 1,
                   alpha = 1) +
         #scale_y_discrete(limits=rev) +
         ylim(rev(aa_order)) +
         geom_tile(data=frames_range, 
                   size=0,
                   height = 1,
                   fill='white', 
                   colour="white") +
         geom_point(aes(colour = logKa_edit),alpha=1) +
         scale_colour_distiller(palette = "RdBu", direction = 1,
                                na.value = '#FCF0F0',
                                #    limits=c(min(inputdata$logKa, na.rm=T),max(inputdata$logKa, na.rm=T))) +
                                limits=c(5.41,6.29)) +
         scale_size(range = c(-0.2, 2)) +
         #The next section adds lightgrey tiles and black points for 
         #the wild type residues
         geom_point(inherit.aes = FALSE, 
                    data = frames_range,
                    aes(site, mutation),
                    shape=16,
                    size = 1.25,
                    colour="black") +
         xlab("Fc Site") + ylab("Mutation") +
         scale_x_discrete(breaks = levels(mut_range$site)[c(rep(F,4),T)]) +
         theme(# Hide panel borders and remove grid lines
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           # Change axis line
           #axis.line = element_line(colour = "black"),
           legend.key.size = unit(0.15, 'inch')
         ),
       width = 16, height = 3, dpi = 300, units = "in", device='png')


################################################################################
#Pymol file for mapping average scores onto the structure
################################################################################
complete_input$site <- complete_input$site+215
complete_input_stripped <- complete_input[c(2,5)]
average_scores <- complete_input_stripped %>% 
  group_by(site) %>% 
  summarise_at(vars(logKa), list(average = mean))

colnames(average_scores) <- c("site", "scores")

################# Need to change Ka of WT for each plot
#Escapees
m_ <- 5.85398 #WT score: 5.85398*19
#m_ <- mean(average_scores$scores, na.rm=T)
sd_ <- sd(average_scores$scores, na.rm=T)
cut_1 <- (m_ + 0.75*sd_)
cut_2 <- (m_ - 0.75*sd_)
escapees_C01 <- subset(average_scores,scores > cut_1)
escapees_C02 <- subset(average_scores,scores < cut_2)

average_scores$scores <- ifelse(average_scores$scores < 0, 0, average_scores$scores)
average_scores$level <- ifelse(average_scores$scores < cut_2, "red_2", 
                             ifelse(average_scores$scores > cut_1, "red_1", "grey" ))
average_scores$level <- as.factor(average_scores$level)

average_scores$scores <- average_scores$scores-5.85398
average_scores$scores[is.na(average_scores$scores)] <- 0
cut_1 <- cut_1-5.85398
cut_2 <- cut_2-5.85398
ggsave(filename = "LogKa_averagescores.png",
       ggplot(average_scores, aes(site, scores, fill = level)) + 
         geom_bar(stat="identity") + 
         geom_segment(aes(x=215,xend=448,y=cut_1, yend = cut_1), linetype = 3, size = 0.15) +
         geom_segment(aes(x=215,xend=448,y=cut_2, yend = cut_2), linetype = 3, size = 0.15) +
         #geom_segment(aes(x=215,xend=448,y=m_, yend = m_), linetype = 1, size = 0.03, alpha = 0.5) +
         scale_fill_manual(values = c("grey" = "#d9d9d9",
                                      "red_1" ="#2166ac",
                                      "red_2"="#b2182b")) +
         xlim(215,448) +
         theme_void() +
         theme(legend.position = "none"),
       width = 12, height = 0.4, dpi = 300, units = "in", device='png')



################################################################################
### Write a script for coloring residues in pymol
################################################################################
#Add "resi " in front of epitope residue names
escapees_C01$resi <- paste("or resi", escapees_C01$site, sep = " ")
reds01 <- paste(escapees_C01$resi)
escapees_C02$resi <- paste("or resi", escapees_C02$site, sep = " ")
reds02 <- paste(escapees_C02$resi)

sink(paste(ab_name, "_pymol.txt", sep = ""))
cat("#Set custom colors with RGB:")
cat("\n")
cat("set ray_shadow, 0")
cat("\n")
cat("#Set sphere radius:")
cat("\n")
cat("set sphere_scale, 0.5")
cat("\n")
cat("set ray_trace_mode, 1")
cat("\n")
cat("set_color red_1, [33,102,172]")
cat("\n")
cat("set_color red_2, [178,24,43]")
cat("\n")
cat("set_color grey_1, [217,217,217]")
cat("\n")
cat("set_color grey_2, [150,150,150]")
cat("\n")
cat("color grey_1, Fc_chA")
cat("\n")
cat("color grey_2, Fc_chB")
cat("\n")
cat("#Select sites to highlight with spheres:")
cat("\n")
cat("create spheresB1, Fc_chB and (resi 0 ", sep = "")
cat(paste(reds01, sep = ""), ")")
cat("\n")
cat("set sphere_color, red_1, spheresB1", sep = "")
cat("\n")
cat("hide everything, spheresB1", sep = "")
cat("\n")
cat("show spheres, spheresB1 and name CA", sep = "")
cat("\n")
cat("create spheresB2, Fc_chB and (resi 0 ", sep = "")
cat(paste(reds02, sep = ""), ")")
cat("\n")
cat("hide everything, spheresB2", sep = "")
cat("\n")
cat("show spheres, spheresB2 and name CA", sep = "")
cat("\n")
cat("set sphere_color, red_2, spheresB2", sep = "")
cat("\n")
cat("create surface_B, Fc_chB", sep = "")
cat("\n")
cat("select surf_B1, surface_B and (resi 0 ", sep = "")
cat(paste(reds01, sep = ""), ")")
cat("\n")
cat("select surf_B2, surface_B and (resi 0 ", sep = "")
cat(paste(reds02, sep = ""), ")")
cat("\n")
cat("color red_1, surf_B1", sep = "")
cat("\n")
cat("color red_2, surf_B2", sep = "")
cat("\n")
cat("show surface, surface_B", sep = "")
cat("\n")
#Now for the second molecule
cat("create spheresA1, Fc_chA and (resi 0 ", sep = "")
cat(paste(reds01, sep = ""), ")")
cat("\n")
cat("set sphere_color, red_1, spheresA1", sep = "")
cat("\n")
cat("hide everything, spheresA1", sep = "")
cat("\n")
cat("show spheres, spheresA1 and name CA", sep = "")
cat("\n")
cat("create spheresA2, Fc_chA and (resi 0 ", sep = "")
cat(paste(reds02, sep = ""), ")")
cat("\n")
cat("hide everything, spheresA2", sep = "")
cat("\n")
cat("show spheres, spheresA2 and name CA", sep = "")
cat("\n")
cat("set sphere_color, red_2, spheresA2", sep = "")
cat("\n")
cat("create surface_A, Fc_chA", sep = "")
cat("\n")
cat("select surf_A1, surface_A and (resi 0 ", sep = "")
cat(paste(reds01, sep = ""), ")")
cat("\n")
cat("select surf_A2, surface_A and (resi 0 ", sep = "")
cat(paste(reds02, sep = ""), ")")
cat("\n")
cat("color red_1, surf_A1", sep = "")
cat("\n")
cat("color red_2, surf_A2", sep = "")
cat("\n")
cat("show surface, surface_A", sep = "")
cat("\n")
sink()


