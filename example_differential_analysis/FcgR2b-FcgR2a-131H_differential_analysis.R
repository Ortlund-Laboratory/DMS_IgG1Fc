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
antibodya <- "FcgR2a-131H"
ab_namea <- antibodya

inputdataa <- read_csv("FcgR2a-131H_mutation_to_Ka.csv")

fun_firsta <- function(x) {
  substring(x, 1,1)
}

fun_lasta <- function(x){
  substr(x, nchar(x), nchar(x))
}

fun_sitea <- function(x){
  substr(x, 2, nchar(x)-1)
}

wildtypea <- data.frame(sapply(inputdataa[1], fun_firsta))
#colnames(wildtype) <- c("wildtype")
inputdataa$wildtypea <- wildtypea[[1]]

mutation <- data.frame(sapply(inputdataa[1], fun_lasta))
colnames(mutation) <- c("mutation")
inputdataa$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(inputdataa[1], fun_sitea))
colnames(site) <- c("site")
inputdataa$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
inputdataa <- inputdataa[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
inputdataa <- inputdataa[
  with(inputdataa, order(site, mutation)),
]

inputdataa <- inputdataa[which(inputdataa$mutation != "*"),]

###############################################################################
###Add missing residues
###############################################################################
Fca <- read.fasta("Fc_prot.fasta")

###First reference data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_rangea <- min(inputdataa$site):max(inputdataa$site)
missing_aaa <- seq_rangea[!seq_rangea %in% unique(inputdataa$site)]
#2. Add the missing residue(s)
missing_dataa <- data.frame(site = missing_aaa, 
                           aa_substitutions = paste(toupper(Fca[["Fc_prot"]][missing_aaa]), missing_aaa, toupper(Fca[["Fc_prot"]][missing_aaa])),
                           mutation = toupper(Fca[["Fc_prot"]][missing_aaa]), 
                           wildtypea = toupper(Fca[["Fc_prot"]][missing_aaa]))
complete_inputa <- rbind.fill(inputdataa, missing_dataa)

#Expand the dataset to include NA values for synonymous amino acids
alla <- complete_inputa %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
inputdataa = complete_inputa %>% group_by_at(vars(wildtypea, site, mutation)) %>% 
  right_join(alla)

inputdataa <- inputdataa[
  with(inputdataa, order(site, mutation)),
]

################################################################################
### Generate a matrix for plotting sequence logos
### This matrix is also used to calculate average Ka scores
################################################################################
#remove the extra rwos that only contain the site number and nothing else 
inputdatafulla <- inputdataa
inputdataa <- inputdataa[which(!is.na(inputdataa$aa_substitutions)),]
logo_matrixa <- matrix(ncol = nrow(inputdatafulla)/20,nrow=20)
row.names(logo_matrixa) <- inputdatafulla$mutation[1:20]
colnames(logo_matrixa) <- seq(216,447, 1)


for(i in 0:ncol(logo_matrixa)-1){
  logo_matrixa[,i+1] <- inputdatafulla$logKa[seq(from = 20*i+1, to=20*i+20, by = 1)]
}

################################################################################
################################################################################
#NOW FOR SECOND SET OF DATA
################################################################################
################################################################################

antibodyb <- "FcgR2b"
ab_nameb <- antibodyb

inputdatab <- read_csv("FcgR2b_mutation_to_Ka.csv")

fun_firstb <- function(x) {
  substring(x, 1,1)
}

fun_lastb <- function(x){
  substr(x, nchar(x), nchar(x))
}

fun_siteb <- function(x){
  substr(x, 2, nchar(x)-1)
}

wildtypeb <- data.frame(sapply(inputdatab[1], fun_firstb))
#colnames(wildtype) <- c("wildtype")
inputdatab$wildtypeb <- wildtypeb[[1]]

mutation <- data.frame(sapply(inputdatab[1], fun_lastb))
colnames(mutation) <- c("mutation")
inputdatab$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(inputdatab[1], fun_siteb))
colnames(site) <- c("site")
inputdatab$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
inputdatab <- inputdatab[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
inputdatab <- inputdatab[
  with(inputdatab, order(site, mutation)),
]

inputdatab <- inputdatab[which(inputdatab$mutation != "*"),]

###############################################################################
###Add missing residues
###############################################################################
Fcb <- read.fasta("Fc_prot.fasta")

###First reference data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_rangeb <- min(inputdatab$site):max(inputdatab$site)
missing_aab <- seq_rangeb[!seq_rangeb %in% unique(inputdatab$site)]
#2. Add the missing residue(s)
missing_datab <- data.frame(site = missing_aab, 
                            aa_substitutions = paste(toupper(Fcb[["Fc_prot"]][missing_aab]), missing_aab, toupper(Fcb[["Fc_prot"]][missing_aab])),
                            mutation = toupper(Fcb[["Fc_prot"]][missing_aab]), 
                            wildtypeb = toupper(Fcb[["Fc_prot"]][missing_aab]))
complete_inputb <- rbind.fill(inputdatab, missing_datab)

#Expand the dataset to include NA values for synonymous amino acids
allb <- complete_inputb %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
inputdatab = complete_inputb %>% group_by_at(vars(wildtypeb, site, mutation)) %>% 
  right_join(allb)

inputdatab <- inputdatab[
  with(inputdatab, order(site, mutation)),
]


################################################################################
### Generate a matrix for plotting sequence logos
### This matrix is also used to calculate average Ka scores
################################################################################
#remove the extra rwos that only contain the site number and nothing else 
inputdatafullb <- inputdatab
inputdatab <- inputdatab[which(!is.na(inputdatab$aa_substitutions)),]
logo_matrixb <- matrix(ncol = nrow(inputdatafullb)/20,nrow=20)
row.names(logo_matrixb) <- inputdatafullb$mutation[1:20]
colnames(logo_matrixb) <- seq(216,447, 1)


for(i in 0:ncol(logo_matrixb)-1){
  logo_matrixb[,i+1] <- inputdatafullb$logKa[seq(from = 20*i+1, to=20*i+20, by = 1)]
}

################################################################################
################################################################################
#Subtract second matrix from the first
################################################################################
################################################################################

antibodyc="FcgR2b-FcgR2a-131H"
inputdatafullc=inputdatafullb
inputdatafullc$logKa=inputdatafullb$logKa-inputdatafulla$logKa
#inputdatac$p_adjb=abs(inputdatab$p_adjb-Ka_trimmeda$p_adja)
logo_matrixc=logo_matrixb-logo_matrixa
ab_namec <- antibodyc


################################################################################
##Plot heat maps
################################################################################

#Generate a dataframe containing the wild type amino acids 
#This will be used to mark wild type in the tiled heat map
tmp <- data.frame(sapply(inputdatafullc[1], fun_firsta), 
                  sapply(inputdatafullc[1], fun_sitea), 
                  sapply(inputdatafullc[1], fun_firsta))
colnames(tmp) <- c("wildtype","site", "mutation")
frames <- distinct(tmp)
frames$site <- as.numeric(frames$site)
inputdatafullc$site <- inputdatafullc$site+215
frames$site <- frames$site+215

#Change data type to integer for "site"
#Required for proper heatmap plotting
frames$site <- as.integer(frames$site)
inputdatafullc$site <- as.integer(inputdatafullc$site)

#Set the order for amino acids in the heatmap (by aa property)
polar <- c("H", "C", "S", "T", "N", "Q")
nonpolar <- c("G", "A", "V", "L", "I", "M", "P")
aromatic <- c("F", "Y", "W")
positive <- c("K", "R")
negative <- c("D", "E")

aa_order <- c(negative,positive, polar, nonpolar, aromatic)

#HeatMap
start = 0
end = 9000
mut_range <- subset(inputdatafullc, site>start & site<end)
mut_range$site <- as.factor(mut_range$site)
frames_range <- subset(frames, site>start & site<end)
frames_range$site <- as.factor(frames_range$site)

maximum <- max(inputdatafullc$logKa,na.rm = NA)
minimum <- min(inputdatafullc$logKa,na.rm = NA)

ggsave(filename = paste(antibodyc,"_LogKa_Heatmap01.png", sep=""), 
       ggplot(mut_range, aes(site, mutation)) + 
         geom_tile(color = "white",
                   fill = '#f7fbff', #azure1 for blueish grey background
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
         geom_point(aes(colour = logKa), alpha=1) +
         scale_colour_distiller(palette = "PiYG", direction = +1,
                                na.value = '#f7fbff',
                                limits=c(min(inputdatafullc$logKa, na.rm=T),max(inputdatafullc$logKa, na.rm=T))) +
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

############
# ADK44 need to normalize with respect to wildtype Kas.
FcgR2a131H <- 6.72974189580845
FcgR2b <-5.85397676784952
diff <- FcgR2b-FcgR2a131H
inputdatafullc$logKa_norm <- inputdatafullc$logKa-diff
maximum <- max(inputdatafullc$logKa_norm,na.rm = NA)
minimum <- min(inputdatafullc$logKa_norm,na.rm = NA)

#Individual Kas:
inputdatafullc$name <- paste(inputdatafullc$wildtypeb, inputdatafullc$site, inputdatafullc$mutation, sep = "")
write.csv(inputdatafullc, file = paste(antibodyc,"_Ka_fractions_na_included.csv", sep = ""), row.names = FALSE)
Ka_fractions_out <- na.omit(inputdatafullc[,c(2,3,4, ncol(inputdatafullc)-1, ncol(inputdatafullc))])
Ka_fractions_out$name <- paste(Ka_fractions_out$wildtypeb, Ka_fractions_out$site, Ka_fractions_out$mutation, sep = "")
write.csv(Ka_fractions_out, file = paste(antibodyc,"_Ka_fractions.csv", sep = ""), row.names = FALSE)

logo_matrixnorm <- logo_matrixc-diff

#ADK we want an average here, NOT colSums, to avoid the case where
#there may be fewer than 19 mutations in a column
average_Ka <- as.data.frame(colMeans(logo_matrixnorm, na.rm=T))
average_Ka$site <- seq(216,447,1)
colnames(average_Ka) <- c("Ka", "site")

#average Kas:
average_Ka_out <- as.data.frame(average_Ka$site)
average_Ka_out$'average Ka Score' <- average_Ka$Ka
colnames(average_Ka_out) <- c("Site", "average Ka Score")
write.csv(average_Ka_out, file = paste(antibodyc,"_average_Ka.csv", sep=""), row.names = FALSE)

#Set limits to 10th highest logKa_norm value rather than highest.
#Therefore, 1.00 and -1.00
inputdatafullc$logKa_norm_2 <- inputdatafullc$logKa_norm
inputdatafullc$logKa_norm_2 <- ifelse(inputdatafullc$logKa_norm_2>1.00,1.00,inputdatafullc$logKa_norm_2)
inputdatafullc$logKa_norm_2 <- ifelse(inputdatafullc$logKa_norm_2< -1.00,-1.00,inputdatafullc$logKa_norm_2)

ggsave(filename = paste(antibodyc,"_Normalized_LogKa_Heatmap01.png", sep=""), 
       ggplot(mut_range, aes(site, mutation)) + 
         geom_tile(color = "white",
                   fill = '#f7fbff', #azure1 for blueish grey background
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
         geom_point(aes(colour = inputdatafullc$logKa_norm_2), alpha=1) + 
         labs(color='logKa') +
         scale_colour_distiller(palette = "PiYG", direction = +1,
                                na.value = '#f7fbff',
                                limits=c(-1.00,1.00)) +
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
### Plot average Ka histogram
################################################################################

ggsave(filename = paste(antibodyc,"_average_Ka.png", sep=""), 
       ggplot(average_Ka, aes(site, Ka)) + 
         geom_bar(stat = "identity") + 
         xlim(215,448) +
         ylim(min(inputdatafullc$logKa-5), max(inputdatafullc$logKa+5)) +
         theme_minimal() + 
         theme_void(),
       width = 12, height = .4, dpi = 300, units = "in", device='png')

################################################################################
#Pymol file for mapping average scores onto the structure
###############################################################################
#Kaes
m_ <- mean(average_Ka$Ka, na.rm=T)
sd_ <- sd(average_Ka$Ka, na.rm=T)
cut_1 <- (m_ - 1*sd_)
cut_2 <- (m_ + 1*sd_)
Kaes_C01 <- subset(average_Ka,Ka < cut_1)
Kaes_C02 <- subset(average_Ka,Ka > cut_2)

#average_Ka$Ka <- ifelse(average_Ka$Ka < 0, 0, average_Ka$Ka)
average_Ka$level <- ifelse(average_Ka$Ka > cut_1 & average_Ka$Ka < cut_2, "grey", 
                             ifelse(average_Ka$Ka > cut_2, "red_1", "red_2" ))
average_Ka$level <- as.factor(average_Ka$level)

ggsave(filename = paste(antibodyc, "_color_average_Ka.png", sep=""), 
       ggplot(average_Ka, aes(site, Ka, fill = level)) + 
         geom_bar(stat="identity") + 
         #geom_segment(aes(x=215,xend=448,y=cut_1, yend = cut_1), linetype = 3, size = 0.15) +
         geom_segment(aes(x=215,xend=448,y=0, yend =0), linetype = 3, size = 0.15) +
         #geom_segment(aes(x=215,xend=448,y=m_, yend = m_), linetype = 1, size = 0.03, alpha = 0.5) +
         scale_fill_manual(values = c("grey" = "#d9d9d9",
                                      "red_1" ="#4dac26",
                                      "red_2"="#d01c8b")) +
         xlim(215,448) +
         theme_void() +
         theme(legend.position = "none"),
       width = 12, height = 0.4, dpi = 300, units = "in", device='png')


################################################################################
### Write a script for coloring residues in pymol
################################################################################
#Add "resi " in front of epitope residue names
Kaes_C01$resi <- paste("or resi", Kaes_C01$site, sep = " ")
reds01 <- paste(Kaes_C01$resi)
Kaes_C02$resi <- paste("or resi", Kaes_C02$site, sep = " ")
reds02 <- paste(Kaes_C02$resi)

sink(paste(ab_namec, "_pymol.txt", sep = ""))
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
cat("set_color red_1, [77,172,38]")
cat("\n")
cat("set_color red_2, [208,28,139]")
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


