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
antibody <- "FcgR2b_Escape"
ab_name <- antibody

escape_counts <- read_tsv("escape_variant_counts.txt")
input_counts <- read_tsv("ref_variant_counts.txt")

escape_counts <- escape_counts[,2:3]
input_counts <- input_counts[,2:3]

################################################################################
### Histogram for read counts
################################################################################
ggsave(filename = paste(antibody, "_Number_of_Reads_escaped.png", sep=""), 
       ggplot(escape_counts, aes(x=count), title = "Mutations reads")+
         geom_histogram(color = "black", fill = "grey", ) +
         xlim(0,quantile(escape_counts$count, probs = .95, na.rm = T)) +
         # geom_vline(aes(xintercept=mu_bc_depth),
         #            linetype="dashed",
         #            color = "red") + 
         #geom_text(aes(label = mu), x = 3, y = 3, vjust = "inward", hjust = "inward") +
         ggtitle("Escape Read Counts") +
         xlab("# Reads") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

ggsave(filename = paste(antibody, "_Number_of_Reads_Myc.png", sep=""), 
       ggplot(input_counts, aes(x=count), title = "Mutations reads")+
         geom_histogram(color = "black", fill = "grey") +
         xlim(0,quantile(input_counts$count, probs = .99, na.rm = T)) +
         # geom_vline(aes(xintercept=mu_bc_depth),
         #            linetype="dashed",
         #            color = "red") + 
         #geom_text(aes(label = mu), x = 3, y = 3, vjust = "inward", hjust = "inward") +
         ggtitle("Reference Read Counts") +
         xlab("# Reads") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

#####################################################################
#Combine duplicated mutations (same mutation with multiple barcodes)
#####################################################################
escape_combined <- escape_counts %>% 
  group_by(mutation) %>% 
  summarise_all(funs(sum))

reference_combined <- input_counts %>% 
  group_by(mutation) %>% 
  summarise_all(funs(sum))

colnames(reference_combined) <- c("name", "depth")
colnames(escape_combined) <- c("name", "depth")

##Adjust all mutations in the reference set with less than the 5th percentile of reads
#depth_cutoff  <- quantile(reference_combined$depth, probs = .01, na.rm = T)
#reference_combined$depth <- ifelse(reference_combined$depth < depth_cutoff, depth_cutoff, reference_combined$depth)

################################################################################
###Convert mutation "names" into: wildtype, site, and mutations
################################################################################
####################################################################
##Renaming for escape
####################################################################

fun_first <- function(x) {
  substring(x, 1,1)
}

fun_last <- function(x){
  substr(x, nchar(x), nchar(x))
}

fun_site <- function(x){
  substr(x, 2, nchar(x)-1)
}

wildtype <- data.frame(sapply(escape_combined[1], fun_first))
#colnames(wildtype) <- c("wildtype")
escape_combined$wildtype <- wildtype[[1]]

mutation <- data.frame(sapply(escape_combined[1], fun_last))
colnames(mutation) <- c("mutation")
escape_combined$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(escape_combined[1], fun_site))
colnames(site) <- c("site")
escape_combined$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
escape_combined <- escape_combined[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
escape_combined <- escape_combined[
  with(escape_combined, order(site, mutation)),
]



####################################################################
##Renaming for reference
####################################################################
wildtype <- data.frame(sapply(reference_combined[1], fun_first))
#colnames(wildtype) <- c("wildtype")
reference_combined$wildtype <- wildtype[[1]]

mutation <- data.frame(sapply(reference_combined[1], fun_last))
colnames(mutation) <- c("mutation")
reference_combined$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(reference_combined[1], fun_site))
colnames(site) <- c("site")
reference_combined$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
reference_combined <- reference_combined[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
reference_combined <- reference_combined[
  with(reference_combined, order(site, mutation)),
]

#Clean-up
remove(mutation)
remove(site)
remove(wildtype)




###Remove stop codons
escape_trimmed <- escape_combined[which(escape_combined$mutation != "*"),]
reference_trimmed <- reference_combined[which(reference_combined$mutation != "*"),]


###############################################################################
###Add missing residues
###############################################################################
Fc <- read.fasta("Fc_prot.fasta")

###First reference data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_range <- min(reference_trimmed$site):max(reference_trimmed$site)
missing_aa <- seq_range[!seq_range %in% unique(reference_trimmed$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           name = paste(toupper(Fc[["Fc_prot"]][missing_aa]), missing_aa, toupper(Fc[["Fc_prot"]][missing_aa])),
                           mutation = toupper(Fc[["Fc_prot"]][missing_aa]), 
                           wildtype = toupper(Fc[["Fc_prot"]][missing_aa]))
complete_reference <- rbind.fill(reference_trimmed, missing_data)

#Expand the dataset to include NA values for synonymous amino acids
all <- complete_reference %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
reference_trimmed = complete_reference %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)


###Next escape data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_range <- min(escape_trimmed$site):max(escape_trimmed$site)
missing_aa <- seq_range[!seq_range %in% unique(escape_trimmed$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           name = paste(toupper(Fc[["Fc_prot"]][missing_aa]), missing_aa, toupper(Fc[["Fc_prot"]][missing_aa])),
                           mutation = toupper(Fc[["Fc_prot"]][missing_aa]), 
                           wildtype = toupper(Fc[["Fc_prot"]][missing_aa]))
complete_escape <- rbind.fill(escape_trimmed, missing_data)

#Expand the dataset to include NA values for synonymous amino acids
all <- complete_escape %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
escape_trimmed = complete_escape %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)



###Order reference_trimmed and escape_trimmed by site and mutation
escape_trimmed <- escape_trimmed[
  with(escape_trimmed, order(site, mutation)),
]
reference_trimmed <- reference_trimmed[
  with(reference_trimmed, order(site, mutation)),
]


###################################################################
### Calculate abundance: n/N
###################################################################
escape_trimmed$abundance <- escape_trimmed$depth/sum(escape_trimmed$depth, na.rm=TRUE)
reference_trimmed$abundance <- reference_trimmed$depth/sum(reference_trimmed$depth, na.rm=TRUE)

###Then adjust the lowest reference abundance values
#Use 95 percentile
#This helps remove exaggerated escape fractions caused by dividing by a very small number
cutoff  <- quantile(reference_trimmed$abundance, probs = .05, na.rm = T)
reference_trimmed$abundance <- ifelse(reference_trimmed$abundance<cutoff, cutoff, reference_trimmed$abundance)

##Fishers exact test
esc_fisher <- escape_trimmed[,1:5]
esc_fisher$ref_depth <- reference_trimmed$depth
esc_sum <- sum(na.omit(esc_fisher$depth))
ref_sum <- sum(na.omit(esc_fisher$ref_depth))

p_values <- data.frame(matrix(ncol = 1, nrow = nrow(esc_fisher)))
colnames(p_values) <- "p"
#Calculate p values for each comparison (using Fisher's exact test)
for (i in 1:nrow(esc_fisher)){
  if (!anyNA(esc_fisher[i,])){
    f_ <- matrix(unlist(c(esc_fisher[i,5], 
                          esc_fisher[i,6],
                          esc_sum-esc_fisher[i,5],
                          ref_sum-esc_fisher[i,6])),2)
    p_ <- fisher.test(f_)
    p_values$p[i] <-p_$p.value
  }
}

p_values$p_adj <- p.adjust(p_values$p, "bonferroni")
p_values$p_adj <- ifelse(p_values$p_adj == 0, 1e-300, p_values$p_adj)

###################################################################
### Calculate escape fractions
### And add adjusted p values
###################################################################
escape_trimmed$Escape <- (escape_trimmed$abundance/reference_trimmed$abundance)
escape_trimmed$p_adj <- p_values$p_adj
max(escape_trimmed$Escape, na.rm=T)
min(escape_trimmed$Escape, na.rm=T)
#Write out raw escape score before normalizing
# escape_raw_out <- na.omit(escape_trimmed[,c(1, ncol(escape_trimmed)-2)])
# colnames(escape_raw_out) <- c("Mutation", "EscapeScore")
# write.csv(escape_raw_out, file = paste(antibody,"_raw_escape.csv", sep = ""), row.names = FALSE)

escape_trimmed$Non_Normalized_Enrich <- escape_trimmed$Escape

###Normalize the data between the 99 and 1 percentiles
max(escape_trimmed$Escape, na.rm=T)
upper_limit  <- quantile(escape_trimmed$Escape, probs = .99, na.rm = T)
lower_limit <- min(escape_trimmed$Escape, na.rm=T)
escape_trimmed$Escape <- ifelse(escape_trimmed$Escape>upper_limit, upper_limit, escape_trimmed$Escape)
escape_trimmed$Escape <- (escape_trimmed$Escape-lower_limit)/(upper_limit-lower_limit)

ggsave(filename = paste(antibody,"_Histogram.png", sep=""), 
       ggplot(escape_trimmed, aes(x=Escape), title = "Escape Fractions")+
         #geom_histogram(aes(y=..density..), color = "black", fill = "grey", binwidth = 0.025) +
         geom_histogram(color = "black", fill = "grey") +
         xlim(0,quantile(escape_trimmed$Escape, probs = .99, na.rm = T)) +
         xlab("Escape Fraction") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

# ##Calculate destabilization instead of escape!!!
# mx <- max(escape_trimmed$Escape, na.rm=T)
# mean_dest <- median(escape_trimmed$Escape, na.rm=T)
# escape_trimmed$destabilization <-abs(mx-escape_trimmed$Escape)
#escape_trimmed <- na.omit(escape_trimmed)


################################################################################
### Generate a matrix for plotting sequence logos
### This matrix is also used to calculate total escape scores
################################################################################
#remove the extra rwos that only contain the site number and nothing else 
escape_trimmed <- escape_trimmed[which(!is.na(escape_trimmed$mutation)),]

logo_matrix <- matrix(ncol = nrow(escape_trimmed)/20,nrow=20)
row.names(logo_matrix) <- escape_trimmed$mutation[1:20]
colnames(logo_matrix) <- seq(216,447, 1)


for(i in 0:ncol(logo_matrix)-1){
  logo_matrix[,i+1] <- escape_trimmed$Escape[seq(from = 20*i+1, to=20*i+20, by = 1)]
}

################################################################################
##Plot heat maps
################################################################################

#Generate a dataframe containing the wild type amino acids 
#This will be used to mark wild type in the tiled heat map
tmp <- data.frame(sapply(escape_trimmed[1], fun_first), 
                  sapply(escape_trimmed[1], fun_site), 
                  sapply(escape_trimmed[1], fun_first))
colnames(tmp) <- c("wildtype","site", "mutation")
frames <- distinct(tmp)
frames$site <- as.numeric(frames$site)
escape_trimmed$site <- escape_trimmed$site+215
frames$site <- frames$site+215

#Change data type to integer for "site"
#Required for proper heatmap plotting
frames$site <- as.integer(frames$site)
escape_trimmed$site <- as.integer(escape_trimmed$site)

write.table(escape_trimmed,file="differential_matrix.txt")

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
mut_range <- subset(escape_trimmed, site>start & site<end)
mut_range$site <- as.factor(mut_range$site)
frames_range <- subset(frames, site>start & site<end)
frames_range$site <- as.factor(frames_range$site)

ggsave(filename = paste(antibody,"_Fraction_heatmap01.png", sep=""), 
       ggplot(mut_range, aes(site, mutation, size = -log(p_adj))) + 
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
         geom_point(aes(colour = Escape), alpha=1) +
         scale_colour_distiller(palette = "RdPu", direction = +1,
                                na.value = '#FCF0F0',
                                limits=c(0,max(escape_trimmed$Escape, na.rm=T))) +
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
###Write files:
################################################################################

################################################################################
### Plot Total Escape histogram
################################################################################
total_escape <- as.data.frame(colSums(logo_matrix, na.rm=T))
total_escape$site <- seq(216,447,1)
colnames(total_escape) <- c("escape", "site")

ggsave(filename = paste(antibody,"_Total_escape.png", sep=""), 
       ggplot(total_escape, aes(site, escape)) + 
         geom_bar(stat = "identity") + 
         xlim(215,448) +
         ylim(0, max(total_escape$escape+5)) +
         theme_minimal() + 
         theme_void(),
       width = 12, height = .4, dpi = 300, units = "in", device='png')


#Individual escapes:
escape_fractions_out <- na.omit(escape_trimmed[,c(2,3,4, ncol(escape_trimmed)-1, ncol(escape_trimmed))])
escape_fractions_out$name <- paste(escape_fractions_out$wildtype, escape_fractions_out$site, escape_fractions_out$mutation, sep = "")
write.csv(escape_fractions_out, file = paste(antibody,"_fractions.csv", sep = ""), row.names = FALSE)

#Total escapes:
total_escape_out <- as.data.frame(total_escape$site)
total_escape_out$'Total Escape Score' <- total_escape$escape
colnames(total_escape_out) <- c("Site", "Total Escape Score")
write.csv(total_escape_out, file = paste(antibody,"_total_escape.csv", sep=""), row.names = FALSE)

################################################################################
#Pymol file for mapping total scores onto the structure
################################################################################
#Escapees
m_ <- mean(total_escape$escape, na.rm=T)
sd_ <- sd(total_escape$escape, na.rm=T)
cut_1 <- (m_ + 1*sd_)
cut_2 <- (m_ + 2*sd_)
escapees_C01 <- subset(total_escape,escape > cut_1 & escape < cut_2)
escapees_C02 <- subset(total_escape,escape > cut_2)

total_escape$escape <- ifelse(total_escape$escape < 0, 0, total_escape$escape)
total_escape$level <- ifelse(total_escape$escape < cut_1, "grey", 
                             ifelse(total_escape$escape < cut_2, "red_1", "red_2" ))
total_escape$level <- as.factor(total_escape$level)

ggsave(filename = paste(antibody, "_color_total_escape.png", sep=""), 
       ggplot(total_escape, aes(site, escape, fill = level)) + 
         geom_bar(stat="identity") + 
         geom_segment(aes(x=215,xend=448,y=cut_1, yend = cut_1), linetype = 3, size = 0.15) +
         geom_segment(aes(x=215,xend=448,y=cut_2, yend = cut_2), linetype = 3, size = 0.15) +
         #geom_segment(aes(x=215,xend=448,y=m_, yend = m_), linetype = 1, size = 0.03, alpha = 0.5) +
         scale_fill_manual(values = c("grey" = "#d9d9d9",
                                             "red_1" ="#dd3497",
                                             "red_2"="#4a2267")) +
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
cat("set_color red_1, [221,52,151]")
cat("\n")
cat("set_color red_2, [73,0,106]")
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


