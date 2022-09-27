## Script for Neighbour-joining clustering
## run for each dataset (neutral, GEA candidates, and RDA candidates)
## Amanda Xuereb, Sept. 2022

rm(list = ls())

library(adegenet)
library(ape)
library(poppr)
library(dendextend)
library(readr)

# convert vcf files to genind with populations as strata
# create the NJ tree using the aboot function 

load("/path/to/folder/genind_BC_450RDAcandidates.RData")
my_genind

pop.map <- read.table("/path/to/folder/01-pop_map_BC_CU.txt", header = F, stringsAsFactors = F)
nrow(pop.map)
pop.map <- pop.map[,-3]
colnames(pop.map) <- c("STRATA","INDIVIDUALS")
length(unique(pop.map$STRATA))

new_pop_cu <- read_tsv("/path/to/folder/pop_cus_colours.txt")
nrow(new_pop_cu)
new_pop_cu <- filter(new_pop_cu, POP %in% unique(pop.map$STRATA))
new_pop_cu
nrow(new_pop_cu)

my_genpop <- genind2genpop(my_genind)

tree3 <- aboot(my_genpop, tree = "nj", distance = "edwards.dist", sample = 1000)
class(tree3)


sites <- "BC"
subset <- "RDAcandidates"
snpset <- "450snps" 

save(tree3, file = paste0("/path/to/folder/NJ_", sites, "_", subset, "_snpset", snpset, "_EDW_PHYLO.RData"))

load(paste0("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/10-cluster_NJ/NJ_", sites, "_", subset, "_snpset", snpset, "_EDW_PHYLO.RData"))

tree3$node.label
tree3$tip.label
tree3$edge

nodelab <- tree3$node.label
nodelab[which(nodelab < 80)] <- NA
nodelab

tiplab <- as.data.frame(tree3$tip.label, stringsAsFactors = F)
tiplab
colnames(tiplab) <- "POP"

new_tips <- merge(tiplab, new_pop_cu, by = "POP", sort = F)
nrow(new_tips)
new_tips$POP == tree3$tip.label

tree3$tip.label <- new_tips$pop_cu
tree3$tip.label

png(file = paste0("/path/to/folder/NJ_unrooted_", sites, "_", subset, "_snpset", snpset, "_EDW_aboot_nobg_rotated.png"), width = 7, height = 7, units = "in", res = 600)
plot.phylo(tree3, type = "unrooted", show.tip.label = TRUE, lab4ut = "axial", cex = 0.5, tip.color = c(new_tips$color), font = 2, label.offset = 0.002, align.tip.label = TRUE, rotate.tree = 160)
dev.off()

# tree with cluster support:
tree3.l <- ladderize(tree3)
tree3.l$tip.label <- new_tips$pop_cu

boots <- tree3.l$node.label
numsites <- length(unique(pop.map$STRATA))
boot.labels <- as.vector(which(boots >= 80) + numsites)
boot.labels

png(file = paste0("/path/to/folder/NJ_phylogram_", sites, "_", subset, "_snpset", snpset, "_EDW_aboot_nobg.png"), width = 8, height = 15, units = "in", res = 600)
plot.phylo(tree3.l, show.tip.label = TRUE, cex = 0.5, tip.color = c(new_tips$color), font = 2, label.offset = 0.001)
nodelabels(node = boot.labels, pch = 21, bg = 'green', cex = 0.6)
dev.off()


# correlation between neutral and adaptive trees 

load("/path/to/folder/NJ_BC_neutral_snpset22992snps_EDW_PHYLO.RData")

tree.neut <- tree3
class(tree.neut)

load("/path/to/folder/NJ_BC_GEAcandidates_snpset1117snps_EDW_PHYLO.RData")

tree.gea <- tree3
class(tree.gea)


library(dendextend)
library(ape)

cophenetic.neut <- cophenetic.phylo(tree.neut)
dim(cophenetic.neut)
cophenetic.neut[1:10,1:10]

cophenetic.gea <- cophenetic.phylo(tree.gea)
cophenetic.gea[1:10,1:10]

colnames(cophenetic.neut) == colnames(cophenetic.gea)
rownames(cophenetic.neut) == rownames(cophenetic.gea)

library(vegan)

mantel(cophenetic.neut, cophenetic.gea, method = 'pearson', permutations = 999)
