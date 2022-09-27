rm(list = ls())

library(adegenet)
library(ade4)
library(factoextra)
library(readr)
library(dplyr)
library(ggplot2)

## FIRST NEED TO CREATE A GENIND OBJECT - SEE script xx-convert_vcf_to_genind.R

load("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/00-vcfs_and_other_data/01-all_coho_inds_filtered_24542snps_GENIND.RData")

na.tot <- sum(is.na(my_genind$tab))
dim(my_genind$tab)
full.dat <- dim(my_genind$tab)[1]*dim(my_genind$tab)[2]
na.tot/full.dat #proportion missing overall 

Y <- scaleGen(my_genind, NA.method = "mean")
class(Y)
Y[1:10,1:10]

pca1 <- dudi.pca(Y, cent = FALSE, scale = FALSE)
pca1
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

save(pca1, file = "/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/02-PCA_full/pca_allSites_24542Snps.RData")


## PLOTTING 
load("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/02-PCA_full/pca_allSites_24542Snps.RData")
load("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/00-vcfs_and_other_data/01-all_coho_inds_filtered_24542snps_GENIND.RData")

## add pop and cu info to genind 
inds.cus <- read_delim("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/00-population_maps/01-pop_map_allSites_CU.txt", delim = "\t", col_names = T)
head(inds.cus)
filter(inds.cus, POP == "MGG")
nrow(inds.cus)


## now set strata in genind to CUs

strata(my_genind) <- inds.cus
setPop(my_genind) <- ~ CU #~CU # or POP
pop(my_genind)


s.class(pca1$li, pop(my_genind))
fviz_eig(pca1, geom = "bar", addlabels = T, main = "proportion of variance explained")

percent= as.data.frame(pca1$eig/sum(pca1$eig)*100)
head(percent)
colnames(percent) <- "Proportion of Variance"
percent$Dimension <- seq(1,nrow(percent),1)
write.table(percent, "/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/02-PCA_full/01-prop_explained_allSites_24542snps.txt", col.names = T, row.names = F, quote = F, sep = "\t")

percent <- read.table("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/02-PCA_full/01-prop_explained_allSites_24542snps.txt", header = T, stringsAsFactors = F)


# # color by Thompson and BC regions
# 
# pca1$li[1:10,1:10]
# pop(my_genind1)
# 
# ## Colour based on Thompson vs. BC
# THO.pops <- read.table("/Volumes/Storage/epic4/00-DEC_2020/01-pop_map_THO.txt", header = T)
# nrow(THO.pops)
# length(unique(THO.pops$STRATA)) # 28 pops
# 
# class(pop(my_genind1))
# genind.pops <- as.data.frame(pop(my_genind1))
# head(genind.pops)
# colnames(genind.pops) <- "POP"
# genind.pops$COLOR <- "color"
# 
# library(dplyr)
# library(adegraphics)
# 
# filter(genind.pops, POP %in% THO.pops$STRATA)
# genind.pops[genind.pops$POP %in% THO.pops$STRATA,]$COLOR <- "Thompson"
# 
# nrow(filter(genind.pops, COLOR == "Thompson"))
# 
# genind.pops[!genind.pops$POP %in% THO.pops$STRATA,]$COLOR <- "Not Thompson"
# nrow(filter(genind.pops, COLOR == "Not Thompson"))
# 
# 
# p1 <- s.class(pca1$li, fac = pop(my_genind1), ppoints.cex = 0, plabels.cex = 0.8, pgrid.text.cex = 0)
# p2 <- s.class(pca1$li, fac = as.factor(genind.pops$COLOR), col = c("purple", "orange"), ellipseSize = 0, starSize = 0, plabels.cex = 0, pgrid.text.cex = 0)
# 
# superpose(p2, p1, plot = TRUE)
# 
# ade4::s.class(pca1$li, pop(my_genind1))
# ade4::s.class(pca1$li, as.factor(genind.pops$COLOR), col = transp(c("purple", "orange", 0.6)), cstar = 0, cellipse = 0, grid = F, add = T)


### Colour points based on CU
pop(my_genind)
genind.pops <- as.data.frame(pop(my_genind))
head(genind.pops)
colnames(genind.pops) <- "Pop"
class(genind.pops$Pop)
genind.pops$CU <- as.character(genind.pops$Pop)
unique(genind.pops$Pop)

#head(popmap.cu.color)

ind.coords <- as.data.frame(pca1$li)
head(ind.coords)

ind.coords$Ind <- rownames(ind.coords)
ind.coords$Pop <- substr(ind.coords$Ind, 1, 3)
ind.coords$CU <- genind.pops$CU
unique(ind.coords$CU)

centroid <- aggregate(cbind(Axis1, Axis2, Axis3, Axis4, Axis5, Axis6) ~ Pop, data = ind.coords, FUN = mean)
centroid

ind.coords <- left_join(ind.coords, centroid, by = "Pop", suffix = c("", ".cen"))
head(ind.coords)

#popmap.cu.color
inds.cus
pops.color <- select(inds.cus, POP, CU, color)
head(pops.color)
nrow(pops.color)

duplicated(pops.color$POP)
pops.color.nodup <- pops.color[-which(duplicated(pops.color$POP) == TRUE), ]
nrow(pops.color.nodup)
head(pops.color.nodup)

cols <- pops.color.nodup$color
cols
length(cols)


library(ggrepel)

g = ggplot(ind.coords, aes(x=Axis1, y=Axis2))+ 
  geom_hline(yintercept = 0, size = 0.5, color = 'lightgrey') +
  geom_vline(xintercept = 0, size = 0.5, color = "lightgrey") +
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, color = Pop), show.legend = F) +
  geom_point(aes(fill = Pop), shape = 21, size = 2.5, alpha = 0.8)+
  #geom_label_repel(data = centroid, aes(label = Pop, color = Pop, group = Pop), size = 2, label.padding = unit(0.05, "lines"), label.size = 0.2, show.legend = F, max.overlaps = 50, fontface = "bold", fill = "white", alpha = 0.8) +
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols) +
  labs(x= paste("PC1 (", round(percent$Proportion_of_Variance[1], 2), "%)", sep = ""))+
  labs(y= paste("PC2 (", round(percent$Proportion_of_Variance[2], 2), "%)", sep = ""))+
  theme_bw() +
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank(), legend.pos = "none")+
  theme(axis.text.y=element_text(colour="black",size=18))+
  theme(axis.text.x=element_text(colour="black",size=18))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title=element_text(size=20,colour="black",family="Helvetica"), 
        axis.text=element_text(size = 18, colour = "black", family = "Helvetica")) 
  
g


ggsave("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/02-PCA_full/pca12_filtered_allSites_24542snps_wlines_Sept2022.pdf", width=10, height=10, dpi=600, units="in", useDingbats=F)
ggsave("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/02-PCA_full/pca34_filtered_allSites_24542snps_wlines_May2022.png", width=10, height=11, dpi=600, units="in")



### PCA with allele frequencies 
load("/Volumes/Storage/epic4/00-DEC_2020/vcf/genind_allSites_allSnps.RData")
my_genind

# convert genind to genpop 
gp <- genind2genpop(my_genind)
gp

# get allele frequencies from genpop object
Xfreq <- makefreq(gp, missing = "mean")
class(Xfreq)
Xfreq[1:10,1:10]

# run pca on allele frequencies: 
pca.freq <- dudi.pca(Xfreq, scale = FALSE, scannf = FALSE)
s.label(pca.freq$li, sub = "Principal Components Analysis", csub = 1.2)
add.scatter.eig(pca.freq$eig, nf = 2, xax = 1, yax=2, posi = "top")

# make better plot: 
library(dplyr)
library(adegraphics)

THO.pops <- read.table("/Volumes/Storage/epic4/00-DEC_2020/vcf/01-pop_map_Thompson.txt", header = TRUE, stringsAsFactors = FALSE)
nrow(THO.pops)
head(THO.pops)
length(unique(THO.pops$STRATA)) # 27 pops

class(pop(my_genind))
genind.pops <- as.data.frame(pop(my_genind))
head(genind.pops)
length(unique(pop(my_genind))) # 145 in the full dataset

colnames(genind.pops) <- "POP"
genind.pops$COLOR <- "color"

filter(genind.pops, POP %in% THO.pops$STRATA)
genind.pops[genind.pops$POP %in% THO.pops$STRATA,]$COLOR <- "Thompson"
nrow(filter(genind.pops, COLOR == "Thompson"))

genind.pops[!genind.pops$POP %in% THO.pops$STRATA,]$COLOR <- "Not Thompson"
nrow(filter(genind.pops, COLOR == "Not Thompson"))

genind.pops <- genind.pops %>% distinct()
nrow(genind.pops)


p1 <- s.class(pca.freq$li, fac = unique(pop(my_genind)), ppoints.cex = 0, plabels.cex = 0.8, pgrid.text.cex = 0)
p2 <- s.class(pca.freq$li, fac = as.factor(genind.pops$COLOR), ppoints.fill = c("purple", "orange"), ppoints.cex = 2, ppoints.pch = 21, ppoints.col = "black", ellipseSize = 0, starSize = 0, plabels.cex = 0, pgrid.text.cex = 0)

superpose(p2, p1, plot = TRUE)

ade4::s.class(pca.freq$li, unique(pop(my_genind)))
ade4::s.class(pca.freq$li, as.factor(genind.pops$COLOR), col = transp(c("purple", "orange", 0.6)), cstar = 0, cellipse = 0, grid = F, add = T)

