## script for performing PCA by individuals on the full dataset
## Amanda Xuereb, Sept. 2022

#rm(list = ls())

library(adegenet)
library(ade4)
library(factoextra)
library(readr)
library(dplyr)
library(ggplot2)

## FIRST CONVERT VCF TO GENIND 

load("/path/to/file/01-all_coho_inds_filtered_GENIND.RData")

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

save(pca1, file = "/path/to/folder/pca_allSites.RData")


## PLOTTING 

## add pop and cu info to genind 
inds.cus <- read_delim("/path/to/folder/01-pop_map_allSites_CU.txt", delim = "\t", col_names = T)
head(inds.cus)
nrow(inds.cus)


## set strata in genind to CUs

strata(my_genind) <- inds.cus
setPop(my_genind) <- ~ CU 
pop(my_genind)

s.class(pca1$li, pop(my_genind))
fviz_eig(pca1, geom = "bar", addlabels = T, main = "proportion of variance explained")

percent= as.data.frame(pca1$eig/sum(pca1$eig)*100)
head(percent)
colnames(percent) <- "Proportion of Variance"
percent$Dimension <- seq(1,nrow(percent),1)

write.table(percent, "/path/to/folder/01-prop_explained_allSites.txt", col.names = T, row.names = F, quote = F, sep = "\t")

percent <- read.table("/path/to/folder/01-prop_explained_allSites.txt", header = T, stringsAsFactors = F)


### Colour points based on CU
pop(my_genind)
genind.pops <- as.data.frame(pop(my_genind))
head(genind.pops)
colnames(genind.pops) <- "Pop"
class(genind.pops$Pop)
genind.pops$CU <- as.character(genind.pops$Pop)
unique(genind.pops$Pop)

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

ggsave("/path/to/folder/pca12_filtered_allSites.pdf", width=10, height=10, dpi=600, units="in", useDingbats=F)
ggsave("/path/to/folder/pca12_filtered_allSites.png", width=10, height=11, dpi=600, units="in")
