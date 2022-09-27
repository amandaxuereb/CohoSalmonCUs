### Script to perform DAPC and make scatterplots with ggplot
### run for each dataset (neutral, GEA candidates, and RDA candidates)
### Amanda Xuereb, Sept. 2022

rm(list = ls())
library(adegenet)
library(ade4)
library(dplyr)
library(readr)
library(stringr)
library(reshape2)
library(ggplot2)
library(vcfR)



## convert vcf to genind with populations as strata, then load genind

load("/path/to/folder/genind_BC_neutralSnps.RData")

my_genind
pop(my_genind)

dapc.pop <- dapc(my_genind, my_genind$pop) 

temp <- optim.a.score(dapc.pop, n.sim = 20)
dapc.pop2 <-dapc(my_genind, my_genind$pop, n.pca = temp$best)

save(dapc.pop2, file = "/path/to/folder/DAPC_BC_neutral_snps.RData")

## plot

rm(list = ls())

load("/path/to/folder/DAPC_BC_neutral_snps.RData")

library(readr)
new_pop_cu <- read_tsv("/path/to/folder/pop_cus_colours.txt")
nrow(new_pop_cu)
new_pop_cu <- filter(new_pop_cu, POP %in% names(dapc.pop2$prior))
head(new_pop_cu)
unique(new_pop_cu$color)

unique(dapc.pop2$grp)
new_pop_cu$POP

ind.coords <- as.data.frame(dapc.pop2$ind.coord)
head(ind.coords)
colnames(ind.coords) <- paste0("Axis", seq(1,ncol(ind.coords),1))

ind.coords$INDIVIDUALS <- rownames(ind.coords)
ind.coords$POP <- substr(ind.coords$INDIVIDUALS, 1, 3)
unique(ind.coords$POP)

centroid <- aggregate(cbind(Axis1, Axis2, Axis3, Axis4, Axis5, Axis6) ~ POP, data = ind.coords, FUN = mean)
centroid

ind.coords <- left_join(ind.coords, centroid, by = "POP", suffix = c("", ".cen"))
head(ind.coords)

ind.coords <- merge(ind.coords, new_pop_cu, by = "POP", sort = F)
centroid <- merge(centroid, new_pop_cu, by = "POP", sort = F)
centroid$color

percent= dapc.pop2$eig/sum(dapc.pop2$eig)*100

library(ggplot2)
library(ggrepel)

# plot
g = ggplot(ind.coords, aes(x=Axis1, y=Axis2, group = "POP"))+ 
  geom_hline(yintercept = 0, color = 'lightgrey') +
  geom_vline(xintercept = 0, color = 'lightgrey') +
  #geom_point(fill = ind.coords$color, col = ind.coords$color, shape = 21, size=1.5, show.legend = F, alpha = 0.3)+ # use this to plot individual points in addition to centroids
  geom_label_repel(data = centroid, aes(label = POP, group = POP), col = centroid$color, size = 2, label.padding = unit(0.1, "lines"), label.size = 0.7, show.legend = F, fontface = "bold", fill = "white", alpha = 0.8, segment.size = 0.3) +
  geom_point(aes(x = Axis1, y = Axis2), fill = centroid$color, color = 'black', shape = 23, size = 2.5, stroke = 0.7, data = centroid) +
  
  labs(x= paste("Axis 1 (", round(percent[1], 2), "%)", sep = ""))+
  labs(y= paste("Axis 2 (", round(percent[2], 2), "%)", sep = ""))+
  theme_bw() +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black",size=14))+
  theme(axis.text.x=element_text(colour="black",size=14))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=2),
        axis.title=element_text(size=16,colour="black",family="Helvetica",face="bold"), legend.position = "none", panel.grid = element_blank())
g

ggsave("/path/to/folder/DAPC_SCATTER_CENTROIDS_AX12_neutralSnps_BC.png", width = 15, height = 15, dpi = 600, units = 'cm')

