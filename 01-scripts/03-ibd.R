### plot for Isolation by Distance analysis; calculating Bst and Hs as a function of geographic distance
### Amanda Xuereb, Sept. 2022

library(readr)
library(hierfstat)
library(reshape2)
library(geosphere)

## get distances 
latlong <- read.table("/path/to/folder/coho_dist_v2.txt", header = T, stringsAsFactors = F)
head(latlong)
latlong$POP_ID

load("/path/to/folder/genind_noThompson.RData")
popmap <- read_tsv("/path/to/folder/01-pop_map_allSites_CU.txt")
strata(my_genind) <- popmap
setPop(my_genind) <- ~POP
pop(my_genind)

latlong <- filter(latlong, POP_ID %in% my_genind@pop)
nrow(latlong)
head(latlong)

filter(latlong, Latitude == min(latlong$Latitude))
# SOO is the southernmost site (for all Sites and coastal BC)
# COL is the southernmost site in Thompson 


coords <- latlong[,c(2,3)]
rownames(coords) <- latlong[,1]
head(coords)

dist.mat <- distm(coords, fun = distGeo)
rownames(dist.mat) <- rownames(coords)
colnames(dist.mat) <- rownames(dist.mat)
dist.mat[1:10, 1:10]
class(dist.mat)


dist.df <- melt(dist.mat)
dist.south <- filter(dist.df, Var1 == "SOO")
#dist.col <- filter(dist.df, Var1 == "COL")
head(dist.south)
nrow(dist.south)
colnames(dist.south) <- c("COL", "POP", "DIST")


hf <- genind2hierfstat(my_genind)
hf
hf[1:10,1:10]

hf.stat <- basic.stats(hf, diploid = TRUE)

class(hf.stat$Hs)

hs <- as.data.frame(hf.stat$Hs)
head(hs)

hs.mean <- as.data.frame(colMeans(hs))
hs.mean$POP <- rownames(hs.mean)
head(hs.mean)
colnames(hs.mean)[1] <- "HS"

hs.mean.dist <- merge(hs.mean, dist.south, by = "POP")
head(hs.mean.dist)

plot(hs.mean.dist$DIST/1000, hs.mean.dist$HS, xlab = "Distance from southernmost site (km)", ylab = "HS", pch = 21, bg = "darkgray")

summary(lm(hs.mean.dist$DIST ~ hs.mean.dist$HS))


h <- ggplot(data = hs.mean.dist,
            aes(x = DIST/1000, y = HS))
h <- h + stat_smooth(method = "lm")
h <- h + geom_point(pch = 21, fill = 'orange' , size = 2)
h <- h + theme_bw()
h <- h + labs(x = "Distance from southernmost site (km)", y = expression(H["S"]) )
h <- h + theme(axis.title.x = element_text(size = 12, family = "Helvetica",face = "bold"),
               axis.text.x = element_text(size = 12,family = "Helvetica",
                                          angle = 90, hjust = 0.5, vjust = 0.5, color = 'black'),
               axis.title.y = element_text(size = 12, family = "Helvetica",
                                           face = "bold",angle = 0, hjust = 0, vjust = 0.5),
               axis.text.y = element_text(size=12,
                                          family = "Helvetica", color = 'black'),
               strip.text.x = element_text(size=12))
h <- h + theme(panel.grid.major = element_blank() )
h <- h + annotate(geom = "text",x = 1350, y = 0.077,
                 label = "paste(italic(R) ^ 2, \" = .48, p < 2e-16\")",
                 parse = TRUE,
                 color = "black", size = 3)
h  <- h + theme(legend.position = "none")

h

ggsave("/path/to/folder/hs_vs_dist_from_southmost_lm_BC.png", width = 6, height = 5, units = "in")


bstat <- betas(hf, nboot = 100)

bstat.df <- as.data.frame(t(rbind(bstat$betaiovl,bstat$ci)))
class(bstat.df)
head(bstat.df)
colnames(bstat.df)[1] <- "bst"
bstat.df$POP <- rownames(bstat.df)

head(dist.south)
dist.south <- dist.south[,c(2,3)]
colnames(dist.south) <- c("POP", "DIST")
head(dist.south)

dist.bst.south <- merge(dist.south, bstat.df, by = "POP")
head(dist.bst.south)
nrow(dist.bst.south)

plot(dist.bst.south$DIST/1000, dist.bst.south$bst, ylab = "BST", xlab = "Distance from southernmost site (km)", pch = 21, bg = "darkgrey")

summary(lm(dist.bst.south$DIST ~ dist.bst.south$bst))

d <- ggplot(data = dist.bst.south,
            aes(x = DIST/1000, y = bst))
d <- d + stat_smooth(method = "lm")
d <- d + geom_point(pch = 21, fill = "orange", size = 2)
d <- d + theme_bw()
d <- d + labs(x = "Distance from southernmost site (km)", y = expression(beta["ST"]) ) 
d <- d + theme(axis.title.x = element_text(size = 12, family = "Helvetica", face = "bold"),
               axis.text.x = element_text(size = 12,family = "Helvetica",
                                          angle = 90, hjust = 0.5, vjust = 0.5, color = 'black'),
               axis.title.y = element_text(size = 12, family = "Helvetica",
                                           face = "bold",angle = 0, hjust = 0, vjust = 0.5),
               axis.text.y = element_text(size=12,
                                          family = "Helvetica", color = 'black'),
               strip.text.x = element_text(size=12))
d <- d + theme(panel.grid.major = element_blank() )
d <- d + annotate(geom = "text",x = 150, y = 0.2,
                 label = "paste(italic(R) ^ 2, \" = 0.47, p < 2e-16 \")",
                 parse = TRUE,
                 color = "black", size = 3)
d  <- d + theme(legend.position = "none")

d

ggsave("/path/to/folder/bst_vs_dist_from_southmost_lm_BC.png", width = 6, height = 5, units = "in")
