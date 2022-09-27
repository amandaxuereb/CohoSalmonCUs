## script for running sNMF with LEA and plotting ancestry barplots 
## Amanda Xuereb, Sept. 2022

#rm(list = ls())

library(vcfR)
library(LEA)
library(adegenet)
library(adegraphics)
library(ape)

vcf <- read.vcfR("/Volumes/Storage/epic4/00-DEC_2020/00-ALL_CLEANED_FOR_MS/00-vcfs_and_other_data/01-all_coho_inds_filtered_m6_p60_x0_S5_mac15_maxmiss95.recode.vcf")

ind.names <- colnames(vcf@gt)[2:ncol(vcf@gt)]
ind.names

# make geno input for snmf

output = vcf2geno("/path/to/input/file/01-all_coho_inds_filtered_m6_p60_x0_S5_mac15_maxmiss95.recode.vcf", "/path/to/output/folder/01-all_coho_inds_filtered_24542snps.geno")

# run snmf with the .geno file
obj.snmf <- snmf("/path/to/input/file/01-all_coho_inds_filtered_24542snps.geno", K = 1:100, project = "new", repetitions = 10, tolerance = 0.00001, entropy = TRUE, ploidy = 2)

obj.snmf <- load.snmfProject("/path/to/folder/01-all_coho_inds_filtered_24542snps.snmfProject")
obj.snmf

plot(obj.snmf, pch = 16, cex = 0.7)
lines(c(11,11), c(0, 0.2595923), col = "blue", lty = "dotted")

k = 11

ce <- cross.entropy(obj.snmf, K = k)
ce
best <- which.min(ce)
best

qmatrix = Q(obj.snmf, K = k, run = best)

# can save qmatrix as .txt file for future use


# rows in qmatrix are the individuals; reorder based on location for plotting
library(readr)

pop.map <- read_tsv("/path/to/file/01-pop_map_allSites_CU.txt")
nrow(pop.map)
head(pop.map)

# assign each pop a number for ordering
popnums <- read.table("/path/to/file/popnames_number.txt", header = T, stringsAsFactors = F)
head(popnums)
nrow(popnums)

class(qmatrix)
nrow(qmatrix)
rownames(qmatrix) <- ind.names
head(qmatrix)

head(pop.map)
qmatrix$INDIVIDUALS <- rownames(qmatrix)
qmatrix.full <- merge(qmatrix, pop.map, by = "INDIVIDUALS", sort = F)
head(qmatrix.full)
nrow(qmatrix.full)

unique(qmatrix.full$CU)
qmatrix.full$CU[qmatrix.full$CU == "CO-1"] <- "CO-01"
qmatrix.full$CU[qmatrix.full$CU == "CO-4"] <- "CO-04"
qmatrix.full$CU[qmatrix.full$CU == "CO-5"] <- "CO-05"
qmatrix.full$CU[qmatrix.full$CU == "CO-7"] <- "CO-07"
qmatrix.full$CU[qmatrix.full$CU == "CO-8"] <- "CO-08"
qmatrix.full$CU[qmatrix.full$CU == "CO-9"] <- "CO-09"

head(qmatrix.full)
colnames(qmatrix.full)[4] <- "SITE"
qmatrix.full.num <- merge(qmatrix.full, popnums, by = "SITE", sort = F)
head(qmatrix.full.num)

qmatrix.ord <- qmatrix.full.num[order(qmatrix.full.num$SITE_NUM),]
head(qmatrix.ord)
unique(qmatrix.ord$CU)

qmatrix.sub <- qmatrix.ord[,3:(k+2)]
head(qmatrix.sub)
head(qmatrix.ord)

barplot(t(qmatrix.sub), col = c("orange", "violet", "green", "blue"), border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients")

# we can see the split between coastal BC and Thompson 
# The next step will be to split panels 

colnames(qmatrix.ord)[3:11] <- c("V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08", "V09")

colnames(qmatrix.ord)[3:4] <- c("V01", "V02")
head(qmatrix.ord)

# set where the panels start and end 
snmf.sorted.1 <- filter(qmatrix.ord, SITE_NUM <= 29)
snmf.sorted.2 <- filter(qmatrix.ord, SITE_NUM >29 & SITE_NUM <= 57)
snmf.sorted.3 <- filter(qmatrix.ord, SITE_NUM >57 & SITE_NUM <= 91)
snmf.sorted.4 <- filter(qmatrix.ord, SITE_NUM >91 & SITE_NUM <= 120)
snmf.sorted.5 <- filter(qmatrix.ord, SITE_NUM >120 & SITE_NUM <= 146)

vector_lines_1 <- vector(mode = "integer", length = 29)
for (i in 1:29) {
  data_filtered_1 <- filter(snmf.sorted.1, SITE_NUM == i)
  if(i ==1) {
    vector_lines_1[i] <- nrow(data_filtered_1)
  } else {
    vector_lines_1[i] <- (nrow(data_filtered_1) + (vector_lines_1[i-1]))
  }
}

vector_lines_2 <- vector(mode = "integer", length = 57)
for (i in 30:57) {
  data_filtered_2 <- filter(snmf.sorted.2, SITE_NUM == i)
  if(i ==1) {
    vector_lines_2[i] <- nrow(data_filtered_2)
  } else {
    vector_lines_2[i] <- (nrow(data_filtered_2) + (vector_lines_2[i-1]))
  }
}

vector_lines_3 <- vector(mode = "integer", length = 91)
for (i in 58:91) {
  data_filtered_3 <- filter(snmf.sorted.3, SITE_NUM == i)
  if(i ==1) {
    vector_lines_3[i] <- nrow(data_filtered_3)
  } else {
    vector_lines_3[i] <- (nrow(data_filtered_3) + (vector_lines_3[i-1]))
  }
}

vector_lines_4 <- vector(mode = "integer", length = 120)
for (i in 92:120) {
  data_filtered_4 <- filter(snmf.sorted.4, SITE_NUM == i)
  if(i ==1) {
    vector_lines_4[i] <- nrow(data_filtered_4)
  } else {
    vector_lines_4[i] <- (nrow(data_filtered_4) + (vector_lines_4[i-1]))
  }
}

vector_lines_5 <- vector(mode = "integer", length = 146)
for (i in 121:146) {
  data_filtered_5 <- filter(snmf.sorted.5, SITE_NUM == i)
  if(i ==1) {
    vector_lines_5[i] <- nrow(data_filtered_5)
  } else {
    vector_lines_5[i] <- (nrow(data_filtered_5) + (vector_lines_5[i-1]))
  }
}


head(snmf.sorted.1)

snmf.sorted.1 <- snmf.sorted.1[,-c(14:16)]
snmf.sorted.2 <- snmf.sorted.2[,-c(14:16)]
snmf.sorted.3 <- snmf.sorted.3[,-c(14:16)]
snmf.sorted.4 <- snmf.sorted.4[,-c(14:16)]
snmf.sorted.5 <- snmf.sorted.5[,-c(14:16)]

snmf.sorted.1 <- snmf.sorted.1[,-c(5:7)]
snmf.sorted.2 <- snmf.sorted.2[,-c(5:7)]
snmf.sorted.3 <- snmf.sorted.3[,-c(5:7)]
snmf.sorted.4 <- snmf.sorted.4[,-c(5:7)]
snmf.sorted.5 <- snmf.sorted.5[,-c(5:7)]

snmf.tidy.1 <- reshape2::melt(snmf.sorted.1, id.vars = c("INDIVIDUALS", "SITE_NUM", "SITE"), variable.name = "ANC", value.name = "PERC")
head(snmf.tidy.1)
snmf.tidy.2 <- reshape2::melt(snmf.sorted.2, id.vars = c("INDIVIDUALS", "SITE_NUM", "SITE"), variable.name = "ANC", value.name = "PERC")
snmf.tidy.3 <- reshape2::melt(snmf.sorted.3, id.vars = c("INDIVIDUALS", "SITE_NUM", "SITE"), variable.name = "ANC", value.name = "PERC")
snmf.tidy.4 <- reshape2::melt(snmf.sorted.4, id.vars = c("INDIVIDUALS", "SITE_NUM", "SITE"), variable.name = "ANC", value.name = "PERC")
snmf.tidy.5 <- reshape2::melt(snmf.sorted.5, id.vars = c("INDIVIDUALS", "SITE_NUM", "SITE"), variable.name = "ANC", value.name = "PERC")

vector_lines_1_start <- c(0, vector_lines_1)

vector_lines_2 <- vector_lines_2[which(vector_lines_2 > 0)]
vector_lines_3 <- vector_lines_3[which(vector_lines_3 > 0)]
vector_lines_4 <- vector_lines_4[which(vector_lines_4 > 0)]
vector_lines_5 <- vector_lines_5[which(vector_lines_5 > 0)]

vector_lines_2_start <- c(0, vector_lines_2)
vector_lines_3_start <- c(0, vector_lines_3)
vector_lines_4_start <- c(0, vector_lines_4)
vector_lines_5_start <- c(0, vector_lines_5)

sites1 <- unique(snmf.tidy.1$SITE)
sites2 <- unique(snmf.tidy.2$SITE)
sites3 <- unique(snmf.tidy.3$SITE)
sites4 <- unique(snmf.tidy.4$SITE)
sites5 <- unique(snmf.tidy.5$SITE)

library(ggplot2)

x_title <- "Individuals"
y_title <- "Admixture coefficient"

display.brewer.pal(n = 11, name = "Paired")
colors <- brewer.pal(n = 11, name = "Paired")
colors

colors <- c("#E6E6FA", "#7FFFD4")
x <- c("V01", "V02")

n <- sprintf('%0.2d', 1:11)
n

x <- paste("V", n, sep = "")
x

anc.colors <- as.data.frame(cbind(colors, x), stringsAsFactors = F)
anc.colors
colnames(anc.colors) <- c("color", "ANC")

head(snmf.tidy.1)
snmf.tidy.1 <- left_join(snmf.tidy.1, anc.colors, by = "ANC")
head(snmf.tidy.1)

snmf.tidy.2 <- left_join(snmf.tidy.2, anc.colors, by = "ANC")
snmf.tidy.3 <- left_join(snmf.tidy.3, anc.colors, by = "ANC")
snmf.tidy.4 <- left_join(snmf.tidy.4, anc.colors, by = "ANC")
snmf.tidy.5 <- left_join(snmf.tidy.5, anc.colors, by = "ANC")


## PANEL 1
p = 1
graph1 <- ggplot(snmf.tidy.1, aes(x = INDIVIDUALS, y = PERC, fill = ANC)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = unique(snmf.tidy.1$color), name = "K", labels = seq(1,2,1)) +
  scale_x_discrete(limits = unique(snmf.tidy.1$INDIVIDUALS)) +
  labs(y = y_title) +
  labs(x = x_title) +
  scale_y_continuous(expand = c(0.1, 0.1), breaks = seq(0, 1, 0.2), labels = seq(0,1,0.2)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), size = 28, family = "Helvetica", face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, family = "Helvetica", face = "bold"),
        axis.text.y = element_text(size = 24, family = "Helvetica", face = "bold"),
        legend.key.size = unit(1, "line"), legend.text = element_text(size = 18, face = "bold"), legend.title = element_text(size = 28, face = "bold"))

for (j in vector_lines_1[-29]) {
  graph1 <- graph1 + geom_segment(x = j, y = 0, xend = j, yend = 1, color = "black", size = 0.3)
}

for (s in 1:length(sites1)) {
  graph1 <- graph1 +
    geom_text(data = NULL, x = (round(vector_lines_1_start[s] + vector_lines_1[s])/2), y = -0.1, label = sites1[s], size = 7, fontface = "bold", angle = 45)
}



ggsave(paste0("/path/to/folder/0",p,"-K", k, ".png"), plot = graph1, width = 30, height = 5)


## PANEL 2
p = 2

graph2 <- ggplot(snmf.tidy.2, aes(x = INDIVIDUALS, y = PERC, fill = ANC)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = unique(snmf.tidy.2$color), name = "K", labels = seq(1,2,1)) +
  scale_x_discrete(limits = unique(snmf.tidy.2$INDIVIDUALS)) +
  labs(y = y_title) +
  labs(x = x_title) +
  scale_y_continuous(expand = c(0.1, 0.1), breaks = seq(0, 1, 0.2), labels = seq(0,1,0.2)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), size = 28, family = "Helvetica", face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, family = "Helvetica", face = "bold"),
        axis.text.y = element_text(size = 24, family = "Helvetica", face = "bold"),
        legend.key.size = unit(1, "line"), legend.text = element_text(size = 18, face = "bold"), legend.title = element_text(size = 28, face = "bold"))

for (j in vector_lines_2[-28]) {
  graph2 <- graph2 + geom_segment(x = j, y = 0, xend = j, yend = 1, color = "black", size = 0.3)
}

for (s in 1:length(sites2)) {
  graph2 <- graph2 +
    geom_text(data = NULL, x = (round(vector_lines_2_start[s] + vector_lines_2[s])/2), y = -0.1, label = sites2[s], size = 7, fontface = "bold", angle = 45)
}

ggsave(paste0("/path/to/folder/0", p, "-K", k, ".png"), plot = graph2, width = 30, height = 5)



## PANEL 3

p = 3

graph3 <- ggplot(snmf.tidy.3, aes(x = INDIVIDUALS, y = PERC, fill = ANC)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = unique(snmf.tidy.3$color), name = "K", labels = seq(1,2,1)) +
  scale_x_discrete(limits = unique(snmf.tidy.3$INDIVIDUALS)) +
  labs(y = y_title) +
  labs(x = x_title) +
  scale_y_continuous(expand = c(0.1, 0.1), breaks = seq(0, 1, 0.2), labels = seq(0,1,0.2)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), size = 28, family = "Helvetica", face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, family = "Helvetica", face = "bold"),
        axis.text.y = element_text(size = 24, family = "Helvetica", face = "bold"),
        legend.key.size = unit(1, "line"), legend.text = element_text(size = 18, face = "bold"), legend.title = element_text(size = 28, face = "bold"))

for (j in vector_lines_3[-34]) {
  graph3 <- graph3 + geom_segment(x = j, y = 0, xend = j, yend = 1, color = "black", size = 0.3)
}

for (s in 1:length(sites3)) {
  graph3 <- graph3 +
    geom_text(data = NULL, x = (round(vector_lines_3_start[s] + vector_lines_3[s])/2), y = -0.1, label = sites3[s], size = 7, fontface = "bold", angle = 45)
}

ggsave(paste0("/path/to/folder/0",p,"-K", k, "_BCSites_allSnps.png"), plot = graph3, width = 30, height = 5)


# PANEL 4 
p = 4
graph4 <- ggplot(snmf.tidy.4, aes(x = INDIVIDUALS, y = PERC, fill = ANC)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = unique(snmf.tidy.4$color), name = "K", labels = seq(1,2,1)) +
  scale_x_discrete(limits = unique(snmf.tidy.4$INDIVIDUALS)) +
  labs(y = y_title) +
  labs(x = x_title) +
  scale_y_continuous(expand = c(0.1, 0.1), breaks = seq(0, 1, 0.2), labels = seq(0,1,0.2)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), size = 28, family = "Helvetica", face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, family = "Helvetica", face = "bold"),
        axis.text.y = element_text(size = 24, family = "Helvetica", face = "bold"),
        legend.key.size = unit(1, "line"), legend.text = element_text(size = 18, face = "bold"), legend.title = element_text(size = 28, face = "bold"))

for (j in vector_lines_4[-29]) {
  graph4 <- graph4 + geom_segment(x = j, y = 0, xend = j, yend = 1, color = "black", size = 0.3)
}

for (s in 1:length(sites4)) {
  graph4 <- graph4 +
    geom_text(data = NULL, x = (round(vector_lines_4_start[s] + vector_lines_4[s])/2), y = -0.1, label = sites4[s], size = 7, fontface = "bold", angle = 45)
}

ggsave(paste0("/path/to/folder/0",p,"-K", k, "_BCSites_allSnps.png"), plot = graph4, width = 30, height = 5)


# PANEL 5 

p = 5

graph5 <- ggplot(snmf.tidy.5, aes(x = INDIVIDUALS, y = PERC, fill = ANC)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = unique(snmf.tidy.5$color), name = "K", labels = seq(1,2,1)) +
  scale_x_discrete(limits = unique(snmf.tidy.5$INDIVIDUALS)) +
  labs(y = y_title) +
  labs(x = x_title) +
  scale_y_continuous(expand = c(0.1, 0.1), breaks = seq(0, 1, 0.2), labels = seq(0,1,0.2)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), size = 28, family = "Helvetica", face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 28, family = "Helvetica", face = "bold"),
        axis.text.y = element_text(size = 24, family = "Helvetica", face = "bold"),
        legend.key.size = unit(1, "line"), legend.text = element_text(size = 18, face = "bold"), legend.title = element_text(size = 28, face = "bold"))

for (j in vector_lines_5[-26]) {
  graph5 <- graph5 + geom_segment(x = j, y = 0, xend = j, yend = 1, color = "black", size = 0.3)
}

for (s in 1:length(sites5)) {
  graph5 <- graph5 +
    geom_text(data = NULL, x = (round(vector_lines_5_start[s] + vector_lines_5[s])/2), y = -0.1, label = sites5[s], size = 7, fontface = "bold", angle = 45)
}

ggsave(paste0("/path/to/folder/0",p,"-K", k, "_allSites_allSnps.png"), plot = graph5, width = 30, height = 5)

