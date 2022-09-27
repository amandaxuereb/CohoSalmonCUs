### After running the PCA on temperature and precipitation variables and writing the coordinates of the relevant axes, this script is used to perform the RDA 
### Amanda Xuereb, Sept. 2022
### modified from scripts by Quentin Rougemont

rm(list = ls())
## load libraries
libs <- c('dplyr','reshape','ade4','data.table', 'magrittr', 'factoextra','vegan', 'cowplot', 'corrplot', 'ggplot2')
invisible(lapply(libs, library, character.only = TRUE))

# Load in the population map 
pop.map <- read.table("/path/to/folder/01-pop_map_BC_CU.txt", header = F, stringsAsFactors = F)
head(pop.map)
colnames(pop.map) <- c("INDIVIDUAL", "SITE")
pop <- as.data.frame(unique(pop.map$SITE), stringsAsFactors = F)
colnames(pop) <- "SITE"

head(pop)

# Load all environmental data 
pop_lat <- read.table("/path/to/folder/coho_dist_v2.txt", T)
pop_lat <- dplyr::select(pop_lat,POP_ID, elevation, dist_max_km, Latitude, Longitude, Region)

#replace altitude of zero by 1
pop_lat$elevation[pop_lat$elevation == 0.00000 ] <- 1

enviro <- read.table("/path/to/folder/climat_epic4_wanted_pop.txt", T)
length(which(enviro$SITE %in% pop$SITE == TRUE)) #check that names match

geol <- read.table("/path/to/folder/era_rocktype_quanti_v2.txt", T)
length(which(geol$SITE %in% pop$SITE == TRUE)) # check that names match

enviro <- merge(enviro, geol, by="SITE", sort = F)
head(enviro)
nrow(enviro)
enviro$SITE

enviro <- merge(pop, enviro, by = "SITE", sort = F)
nrow(enviro)
enviro$SITE
enviro$SITE == pop$SITE

env1 <- dplyr::select(enviro, SITE, ROCK)
head(env1)
env1$SITE


# Get temp and prec coordinates 
temp <- read.table("/path/to/folder/temperature_coordinates_pca_BC.txt", header = T)
prec <- read.table("/path/to/folder/precipitation_coordinates_pca_BC.txt", header = T)

temp <- filter(temp, SITE %in% env1$SITE)
prec <- filter(prec, SITE %in% env1$SITE)

env1 <- left_join(env1, temp)
env1 <- left_join(env1, prec)

head(pop_lat)
colnames(pop_lat)[1] <- "SITE"

env2 <- merge(env1, pop_lat, by = "SITE", sort = F)
env2$SITE == pop$SITE # check that names match

#now take into account the standardized elevation * distance interaction
#standardization function:
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
env2$normalized_distance =  env2$elevation * env2$dist_max_km
env2$normalized_distance <- range01(env2$normalized_distance)

nrow(env2)
head(env2)
env2$SITE

# remove unwanted variables:
env <- dplyr::select(env2, -elevation, -Region, -dist_max_km, -Longitude)
head(env)
colnames(env)[2] <- "Geology"

# save the environmental dataframe 
write.table(env, file = "/path/to/folder/ENV_df_for_RDA_BC.txt", quote = F, col.names = T, row.names = F, sep = "\t")

### start here with the cleaned up env dataframe and the SNP frequency data to perform the RDA 
rm(list = ls())

# load the env dataframe 
env <- read.table("/path/to/folder/ENV_df_for_RDA_BC.txt", header = T, stringsAsFactors = F)

# Load the SNP frequency data 
# the freq.strat file is obtained from plink
freq <- fread("/path/to/folder/BC_coho_inds_filtered_freq.frq.strat", header = T) 
freq$CLST[1:120]
length(unique(freq$SNP))
length(unique(freq$CLST))

# keep relevant columns
freq2 <- dplyr::select(freq,SNP,CLST,MAF)
head(freq2)
unique(freq2$CLST)

freq3 <- reshape2::dcast(freq2,CLST~SNP)
dim(freq3)
freq3[1:10,1:10]
freq3$CLST

clst.ord <- as.data.frame(freq3$CLST)
head(clst.ord)
colnames(clst.ord) <- "SITE"
nrow(clst.ord)

env <- merge(clst.ord, env, by = "SITE", sort = F)
head(env)
nrow(env)

# check that the order of SITE in the env dataframe is the same as CLST in freq3
env$SITE == freq3$CLST # should be all TRUE !! 

freq4 <- freq3[,-1] #remove the site label column 
dim(freq4)

snp <- colnames(freq4)
head(snp)
length(snp)


## NOW DO THE RDA
# with latitude condition 
rda1 <- rda(freq4 ~ Geology +
              TemperaturePC1 + TemperaturePC2 + TemperaturePC3 + PrecipitationPC1 +
              PrecipitationPC2 + PrecipitationPC3 + normalized_distance +
              Condition(Latitude) ,
            data=env, scale=T)


# remove latitude condition for Thompson
# rda1 <- rda(freq4 ~ Geology + 
#               TemperaturePC1 + TemperaturePC2 + TemperaturePC3 + PrecipitationPC1 + 
#               PrecipitationPC2 + normalized_distance, 
#             data=env, scale=T)

# test significance
signif.axis <- anova.cca(rda1, by="axis", parallel=5 ) 
signif.axis

signif.marg <- anova.cca(rda1, by="margin", parallel=5)
signif.marg

signif.full <- anova.cca(rda1, parallel=5)
signif.full

# save the rda object to file.. if needed later without having to re-run 
save(rda1, file = "/path/to/folder/rda_LatCorr_BC_24109snps.RData")

RsquareAdj(rda1)
summary(eigenvals(rda1, model= "constrained"))
#check the vif:
vif.cca(rda1)

## Outlier detection
env <- select(env, -Latitude) 
head(env)

signif <- 6 # 6 axes were significant
load.rda <- scores(rda1, choices=c(1:signif), display="species") 

# outlier function (from Forester et al.):
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}

thresh = 3 # detect outliers based on a threshold of +/- 3 SD from mean loading on each axis 

for(i in 1:signif){
  nam <- paste("cand",i,sep="")
  assign(nam, outliers(load.rda[,i],thresh) )
}

#get total number of outliers 
tot = ls(pattern="cand") 
total=NULL
tmp=NULL
for(i  in 1:signif){
  tmp <- length(get(tot[i]))
  total <-rbind(total, tmp)
}
ncand <- sum(total[,1])
ncand

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5,times=length(cand5)), names(cand5), unname(cand5))
cand6 <- cbind.data.frame(rep(6,times=length(cand6)), names(cand6), unname(cand6))

colnames(cand1) <- colnames(cand2)<- c("axis","snp","loading")
colnames(cand3) <- colnames(cand4) <- colnames(cand5) <- colnames(cand6) <- colnames(cand1)
cand <- rbind(cand1,cand2,cand3,cand4,cand5,cand6)

cand$snp <- as.character(cand$snp)
foo <- matrix(nrow=(ncand), ncol=ncol(env[,-1]))  
colnames(foo) <- colnames(env[,-1])

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- freq4[,nam]
  foo[i,] <- apply(env[,-1],2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)
length(cand$snp[duplicated(cand$snp)])  # no. of duplicate detections across axes

foo <- cbind(cand$axis, duplicated(cand$snp))
table(foo[foo[,1]==1,2]) 
table(foo[foo[,1]==2,2])  

cand <- cand[!duplicated(cand$snp),]  #remove duplicates
col<-ncol(cand)

#correlation extraction:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,col+1] <- names(which.max(abs(bar[4:col]))) # gives the variable
  cand[i,col+2] <- max(abs(bar[4:col]))              # gives the correlation
}
colnames(cand)[col+1] <- "predictor"
colnames(cand)[col+2] <- "correlation"
d = data.frame(table(cand$predictor))

d

sum(d$Freq) # number of outliers associated with each predictor 

# write outlier info to file (the number of outliers associated with each variable, and all correlations)
s = "BC"
cond = "LatCorr"
snpset = "24109SNPs"

write.table(d, paste0("/path/to/folder/outlier_per_predictor_", cond, "_", s, "_snpset", snpset,".txt"), quote =F , row.names = F, col.names = F)

write.table(cand,
            paste0("/path/to/folder/candidate_outliers_with_var",thresh,"_", cond, "_", s, "_snpset", snpset, ".txt"),
            quote=F, row.names=F, col.names=T)


# PLOT RDA RESULTS

library(vegan)
library(tidyverse)

# load the saved rda object if not already in environment
load("/path/to/folder/rda_LatCorr_BC_24109snps.RData")
env <- read.table("/path/to/folder/ENV_df_for_RDA_BC.txt", header = T, stringsAsFactors = F)

# get the % for each axis:
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[2]*100,2)   
axis1 <- paste("RDA1 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[5]*100,2)   
axis2 <- paste("RDA2 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[8]*100,2)   
axis3 <- paste("RDA3 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[11]*100,2)   
axis4 <- paste("RDA4 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[14]*100,2)   
axis5 <- paste("RDA5 ",tmp,"%",sep="")                                
tmp <- round(summary(eigenvals(rda1, model= "constrained"))[17]*100,2)   
axis6 <- paste("RDA6 ",tmp,"%",sep="")     

## Outliers: 

cand <- read.table("/path/to/folder/candidate_outliers_with_var3_LatCorr_BC_snpset24109SNPs.txt", T, stringsAsFactors = F)

sel <- data.frame(cbind(cand$snp, cand$predictor, cand$axis))
colnames(sel) <- c("SNP", "predictor", "axis")
nrow(sel)

# colour by predictor: 
col.pred <- rownames(rda1$CCA$v)  %>% 
  as.data.frame() %>% 
  set_colnames(., "SNP") # get the SNP names
head(col.pred)
nrow(col.pred)

col.pred <- left_join(col.pred , sel) #, by=c("SNP"="V1"))
head(col.pred)

# set colour for each predictor
col.pred$color <- ifelse(col.pred$predictor=="Geology",'#1f78b4',
                         ifelse(col.pred$predictor=="TemperaturePC1", '#a6cee3' ,
                                ifelse(col.pred$predictor=="TemperaturePC2", "#33a02c",
                                       ifelse(col.pred$predictor=='TemperaturePC3','#ffff33',
                                              ifelse(col.pred$predictor=="PrecipitationPC1",'#fb9a99',
                                                     ifelse(col.pred$predictor=="PrecipitationPC2",'#6a3d9a',
                                                            ifelse(col.pred$predictor=="PrecipitationPC3",'#e31a1c',
                                                                   ifelse(col.pred$predictor=="normalized_distance",'#b2df8a',
                                                                          "#f1eef6") )))))))

col.pred[is.na(col.pred)] <-"#f1eef6" 

empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
head(col.pred)

bg <- c('#1f78b4','#a6cee3', '#33a02c', '#ffff33','#fb9a99',
        '#6a3d9a','#e31a1c', '#b2df8a') 
leg <- colnames(env)[-c(1,9)] #remove SITE and Latitude

nrow(col.pred)
head(col.pred)

col.pred12 <- col.pred
col.pred12$predictor <- as.character(col.pred12$predictor)
col.pred12$axis <- as.numeric(col.pred12$axis)
col.pred12$predictor[!col.pred12$axis %in% c(1:2)] <- "none"
head(col.pred12)

col.pred12$color[col.pred12$predictor == "none"] <- "#f1eef6"
filter(col.pred12, axis == 1)

leg.bg <- cbind.data.frame(leg, bg)
head(leg.bg)


lat <- 'LatCorr'
site.n <- "BC"
snpset <- "24109SNPs"
thresh <- 3
ax <- 12 # axes 1 and 2

png(paste0("/path/to/folder/snp_rda_", ax, "_var", thresh, "_", site.n, "_", lat, "_snpset", snpset,".png"),900,900)
par(mar=c(5,6,4,2)+0.1)
plot(rda1, type="n", scaling=3,
     cex.lab=3, cex.axis=2.5,
     xlab=axis1, ylab=axis2,
     xlim=c(-1.2,1), xaxs = "i", ylim=c(-1,1), choices = c(1,2))
points(rda1, display="species",
       pch=21, cex=2, col="gray32", bg=col.pred12$color, scaling=3, choices=c(1,2))
text(rda1, scaling=3, display="bp", col="#0868ac", cex=1.5, choices=c(1,2))
legend("topright", legend=leg[which(leg %in% col.pred12$predictor)],
       bty="n", col="gray32", pch=21, cex=2.5, pt.bg=as.character(leg.bg$bg[which(leg.bg$leg %in% col.pred12$predictor)]))
dev.off()


#### for axes 3 and 4 

col.pred34 <- col.pred
col.pred34$predictor <- as.character(col.pred34$predictor)
col.pred34$axis <- as.numeric(col.pred34$axis)
col.pred34$predictor[!col.pred34$axis %in% c(3:4)] <- "none"
head(col.pred34)

col.pred34$color[col.pred34$predictor == "none"] <- "#f1eef6"
filter(col.pred34, axis == 1)


ax <- 34

png(paste0("/path/to/folder/snp_rda_", ax, "_var", thresh, "_", site.n, "_", lat, "_snpset", snpset,".png"),900,900)
par(mar=c(5,5,4,2)+0.1)
plot(rda1, type="n", scaling=3,
     cex.lab=3, cex.axis=2.5,
     xlab=axis3, ylab=axis4,
     xlim=c(-1,1), ylim=c(-1,1), choices = c(3,4))
points(rda1, display="species",
       pch=21, cex=2, col="gray32", bg=col.pred34$color, scaling=3, choices=c(3,4))
text(rda1, scaling=3, display="bp", col="#0868ac", cex=1.5, choices=c(3,4))
legend("bottomright", legend=leg[which(leg %in% col.pred34$predictor)],
       bty="n", col="gray32", pch=21, cex=2.5, pt.bg=as.character(leg.bg$bg[which(leg.bg$leg %in% col.pred34$predictor)]))
dev.off()


# for axes 5 and 6
col.pred56 <- col.pred
col.pred56$predictor <- as.character(col.pred56$predictor)
col.pred56$axis <- as.numeric(col.pred56$axis)
col.pred56$predictor[!col.pred56$axis %in% c(5:6)] <- "none"
head(col.pred56)

col.pred56$color[col.pred56$predictor == "none"] <- "#f1eef6"
filter(col.pred56, axis == 1)


ax <- 56

png(paste0("/path/to/folder/snp_rda_", ax, "_var", thresh, "_", site.n, "_", lat, "_snpset", snpset,".png"),900,900)
par(mar=c(5,5,4,2)+0.1)
plot(rda1, type="n", scaling=3,
     cex.lab=3, cex.axis=2.5,
     xlab=axis5, ylab=axis6,
     xlim=c(-1,1), ylim=c(-1,1), choices = c(5,6))
points(rda1, display="species",
       pch=21, cex=2, col="gray32", bg=col.pred56$color, scaling=3, choices=c(5,6))
text(rda1, scaling=3, display="bp", col="#0868ac", cex=1.4, choices=c(5,6))
legend("bottomleft", legend=leg[which(leg %in% col.pred56$predictor)],
       bty="n", col="gray32", pch=21, cex=2.5, pt.bg=as.character(leg.bg$bg[which(leg.bg$leg %in% col.pred56$predictor)]))
dev.off()                                                                                

### RDA plot by SITES 

# load RDA object if not loaded already
load("/path/to/folder/rda_LatCorr_BC_24109snps.RData")

pop.cu.col <- read_tsv("/path/to/folder/pop_cus_colours.txt")
head(pop.cu.col)
nrow(pop.cu.col)
pop.cu.col$POP

pop.cu.col <- filter(pop.cu.col, POP %in% freq3$CLST)
nrow(pop.cu.col)
pop.cu.col$POP == freq3$CLST

site.ord <- as.data.frame(freq3$CLST)
colnames(site.ord) <- "POP"
head(site.ord)
head(pop.cu.col)

pop.cu.col <- merge(site.ord, pop.cu.col, by = "POP", sort = F)
pop.cu.col$POP == freq3$CLST
pop.cu.col$POP == env$SITE
pop.cu.col

plot(rda1, type = "n", scaling = 2)
points(rda1, pch = 21, display = 'sites', bg = pop.cu.col$color)

ax <- 12
site.n = "BC"
lat = "LatCorr"
snpset = "24109snps"

png(paste0("/path/to/folder/", ax, "_bySITE_", site.n, "_", lat, "_snpset", snpset,".png"),900,900)
#par(mar=c(5,6,4,2)+0.1)
plot(rda1, type="n", scaling=2,
     cex.lab=3, cex.axis=2.5,
     xlab=axis1, ylab=axis2, xlim = c(-12,10))
#xlim=c(-1,1), ylim=c(-1,1), choices = c(1,2))
points(rda1, display="sites",
       pch=21, cex=3, col="black", bg=pop.cu.col$color, scaling=2, choices=c(1,2))
text(rda1, scaling=2, display="bp", col="#0868ac", cex=1.6, choices=c(1,2))
dev.off()

