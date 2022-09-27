### Script for performing LFMM 
### Amanda Xuereb, Sept. 2022
### modified from scripts by Quentin Rougemont 

rm(list = ls())

# load libraries
libs <- c('dplyr','reshape2','ade4','data.table', 'magrittr', 'factoextra','vegan', 'lfmm', 'ade4')
invisible(lapply(libs, library, character.only = TRUE))
library(LEA)

vcf2lfmm("/path/to/folder/BC_coho_inds_filtered_m6_p60_x0_S5_miss0.95_mac15.recode.vcf", "/path/to/folder/BC_coho_inds_filtered_m6_p60_x0_S5_miss0.95_mac15.lfmm")

dat <- fread("/path/to/folder/BC_coho_inds_filtered_m6_p60_x0_S5_miss0.95_mac15.lfmm")
dat[1:10,1:10]


pop.map <- read.table("/path/to/folder/01-pop_map_BC_CU.txt", header = F, stringsAsFactors = F)
head(pop.map)
pops <- as.data.frame(unique(pop.map$V1), stringsAsFactors = F)
colnames(pops) <- "SITE"
head(pops)
nrow(pops)


### DOWNLOAD ENV DATA 
pop_lat <- read.table("/path/to/folder/coho_dist_v2.txt",T, stringsAsFactors = F)
pop_lat <- dplyr::select(pop_lat,POP_ID, elevation, dist_max_km, Latitude, Longitude, Region)

#replace altitude of zero by 1
pop_lat$elevation[pop_lat$elevation == 0.00000 ] <- 1

enviro <-read.table("/path/to/folder/climat_epic4_wanted_pop.txt",T, stringsAsFactors = F)
geol   <- read.table("/path/to/folder/era_rocktype_quanti_v2.txt",T)

enviro <- merge(enviro, geol, by="SITE", sort = F)

enviro <- merge(pops, enviro, by = "SITE", sort = F)
enviro$SITE == pops$SITE
enviro$SITE

# Get PC axes for temp and precip
prec <- read.table("/path/to/folder/precipitation_coordinates_pca_BC.txt", header = T, stringsAsFactors = F)
temp <- read.table("/path/to/folder/temperature_coordinates_pca_BC.txt", header = T, stringsAsFactors = F)

env1 <- dplyr::select(enviro, SITE, ROCK)
head(env1)
env1 <- left_join(env1, temp)
env1 <- left_join(env1, prec)
nrow(env1)
head(env1)
env1$SITE

head(pop_lat)
colnames(pop_lat)[1] <- "SITE"
pop_lat$SITE

env2 <- merge(env1, pop_lat, by = "SITE", sort = F)
nrow(env2)
env2$SITE == pops$SITE

#now take into accont the standardized elevation * distance interaction
#standardization function:
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
env2$normalized_distance =  env2$elevation * env2$dist_max_km
env2$normalized_distance <- range01(env2$normalized_distance)

env2 <- left_join(pops, env2)
head(env2)
env2$SITE

# remove unwanted variables:
env <- dplyr::select(env2, -elevation, -Region, -Longitude, -dist_max_km)
head(env)
colnames(env)[2] <- "Geology"
env$SITE


#set them to an individual level:
ind <- pop.map[,c(1,2)]
head(ind)
colnames(ind) <- c("SITE", "IND")
unique(ind$SITE)
nrow(ind)

X <- merge(ind, env, by = "SITE", sort = F) 
head(X)
unique(X$SITE)
nrow(X)
X <-  dplyr::select(X, -SITE, -IND) 

### data imputation for LFMM1:
dat[1:10,1:10]
dat[dat == 9 ] <- NA   
Y <- apply(dat, 
           2, 
           function(x) {
             replace(x, is.na(x), as.numeric(names(which.max(table(x)))))
           })

#PCA:
pc <- prcomp(Y) 
pc

pdf(file=paste0("/path/to/folder/BC_24109snps_pca_dev_withlines.pdf")) 
plot(pc$sdev[1:40]^2, xlab = 'PC', ylab = "Variance explained")
points(6,pc$sdev[2]^2, type = "h", lwd = 1, col = "blue")             
dev.off() 


k = 10
loci <- read.table("/path/to/folder/BC_coho_inds_filtered_m6_p60_x0_S5_miss0.95_mac15.vcfsnp")
dim(loci)
loci[1:10,]

mod.lfmm <- lfmm_ridge(Y = Y,
                       X = X,
                       K = k)

pv <- lfmm_test(Y = Y,
                X = X,
                lfmm = mod.lfmm,
                calibrate = "gif")
pvalues <- pv$calibrated.pvalue

#look at qqplot

qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)


#associate pvalues to loci
pval_loc <- cbind(loci[,c(1:3)],pvalues)                              
pvadj <- matrix(ncol = ncol(pvalues),nrow = nrow(pvalues), NA) 

#correct value with BH method
for(i in 1:ncol(pvalues)){
  pvadj[,i] <- p.adjust(pvalues[,i], method = 'BH')  
}

colnames(pvadj) <- colnames(pvalues)
padj_loc <- cbind(loci[,c(1:3)],pvadj)

head(padj_loc)

#€xport value
write.table(pval_loc,
            paste0("/path/to/folder/pvalues_lfmm_K",k,"_BC_24109snps.txt"),
            sep = "\t", quote = F, row.names = F)

#corrected value
write.table(padj_loc,
            paste0("/path/to/folder/adjust_pvaluesBH_lfmm_K",k,"_BC_24109snps.txt"),
            sep = "\t", quote = F, row.names = F)

# Extract significiant variables using p-value cutoff
Geology <- filter(padj_loc, Geology < 0.01) %>% select(V1,V2,V3,Geology) %>% 
  set_colnames(.,c("CHR","POS","SNP","BH"))
TemperaturePC1 <- filter(padj_loc,TemperaturePC1 < 0.01) %>% select(V1,V2,V3,TemperaturePC1) %>% 
  set_colnames(.,c("CHR","POS","SNP","BH")) 
TemperaturePC2 <- filter(padj_loc,TemperaturePC2 < 0.01) %>% select(V1,V2,V3,TemperaturePC2) %>% 
  set_colnames(.,c("CHR","POS","SNP","BH"))
TemperaturePC3 <- filter(padj_loc,TemperaturePC3 < 0.01) %>% select(V1,V2,V3,TemperaturePC3) %>% 
  set_colnames(.,c("CHR","POS","SNP","BH"))
PrecipitationPC1 <- filter(padj_loc,PrecipitationPC1 < 0.01) %>% select(V1,V2,V3,PrecipitationPC1) %>% 
  set_colnames(.,c("CHR","POS","SNP","BH"))
PrecipitationPC2 <- filter(padj_loc,PrecipitationPC2 < 0.01) %>% select(V1,V2,V3,PrecipitationPC2) %>% 
  set_colnames(.,c("CHR","POS","SNP","BH"))
PrecipitationPC3 <- filter(padj_loc,PrecipitationPC3 < 0.01) %>% select(V1,V2,V3,PrecipitationPC3) %>% 
  set_colnames(.,c("CHR","POS","SNP","BH"))
Latitude <- filter(padj_loc,Latitude < 0.01) %>% select(V1,V2,V3,Latitude) %>% 
  set_colnames(.,c("CHR","POS","SNP","BH"))
Dist <- filter(padj_loc,normalized_distance < 0.01) %>% 
  select(V1,V2,V3,normalized_distance) %>% 
  set_colnames(.,c("CHR","POS","SNP","BH")) 

Geology$var = "Geology"
TemperaturePC1$var = "Temp1"
TemperaturePC2$var = "Temp2"
TemperaturePC3$var = "Temp3"
PrecipitationPC1$var = "Prec1"
PrecipitationPC2$var = "Prec2"
PrecipitationPC3$var = "Prec3"
Latitude$var = "Latitude"
Dist$var = "Dist"

all <- rbind(Geology, TemperaturePC1, TemperaturePC2, TemperaturePC3, 
             PrecipitationPC1, PrecipitationPC2, PrecipitationPC3, Dist)
all

#remove variable that covary with latitude:
all_corrected <- anti_join(all, Latitude, c("SNP" = "SNP"))   

write.table(all_corrected,
            paste0("/path/to/folder/significant_outlier_control_forLatitudeK", k, "_BC_24109snps.txt"),
            quote = F, row.names = F, sep = "\t")

write.table(all,
            paste0("/path/to/folder/significant_outlierK",k,"_BC_24109snps.txt")
            ,quote = F,row.names = F, sep = "\t")

head(all_corrected)
table(all_corrected$var)
nrow(all_corrected)

length(which(duplicated(all_corrected$SNP)==TRUE)) # check duplicates

