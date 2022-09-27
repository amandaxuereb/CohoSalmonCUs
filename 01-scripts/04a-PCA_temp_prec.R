## This script is used to perform a PCA on temperature and precipitation variables for use in the RDA 
## Amanda Xuereb, Sept. 2022
## modified from scripts by Quentin Rougemont 

rm(list = ls())
## load libraries
libs <- c('dplyr','reshape','ade4','data.table', 'magrittr', 'factoextra','vegan', 'cowplot', 'corrplot')
invisible(lapply(libs, library, character.only = TRUE))

# load the population map 
pop.map <- read.table("/path/to/file/01-pop_map_BC_CU.txt", header = F, stringsAsFactors = F)
head(pop.map)
pop.map <- pop.map[,-3]
colnames(pop.map) <- c("SITE", "INDIVIDUAL")
pop <- as.data.frame(unique(pop.map$SITE), stringsAsFactors = F)
colnames(pop) <- "SITE"

# load the climate data 
bioclim <- read.table("/path/to/folder/climat_epic4_wanted_pop.txt", header = T, stringsAsFactors = F)
head(bioclim)

bioclim <- merge(pop, bioclim, by = "SITE", sort = F)
head(bioclim)
nrow(bioclim)
bioclim$SITE 
bioclim$SITE == pop$SITE

## PERFORM PCA ON TEMPERATURE AND PRECIPITATION VARIABLES

X.temp <- dudi.pca(df = bioclim[, 3:57], center = T, scale = T, scannf = F)
X.prec <- dudi.pca(df = bioclim[,58:97], center = T, scale = T, scannf = F)

eig.val <- get_eigenvalue(X.temp)
eig.val$eig <- eig.val$variance.percent/100
expected <- bstick(length(eig.val$eig))
signif <- eig.val$eig > expected
signif  # number of significant axes

eig.val.p <- get_eigenvalue(X.prec)
eig.val.p$eig <- eig.val.p$variance.percent/100
expected <- bstick(length(eig.val.p$eig))
signif.p <- eig.val.p$eig > expected
signif.p  

# Run PCA again, this time choosing 3 axes for both precipitaiton and temperature (nf = 3), for BC region
# for Thompson region, choose 3 axes for temperature and 2 axes for precipitation (change nf)

X.temp <- dudi.pca(df = bioclim[, 3:57], center = T, scale = T, scannf = F, nf = 3)
X.prec <- dudi.pca(df = bioclim[,58:97], center = T, scale = T, scannf = F, nf = 3)

fviz_pca_ind(X.prec, col.ind="cos2", geom = "point") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red", midpoint = 0.6) +
  theme_minimal()

p12 <- fviz_pca_var(X.prec, col.var = "steelblue") +
  theme_minimal() +
  ggtitle("Precipitation PC axis 12\nvariables") +
  theme(plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5))

p23 <- fviz_pca_var(X.prec, axes = c(2,3), col.var = "steelblue") +
  theme_minimal() +
  ggtitle("Precipitation PC axis 23\nvariables") +
  theme(plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5))

t12 <- fviz_pca_var(X.temp, col.var = "steelblue") +
  theme_minimal() +
  ggtitle("Temperature PC axis 12\nvariables") +
  theme(plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5))


t23 <- fviz_pca_var(X.temp, axes = c(2,3), col.var = "steelblue") +
  theme_minimal() +
  ggtitle("Temperature PC axis 23\nvariables") +
  theme(plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5))


t_eig <- fviz_eig(X.temp, addlabels = TRUE) +
  ggtitle("Temperature eigenvalues") +
  theme(plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5))

p_eig <- fviz_eig(X.prec, addlabels = TRUE, ylim = c(0,50)) +
  ggtitle("Precipitation eigenvalues") +
  theme(plot.title = element_text(lineheight = .8, face = "bold", hjust = 0.5))


t_eig
p_eig

pdf(file = "/path/to/folder/temperature_and_precipitation_pca_BC.pdf", 18, 12)
plot_grid(t12, t23, p12, p23, t_eig, p_eig, labels = "auto", ncol = 4)  # remove p23 for Thompson
dev.off()

# keep the 3 significant axes
prec.ind <- get_pca_ind(X.prec)
temp.ind <- get_pca_ind(X.temp)

colnames(prec.ind$coord) <- paste0("PrecipitationPC", seq(1,3,1))
colnames(temp.ind$coord) <- paste0("TemperaturePC", seq(1,3,1))

prec <- cbind(bioclim$SITE, prec.ind$coord)
temp <- cbind(bioclim$SITE, temp.ind$coord)

colnames(prec)[1] <- colnames(temp)[1] <- "SITE"
head(prec)
head(temp)

prec$SITE == pop$SITE
temp$SITE == pop$SITE
pop$SITE

# write the temperature and precipitation PC coordinates data to file. These will be the variables to use in the RDA.
wd <- "../path/to/folder/"
write.table(temp, paste0(wd, "temperature_coordinates_pca_BC.txt"), quote = F, row.names = F)
write.table(prec, paste0(wd, "precipitation_coordinates_pca_BC.txt"), quote = F, row.names = F)
