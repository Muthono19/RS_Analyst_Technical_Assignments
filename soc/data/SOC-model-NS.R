#-------------------------------------------------------------------------------
# This code was created by Stella Gachoki as part of remote sensing analyst
#role for the Natural state organisation on 22nd April 2025. The final version 
# was commited to Github on XXXX.
#-------------------------------------------------------------------------------
# The scripts reads the soil organic carbon samples, exploratory data analyis,
# reads three raster files (NDVI, LST and evapotranspiration), extracts the raster data
# to the samples, Runs a ranger (random forest) mode, evaluates perfomance and makes
# spatial predictions
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function to install and load essential packages
load_libraries <- function() {
  # Define a vector of required packages
  libs <- c(
    "dplyr",         # Data manipulation
    "ggplot2",       # Data visualization
    "readr",         # Reading CSV files
    "sf",            # Shapefiles loading
    "ggthemes",      # Extra themes for ggplot2
    "RColorBrewer",   # Color palettes for visualizations
    "spatstat.geom",          # for NN
    "spatstat.explore",    # test spatial randomness
    "terra",
    "car" ,           #for vif
    "VSURF" ,      #for feature elimination
    "caret"  ,      #for model training RF
    "ranger" ,      ### to call ranger
    "lattice",
    "rasterVis"
  )
  to_install <- setdiff(libs, rownames(installed.packages())) #uninstalled packages
  if (length(to_install)) install.packages(to_install)
  invisible(lapply(libs, library, character.only = TRUE))
}

# Run the function to load the required packages
load_libraries()

setwd("C:/Users/gstel/OneDrive/Desktop/My documents/Others/Natural State Assignment/Soil-Organic-Carbon-NS")

#-------------------------------------------------------------------------------
#Exploratory data analysis
#-------------------------------------------------------------------------------

#Load the CSV and the shapefile
soc <- read_csv("SOC_samples.csv")
names(soc)
soc_sf <- st_read("SOC_samples_NS.shp")
soc_proj <- st_transform(soc_sf, crs = 32637)
aoi <- st_read("Study_Extent_NS.shp")
aoi_proj <-  st_transform(aoi, crs = 32637)

# Check the distribution of SOC
ggplot(soc, aes(x = MgC_per_ha)) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "white") +
  theme_minimal() +
  labs(title = "", x = "MgC per ha")+
  theme(legend.position = "right",legend.text = element_text(),
        legend.key.height = unit(0.8, 'cm'), legend.key.width = unit(0.5, 'cm'), 
        legend.background = element_rect(fill = "white"), 
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"))

#make a spatial map of SOC sized by SE. High SE=high variablity or low samples
ggplot() + geom_sf(data = aoi, fill = "lightgray") +
  geom_sf(data = soc_sf, aes(color = MgC_per_ha,size = MgC_SE)) +
  geom_sf_text(data = aoi, aes(label = FIRST_IEBC), size = 3, fontface = "bold") +
  scale_color_viridis_c() + theme_minimal() +
  labs(title = "",x = "Longitude", y = "Latitude",color = "MgC_per_Ha") +
  theme(legend.position = "right",legend.text = element_text(),
        legend.key.height = unit(0.8, 'cm'), legend.key.width = unit(0.5, 'cm'), 
        legend.background = element_rect(fill = "white"), 
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"))

#Test for spatial randomness. R < 1 (clustered); R = 1 (random); R > 1 (evenly)
coords_soc <- st_coordinates(soc_proj) # Coordinates extraction
win <- owin(xrange = range(coords_soc[,1]),yrange = range(coords_soc[,2])) # Bounding box
pp_soc <- ppp(x = coords_soc[,1], y = coords_soc[,2], window = win) #Point pattern
nn_dist_soc <- nndist(pp_soc) #nearest neighbor distances

hist(nn_dist_soc, breaks = 20, col = "lightblue",
     main = "Nearest Neighbor Distances",
     xlab = "Distance (m)")

ce_soc <- clarkevans(pp_soc) # Clark and Evans test
print(ce_soc)

ce_test <- clarkevans.test(pp_soc, correction = "Donnelly")
print(ce_test) # There is significant spatial clustering

sum(nn_dist_soc <= 300) # those less than 300 m apart
close_indices <- which(nn_dist_soc <= 300)
close_sites <- soc_proj[close_indices, ]
close_sites$plot_no 

#-------------------------------------------------------------------------------
#merging samples with the predictor variables min_NDVI, max_NDVI, mean_LST, mean_MNDWI
#mean_evapo
# GEE link for code used to generate the rasters: 
#-------------------------------------------------------------------------------
cov_var <- rast("preds_NS.tif")
names(cov_var)

#scale various rasters as they were extracted from GEE
#NDVI scale factor 0.0001
ndvi_layers <- c("mam_ndvi", "jja_ndvi", "ond_ndvi", "jf_ndvi")
cov_var[[ndvi_layers]] <- cov_var[[ndvi_layers]] * 0.0001

#LST scale factor 0.02, then degrees - 273.15
lst_layers <- c("maxlst","minlst","medlst" )
cov_var[[lst_layers]] <- (cov_var[[lst_layers]] * 0.02) - 273.15 #degrees celcius

#evapotranspiration scale factor 0.02
evapo_layers <- c("maxevapo" , "minevapo" , "medevapo")
cov_var[[evapo_layers]] <- cov_var[[evapo_layers]] * 0.1 #kg/m^2

plot(cov_var)

#project predictors to the project samples coordinates
cov_var.prj <- project(cov_var,soc_proj)

### Extract variables to SOC sites and save as CSv#####
rasterValues_variables <- extract(cov_var.prj, soc_proj)
rasterValues_variables.combine <- cbind(soc_proj,rasterValues_variables)
write.table(rasterValues_variables.combine,file='dBase_SOC_NS_raw.csv', 
            append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)

#-------------------------------------------------------------------------------
#modelling Part after ensuring CSV saved is good.
# this case i just deleted geometry column from the csv
#-------------------------------------------------------------------------------
d <- read.csv("dBase_SOC_NS_final.csv")
names(d)

# define the y (response) and x(predictors)
y.d <- d$MgC_per_ha
x.d <- d[, c(6:17)]

## ***************** GLM Model ************************
#check for multicollinerity between predictors drop predictors when VIF >10
glm_col <- lm(MgC_per_ha ~ ., data = cbind(MgC_per_ha = d$MgC_per_ha, x.d))
glm_step <- step(glm_col, direction = "both", trace = FALSE)
summary(glm_step)
vif(glm_step)

par(mfrow = c(2, 2))
plot(glm_step)

#spatial predictions
soc_pred_glm <- predict(cov_var.prj, model = glm_step, type = "response")
plot(soc_pred_glm, main = "GLM_based prediction of soil organic carbon (Mg/ha)")

#*********************RANGER (random forest)**********************************
##### Feature elimination using VSURF  #####
set.seed(123) # for reproducibility
vsurf.soc <- VSURF(x=x.d,y=y.d ,mtry=6,ntree.thres = 500,
                   nfor.thres = 30,ntree.interp = 300,nfor.interp = 30)
vsurf.soc$varselect.pred # variables improving predictions

pred.var <- data.frame(x.d[,c(7, 4 ,1 ,8)])
names(pred.var)

# formula for RF training
fm.soc <- MgC_per_ha ~ mam_ndvi+maxevapo+maxlst+jja_ndvi

ctrl.soc<- trainControl(method = "repeatedcv", number = 5, repeats=3)

# check which model is best/optimal
ranger.soc = caret::train(fm.soc, data = d, method = "ranger", trControl = ctrl.soc,importance="permutation",
                          tuneGrid = expand.grid(.mtry =seq(2, 4, 1), .splitrule = "variance", .min.node.size = seq(1, 5, 2)),
                          num.trees = 500,max.depth=6,min.bucket=1)

ranger.soc.opt = caret::train(fm.soc, data = d, method = "ranger", trControl = ctrl.soc,importance="permutation",
                          tuneGrid = expand.grid(.mtry =2, .splitrule = "variance", .min.node.size = 5),
                          num.trees = 500,max.depth=6,min.bucket=1)

r2_model_soc <- round(max(ranger.soc$results$Rsquared), 2)
df_pred_act <- data.frame(Predicted = predict(ranger.soc.opt),Observed = d$MgC_per_ha)
rsq_pred_act <- round(cor(df_pred_act$Predicted, df_pred_act$Observed)^2, 2)

ggplot(df_pred_act, aes(x = Observed, y = Predicted)) +
  geom_point(color = "#2C7BB6", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  annotate("text", x = quantile(df_pred_act$Observed, 0.05), 
           y = quantile(df_pred_act$Predicted, 0.95), 
           label = paste("R² (CV) =", round(r2_model_soc, 2), "\nPearson R² =", round(rsq_pred_act, 2)),
           size = 5, fontface = "italic", hjust = 0) +
  labs(x = "Observed SOC (MgC/ha)",y = "Predicted SOC (MgC/ha)",
    title = "") + theme_minimal()+
  theme(legend.position = "right",legend.text = element_text(),
        legend.key.height = unit(0.8, 'cm'), legend.key.width = unit(0.5, 'cm'), 
        legend.background = element_rect(fill = "white"), 
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"))

# variable importance and PDPS
imp.soc <- varImp(ranger.soc.opt)
imp.soc2 <- as.data.frame(imp.soc$importance)
imp.soc2$varnames <- rownames(imp.soc2)

#png(file = "VarImp_Aflatoxin.png", width = 6000, height =3000, units = "px", res = 750, type = "cairo")
ggplot(imp.soc2, aes(x=reorder(varnames, Overall), y=Overall)) +  geom_point(color="blue",size=4)+
  ggtitle("Random forest (ranger)")+ xlab("") + ylab("")+ coord_flip()+theme_classic() + 
  theme(plot.title = element_text(size=16, hjust = 0.5, color="black"),
        text = element_text(size = 16, face="bold",color = "black"))
#dev.off()

#partial dependence plots
library(pdp)
soc.top = topPredictors(ranger.soc.opt,n=4)
pd.soc <- NULL
for (i in soc.top) {
  tmp <- partial(ranger.soc.opt, pred.var = i, data = d)
  names(tmp) <- c("x", "y")
  pd.soc <- rbind(pd.soc, cbind(tmp, predictor = i))
}
pd.soc$predictor <- factor(pd.soc$predictor, levels = unique(pd.soc$predictor))

ggplot(pd.soc, aes(x = x, y = y)) +
  geom_line(color = "steelblue", size = 1) +
  facet_wrap(~ predictor, scales = "free_x") +
  labs(
    title = "Partial Dependence of Predictors on Relative SOC",
    x = "Predictor Value",
    y = "Predicted Relative SOC"
  ) +
  theme_minimal(base_size = 14)

#spatial prediction
ret.cov_var <- c("mam_ndvi", "maxevapo", "maxlst", "jja_ndvi")
ranger.ras <- cov_var.prj[[ret.cov_var]]
plot(ranger.ras)
pred.ranger <- predict(object=ranger.ras,model=ranger.soc.opt,na.rm=T)
gc()
plot(pred.ranger)

#Merge the two predictions for level plot
# Convert SpatRaster to RasterStack
raster_glm <- as(raster::raster(soc_pred_glm), "RasterLayer")
raster_rf <- as(raster::raster(pred.ranger), "RasterLayer")

print(class(raster_glm))
print(class(raster_rf))

# Merge the two predictions into a RasterStack
pred_soc <- raster::stack(raster_glm, raster_rf)
names(pred_soc) <- c("GLM-SOC", "RF-SOC")

# Function to convert SpatialPolygons to coordinates
convert_to_coords <- function(sp_poly) {
  polys <- slot(sp_poly, "polygons")
  coords <- matrix(nrow = 0, ncol = 2)
  for (i in seq_along(polys)) {
    for (j in seq_along(slot(polys[[i]], "Polygons"))) {
      crds <- rbind(slot(slot(polys[[i]], "Polygons")[[j]], "coords"), c(NA, NA))
      coords <- rbind(coords, crds)
    }
  }
  coords
}

aoi_sp <- as(aoi_proj, "Spatial")
coords_aoi <- convert_to_coords(aoi_sp)

clorpalet <- colorRampPalette(brewer.pal(9,"RdYlGn"))
breaks <- c(-40,-30,-20,-10, 0, 10,20,30,40)

levelplot(pred_soc, at = breaks, col.regions = clorpalet(9), layout = c(2, 1),
          panel = function(...) {
            panel.levelplot(...)
            panel.polygon(coords_aoi, border = "black", lwd = 1.5) # Add AOI outline
          },
          margin = FALSE)

###---------------------------- THE END--------------------------------------------








