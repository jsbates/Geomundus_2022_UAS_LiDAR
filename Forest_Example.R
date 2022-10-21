
"
/***************************************************************************
 Tree Metrics Example - Geomundus 2022 Worhskop tutorial
                                
                                
 Reference: LiDR Package https://r-lidar.github.io/lidRbook/
                              -------------------
        begin                : 2022-10-20
        Author              : Jordan Bates
        email                : j.bates@fz-juelich.de
 ***************************************************************************/
 "

library(lidR)
library(shiny)
library(ggplot2)
library(rlas)
library(raster)
library(rgdal)
library(plot3D)
library(plotly)
library(rgeos)
library(EBImage)
library(ggpubr)


las <- readLAS("D:/Wuestebach/LIDAR/Wuestebach_20220310_Geomundus_Workshop_2022.las", select = "xyzcarni")

#view 3D point cloud of forest segment
plot(las)

# establishing the resolution of any rasters created 
resolution =  1 #define resolution in meters for rasters

#Showing how height when rasterized is currently in mean sea level and needs to be normalized (file, mean height, resolution)
msl_height <- grid_metrics(las, ~mean(Z), resolution) 
plot(msl_height, col = height.colors(50))

#using cloth simulation filter to classify the ground w/ threshold = points to include above the ground, resolution = distance to searched around each point
mycsf <- csf(sloop_smooth = FALSE, class_threshold = .1, cloth_resolution = .5, time_step = .65) #classify ground in the bar soil las
las_clas1 <- classify_ground(las, mycsf)

# check classificaiton in 3D point cloud
plot(las_clas1, color = "Classification", bg="white")

#we can also look at the side profile of the segmentation
plot_crossection <- function(las_clas1,
                             p1 = c(310890, 5598000),
                             p2 = c(310800,5598090),
                             width = 4, colour_by = NULL)
{
  colour_by <- enquo(colour_by)
  data_clip <- clip_transect(las_clas1, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = 1.6) + coord_equal() + theme_minimal()
  
  if (!is.null(colour_by))
    p <- p + aes(color = !!colour_by) + labs(color = "")
  
  return(p)
}

plot_crossection(las_clas1, colour_by = factor(Classification))


#can then create a Digital Terrain Model raster of the below canopy surface w/ interpolating missing cells with using a tin
dtm_tin <- grid_terrain(las_clas1, res = resolution, algorithm = tin()) #creates DTM
plot(dtm_tin)

#can normalize the height in the point cloud by subtracting the DTM
height_las <- las - dtm_tin
plot(height_las, color = "Z", bg="white", axis = TRUE, legend = TRUE)

#Tree height Raster
height <- grid_canopy(height_las, resolution, pitfree(subcircle = 0.01))
plot(height, main = paste("Tree Height (m) -", resolution, "m GSD")) 

#we can also do a profile of the tree height
plot_crossection <- function(height_las,
                             p1 = c(310890, 5598000),
                             p2 = c(310800,5598090),
                             width = 4, colour_by = NULL)
{
  colour_by <- enquo(colour_by)
  data_clip <- clip_transect(height_las, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z, color = Z)) + geom_point(size = 2) + coord_equal() + theme_minimal() + scale_color_gradientn(colours = height.colors(50))
  
  return(p)
}

plot_crossection(height_las)


#save the raster 
writeRaster(height, filename = "D:/Geomundus_2022/Tree_Metrics_Results/Tree_Height", format = "GTiff", overwrite=TRUE)

#tree detection using a variable window size based on tree height

"
f <- function(x) {x * 0.08 + 2}
heights <- seq(5,30,5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,6))
"

#ttops <- find_trees(height_las, lmf(f, hmin=15))
ttops <- locate_trees(las, lmf(ws = 3.5, hmin = 615))

# Plot Heihgt raster with tree detection markers
plot(height, col = height.colors(50))
plot(ttops, add = TRUE)

# View 3D point cloud with 3D sphere markers
x <- plot(height_las, bg = "white", size = 3)
add_treetops3d(x, ttops)
plot(height)

# tree segmentation algorithm
algo <- silva2016(height, ttops)
seg <- segment_trees(las, algo) # segment point cloud
plot(seg, bg = "white", size = 4, color = "treeID") # visualize trees

#Create convex hull around tree crown
crowns <- crown_metrics(seg, func = .stdtreemetrics, geom = "convex")
crowns
plot(height, col = height.colors(50))
plot(crowns["convhull_area"], main = "Crown area (convex hull)")

#Process crown length as new table
crowns_new <- crowns
print(crowns_new)
cat('\n\n')
crowns_new$crown_length <- (sqrt(crowns_new$convhull_area/3.14159))*2
print(crowns_new)
plot(crowns_new["crown_length"], main = "Average crown width (convex hull)")

# can try concave for more precise area and length
metric = tree_metrics(seg, .stdtreemetrics)
hulls  = delineate_crowns(seg, type = c("concave"), concavity = 1)
hulls@data = dplyr::left_join(hulls@data, metric@data)
spplot(hulls, "Z")


poly <- shapefile("D:/Wuestebach/crownwidth_new.shp")
poly 

#inport NDVI raster from multispectral data
NDVI = raster("D:/Geomundus_2022/Prep/Forest_NDVI_clip.tif")
plot(NDVI)

# Use the tree crown polygons to extract the average NDVI of the tree tops and plot with histogram
ex_NDVI <- extract(NDVI, poly, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
ex_NDVI
ex2_NDVI = ex_NDVI[,-1]
ex2_NDVI
ex_num <- as.numeric(unlist(ex2_NDVI))
ex_num
hist(ex_num,xlab = "ID",col = "yellow",border = "blue")

## Use the tree crown polygons to extract the average tree height  and plot with histogram
tree_height = raster("D:/Geomundus_2022/Tree_Metrics_Results/Tree_Height.tif")
plot(tree_height)
ex_height <- extract(tree_height, poly, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
ex_height
ex_height2 = ex_height[,-1]
ex_height2
ex_num_height <- as.numeric(unlist(ex_height2))
ex_num_height
hist(ex_num_height,xlab = "ID",col = "yellow",border = "blue")

#Plot correlation between tree height and NDVI
joined_df <- merge(ex_height, ex_NDVI, by.x = "ID", 
                   by.y = "ID", all.x = TRUE, all.y = FALSE)
head(joined_df)
ggplot( joined_df, aes( x=Tree_Height, y=Forest_NDVI_clip ))+
  geom_point()+
  stat_cor(method = "pearson", label.x = -5, label.y = 30)

#normalize tree height to see trend in graph
maxs_int <- 30
mins_int <- 10
joined_df$Height_norm <- scale(joined_df$Tree_Height, center = mins_int, scale = maxs_int - mins_int)
joined_df
ggplot( joined_df, aes( x=Height_norm, y=Forest_NDVI_clip ))+
  geom_point()+  geom_smooth(method=lm) +
  stat_cor(method = "pearson")