
"
/***************************************************************************
 Winter Wheat Agriculture Example - Geomundus 2022 Worhskop tutorial
                                
                                
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

Location = "CKA"
Date = "20210614"

# input las at X date during the growing season
las <- readLAS("D:/Geomundus_2022/Winter_Wheat_Plots/20210705.las", select = "xyzcarni")

# input las of bare soil before growing season 
#las2 <- readLAS("E:/CKA/LiDAR/Aligned_las/20210419_CKA_WinterWheat_2.las", select = "xyzcarni")  

resolution =  .15 #define resolution in meters for rasters
k <- .35       #extinction coefficient

"
----------------------------------------------------------------------------------------------------------
"
#Segment ground  
  
mycsf <- csf(sloop_smooth = FALSE, class_threshold = .1, cloth_resolution = .1, time_step = .65) #classify ground in the bar soil las
las_clas1 <- classify_ground(las, mycsf)
gnd <- filter_ground(las_clas1)
plot(gnd)
"
---------------------------------
#Crop Height Model
-------------------------------
"
  
dtm_tin <- rasterize_terrain(las_clas1, res = resolution, algorithm = tin()) #creates DTM
dsm <- rasterize_canopy(las, res = resolution, algorithm = p2r()) #creates dsm
lidar_chm <- dsm - dtm_tin #creates CHM


#plot the resulting chm
plot(lidar_chm,
     main = paste("LiDAR Canopy Height Model (CHM)-", Location, Date))

#retrieve average height
avg_height <- cellStats(lidar_chm, 'mean') 
avg_height

-------------------------------------------
#LiDAR Intensity
-------------------------------------------

canopy = filter_poi(las_clas1, Classification != 2L & ReturnNumber == 1L) 

canopy_intensity <- grid_metrics(canopy, ~mean(Intensity), resolution) # calculate Intensity
plot(canopy_intensity, col = gray.colors(50,0,1),  main = paste("Canopy Intensity-", resolution, "Grid Size"))

-------------------------------------------
#  LiDAR GF & PAI Estimation
-------------------------------------------

#parameterize the Grid Resolution (GR) based on the crop height  
if(avg_height < 0.4){gr <- 0.1
} else if (avg_height < 0.6){gr <- 3
} else if (avg_height < 0.8){gr <- 4
} else if (avg_height < 1){gr <- 6}

gr

# ground point segmentation
mycsf <- csf(sloop_smooth = FALSE, class_threshold = .1, cloth_resolution = .1, time_step = .65)
las_clas2 <- classify_ground(las, mycsf) 


#plot segmentation to verify correct classification
plot_crossection <- function(las_clas2,
                             p1 = c(min(las_clas2@data$X), mean(las_clas2@data$Y)),
                             p2 = c(max(las_clas2@data$X), mean(las_clas2@data$Y)),
                             width = 2, colour_by = NULL)
{
  colour_by <- enquo(colour_by)
  data_clip <- clip_transect(las_clas2, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = .6) + coord_equal() + theme_minimal()
  
  if (!is.null(colour_by))
    p <- p + aes(color = !!colour_by) + labs(color = "")
  
  return(p)
}

#plot crossection of ground to other point segmentation
plot_crossection(las_clas2, colour_by = factor(Classification))

#filter points to only ground points
gnd <- filter_ground(las_clas2)
#plot(gnd, size = 3, bg = "white")  

#grid ground returns
grnd_returns <- grid_metrics(gnd, ~length(Z), resolution) # calculate density
x <- reclassify(grnd_returns, cbind(-Inf, NA, 1), right=TRUE)
plot(x, col = gray.colors(50,0,1),  main = paste("Ground Return Count-", resolution, "Grid Size")) # some plotting

#grid all returns
all_returns <- grid_metrics(las, ~length(Z), resolution) # calculate density
plot(all_returns, col = gray.colors(50,0,1),  main = paste("All Returns Count-", resolution, "Grid Size")) # some plotting

#calculate Gap Fraction
GF <- x/all_returns
plot(GF , col = gray.colors(50,0,1), main = paste("Gap Fraction -", resolution, "Grid Size")) # some plotting

#process for extracting average scan angle
#Grid the mean scan angle
angl_mean <- grid_metrics(las, mean(abs(ScanAngleRank)), resolution)
angl_mean
plot(angl_mean, col = gray.colors(50,0,1), main = paste("Average Scan Angle", resolution,"m", "Grid Size"))
layer_angl_mean <- cellStats(angl_mean, stat='mean', na.rm=TRUE) #find mean scan angle of all cells 
layer_angl_mean
angl_rad <- layer_angl_mean * pi/180 # convert to radians
angl_rad

PAI_LiDAR_r = -((cos(angl_rad)*log(GF))/k)
plot(PAI_LiDAR_r)

#3D Plots
#plot(nongnd, color = "Intensity", bg="white", axis = TRUE, legend = TRUE)
#plot(las_clas2, color = "Classification", bg="white", axis = TRUE)
#plot(chm_las, color = "Z", bg="white", axis = TRUE, legend = TRUE)
#plot(gnd)
#plot(gnd_chm)

#save results to file location
writeRaster(lidar_chm, "C:/Users/jsbat/Documents/LiDAR_output/LiDAR_CHM", options=c('TFW=YES'))
writeRaster(GF, "C:/Users/jsbat/Documents/LiDAR_output/LiDAR_GF", options=c('TFW=YES'))
writeRaster(Int, "C:/Users/jsbat/Documents/LiDAR_output/LiDAR_Int", options=c('TFW=YES'))