##ArcGIS_Acorn_Ant ##

#bridge R and ArcGISPro
library(arcgisbinding)
install.packages("raster")
library(raster)
arc.check_product()

setwd("C:/Users/prile/Box/Research/AcornAntOverwintering2021/NLCD_NV35fSbJOw5X2dG0XKy7")

library(raster)
install.packages("rgdal")
library(rgdal)
# Follow SED's instructions / code for raster path for ISA
# set to your own file path (a copy is in box; downloaded cleveland + farm tile from mrlc online https://www.mrlc.gov/viewer/)

tmp<- raster("NLCD_2019_Impervious_L48_20210604_NV35fSbJOw5X2dG0XKy7.tiff")
# Urban: CWRU S. Campus lat  / lon
#41.502859, -81.598711

# Urban: Ambler / Shaker / Doan Brook
#41.490437, -81.585728

# Urban: Forest Hills
#41.519641, -81.576117

# Rural: Manor House woods
#41.49649580062004, -81.4208774786056

# Rural: Western Land Conservancy
# 41.452935, -81.411649

lat.lons <- SpatialPoints(
  cbind(c(-81.598711, -81.585728, -81.57777, -81.4208774786056,-81.411649),
        c(41.502859, 41.490437, 41.519711, 41.49649580062004, 41.452935)),
  proj4string = CRS("+proj=longlat +datum=WGS84"))

reproj.coords <- spTransform(lat.lons,crs(tmp))

isa.vals<- extract(tmp, reproj.coords)

isa.vals60<- extract(tmp, reproj.coords, buffer=60)
lapply(isa.vals60, median)
lapply(isa.vals60, mean)

isa.vals360 <- extract(tmp, reproj.coords, buffer = 360)
lapply(isa.vals360, median)
lapply(isa.vals360, mean)
#median 120
#[[1]] CWRU
#[1] 35

#[[2]] Shaker / Ambler
#[1] 11

#[[3]] Forest  Hills
#[1] 6

#[[4]] Farm
#[1] 0

#[[5]] WLC
#[1] 0

#means 120
#[[1]] CWRU S
#[1] 35.98039

#[[2]] Shaker / Ambler
#[1] 13.96078

#[[3]] forest hills
#[1] 7.384615

#[[4]] Farm
#[1] 0

#[[5]] #western Land conservancy
#[1] 2.14

#median 300
#[[1]]
#[1] 49

#[[2]]
#[1] 24

#[[3]]
#[1] 10

#[[4]]
#[1] 0

#[[5]]
#[1] 0

#means 300
#[[1]]
#[1] 47.85489

#[[2]]
#[1] 25.74679

#[[3]]
#[1] 22.00321

#[[4]]
#[1] 0.1012658

#[[5]]
#[1] 5.987302


#median 360
#[[1]]
#[1] 51

#[[2]]
#[1] 26

#[[3]]
#[1] 12

#[[4]]
#[1] 0

#[[5]]
#[1] 0

#means
#[[1]]
#[1] 49.65188

#[[2]]
#[1] 28.30022

#[[3]]
#[1] 24.66593

#[[4]]
#[1] 0.3303965

#[[5]]
#[1] 6.298013






















#load nlcd map data into R from GIS
arc_map <- arc.open("nlcd_2019_impervious_descriptor_l48_20210604.img")
dim(arc_map)

#adjust spatial resolution to 1/4 (to match 30 m change to 120 m)
arc_map_res <- as.raster(arc.raster(arc_map, nrow = 26106, ncol = 40297, na.rm = TRUE))

#plot
raster::plot(arc_map_res)


