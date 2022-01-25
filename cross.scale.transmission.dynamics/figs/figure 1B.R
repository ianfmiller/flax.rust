# from https://www.tylermw.com/a-step-by-step-guide-to-making-3d-maps-with-satellite-imagery-in-r/
options(rgl.useNULL = FALSE)
library(rgl)
library(rayshader)
library(sp)
library(raster)
library(scales)

elevation1 = raster::raster("~/Documents/GitHub/flax.rust/data/site.metadata/N38W107.hgt")
elevation2 = raster::raster("~/Documents/GitHub/flax.rust/data/site.metadata/N38W108.hgt")

elevation = raster::merge(elevation1,elevation2)

r = raster::raster("~/Documents/GitHub/flax.rust/data/site.metadata/LC08_L1TP_034033_20210907_20210916_01_T1_B1.TIF")
g = raster::raster("~/Documents/GitHub/flax.rust/data/site.metadata/LC08_L1TP_034033_20210907_20210916_01_T1_B3.TIF")
b = raster::raster("~/Documents/GitHub/flax.rust/data/site.metadata/LC08_L1TP_034033_20210907_20210916_01_T1_B2.TIF")

rbg = raster::stack(r, g, b)
rbg_corrected = sqrt(raster::stack(r, g, b))

raster::crs(r)
raster::crs(elevation)
elevation_utm = raster::projectRaster(elevation, crs = crs(r), method = "bilinear")

bottom_left = c(y=-107.059001, x=38.796994)
top_right   = c(y=-106.85030, x=39.000000)

extent_latlong = sp::SpatialPoints(rbind(bottom_left, top_right), proj4string=sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
extent_utm = sp::spTransform(extent_latlong, raster::crs(elevation_utm))

e = raster::extent(extent_utm)

rgb_cropped = raster::crop(rbg_corrected, e)
elevation_cropped = raster::crop(elevation_utm, e)

names(rgb_cropped) = c("r","g","b")

r_cropped = rayshader::raster_to_matrix(rgb_cropped$r)
g_cropped = rayshader::raster_to_matrix(rgb_cropped$g)
b_cropped = rayshader::raster_to_matrix(rgb_cropped$b)

el_matrix = rayshader::raster_to_matrix(elevation_cropped)

rgb_array = array(0,dim=c(nrow(r_cropped),ncol(r_cropped),3))

rgb_array[,,1] = r_cropped/255 #Red layer
rgb_array[,,2] = g_cropped/255 #Blue layer
rgb_array[,,3] = b_cropped/255 #Green layer

rgb_array = aperm(rgb_array, c(2,1,3))

#plot_map(rgb_array)

rgb_contrast = scales::rescale(rgb_array,to=c(.15,1))

#plot_map(rgb_contrast)

plot_3d(rgb_contrast, el_matrix, windowsize = c(1000,510), zscale = 15,
        zoom=.50, phi=32,theta=220,fov=1, background = 'white')
render_scalebar(limits=c(0, 5, 10),label_unit = "km",position = "W", y=50,offset=50,text_x_offset = 50,
                scale_length = 1-10/(dim(el_matrix)[2]*30/1000)) # set scale based on 30m^2 pixel size


render_compass(position="N",scale_distance = 1)

site.cols<-viridis_pal(alpha=1)(20)[c(20,15,6,1)]
points<-sp::SpatialPoints(data.frame(y=c(-106.868884,-106.996155,-107.019382,-107.021842),x=c(38.821476,38.971273,38.979707,38.967786)),proj4string=sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")) %>% sp::spTransform(raster::crs(elevation_utm))
render_label(heightmap=el_matrix,text="CC",lat=points@coords[1,2],long=points@coords[1,1],extent=e,textcolor=site.cols[1],linecolor = site.cols[1],z=200,zscale=15,textsize = 2,linewidth = 4)
render_label(heightmap=el_matrix,text="BT",lat=points@coords[2,2],long=points@coords[2,1],extent=e,textcolor=site.cols[2],linecolor = site.cols[2],z=750,zscale=15,textsize = 2.5,linewidth = 4)
render_label(heightmap=el_matrix,text="GM",lat=points@coords[3,2],long=points@coords[3,1],extent=e,textcolor=site.cols[3],linecolor = site.cols[3],z=2000,zscale=15,textsize = 2.5,linewidth = 4)
render_label(heightmap=el_matrix,text="HM",lat=points@coords[4,2],long=points@coords[4,1],extent=e,textcolor=site.cols[4],linecolor = site.cols[4],z=2500,zscale=15,textsize = 2.5,linewidth = 4)

render_snapshot() #save as ~/Documents/GitHub/flax.rust/cross.scale.transmission.dynamics/figs/figure.1.files/B.jpg
# export at dimensions 1219x842, then crop to minimum bounding dimensions