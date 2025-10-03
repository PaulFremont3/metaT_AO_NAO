setwd("/env/cns/scratch_TaraOcean/BioAdvection_II/MetaT_4/Groups_metaT/Polar_front/")

library('ncdf4')
library("rgdal")                                                                                                      
library("raster")
library("ggplot2")
library('readxl')
library('gridExtra')
library('tidyr')
library('dplyr')
library('RColorBrewer')
library('matlab')
library('maptools')
library('mapproj')
library('ggrepel')

u=ncdf4::nc_open('sst.day.mean.2013.nc')

sst=ncdf4::ncvar_get(u, 'sst')

june1=31+28+31+30+31+1
june30=31+28+31+30+31+30

sst_sub=sst[,,june1:june30]

sst_sub_mean=apply(sst_sub, c(1,2), mean)

sst_sub_mean_arc=sst_sub_mean[,561:720]

lat=ncdf4::ncvar_get(u, 'lat')
lon=ncdf4::ncvar_get(u, 'lon')
lon[lon>180]=lon[lon>180]-360
temp_dat=sst_sub_mean
temp_dat[,lat<50] <- NA
temp_dat_plot <- as.vector(temp_dat)

temp_dat=sst_sub_mean
temp_dat[,lat<50] <- NA
temp_dat_plot <- as.vector(temp_dat)




grid_coords <- expand.grid(lon=lon, lat=lat)
lons <- grid_coords$lon
lats <- grid_coords$lat
colors <- (temp_dat_plot-
             min(temp_dat_plot, na.rm = T))*99/(max(temp_dat_plot, na.rm = T)-
                                                  min(temp_dat_plot, na.rm = T))+1
col_vec <- jet.colors(100)
colors_plot <- col_vec[colors]



proj <- "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"

data("wrld_simpl", package = "maptools")                                                                            
wm_ggplot <- crop(wrld_simpl, extent(-180, 180, 50, 90)) 
# Defines the x axes required
x_lines <- seq(-120,180, by = 60)

arctic_map <- function(longs, lats, cols, sizes,  colo, lbs, lims, brks,na, longs1, lats1, cols1,sizes1, labels1, datal, datal_bis){
  v <- ggplot() +
  
    geom_tile(aes(x = longs, y = lats,fill=cols))+
    scale_fill_gradientn(colours=colo,labels=lbs, limits=lims, breaks=brks,name=na, guide="none")+
    geom_contour(color="black", alpha=0.5) +
    stat_contour(aes(x = longs, y = lats, z = cols), binwidth=1, colour='black') +

      geom_polygon(data = wm_ggplot, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey", alpha = 1) +
      geom_line(data=datal, aes(x=x, y=y, group=gp),size = 0.1, linetype = 'dashed')+
      geom_line(data=datal_bis, aes(x=x, y=y, group=gp),size = 1) +      

      # Convert to polar coordinates
      coord_map("ortho", orientation = c(90, 0, 0)) +
      scale_y_continuous(breaks = seq(50, 90, by = 5), labels = NULL) +
      
      # Removes Axes and labels
      scale_x_continuous(breaks = NULL) +
      xlab("") + 
      ylab("") +
      
      # Adds labels
      geom_text(aes(x = 180, y = seq(55, 85, by = 10), hjust = -0.2, label = paste0(seq(55, 85, by = 10), "°N"))) +
      geom_text(aes(x = x_lines, y = 60, label = c("120°W", "60°W", "0°", "60°E", "120°E", "180°W"))) +
      
      # Adds axes
      #geom_hline(aes(yintercept = 50), size = 1)  +
      geom_segment(aes(y = 50, yend = 90, x = x_lines, xend = x_lines),size=0.3, linetype = "dashed") +
      
      # Change theme to remove axes and ticks
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks=element_blank())  +
    geom_point(aes(x=longs1, y=lats1), shape=16, fill=cols1, color=cols1, 
               size=sizes1) 

  return(v)
}

mi <- min(temp_dat_plot, na.rm = T)
mx <- max(temp_dat_plot, na.rm = T)
sel <- lats>50 


env <- read.table('../env_arctic_3.txt')
env <- env[env$Station %in% c(155:168),]
lats1 <- env$Latitude
lons1 <- env$Longitude
labels1 <- paste(env$Station, '', sep='')
sel1=!is.infinite(temp_dat_plot) & !is.na(temp_dat_plot)
data_lines=data.frame(x=rep(seq(-180, 180, 360/999), 4), y=c(rep(55, 1000), rep(65, 1000), rep(75, 1000), rep(85, 1000)),
                      gp=c(rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000)))

data_lines_bis=data.frame(x=seq(-180, 180, 360/999), y=50, gp=1)
pdf(paste('polar_front_temperature_noaa_june.pdf', sep=''))
a <-arctic_map(longs = lons[sel & sel1], lats = lats[sel & sel1], cols = temp_dat_plot[sel & sel1],
               sizes = rep(0.01, length(lats[sel & sel1])),
               colo = col_vec,
               lbs = round(c(mi,  mx),1), lims = c(mi, mx),
               brks = c(mi, mx),
               na = 'Temperature (°C)',
               longs1 = lons1, lats1 = lats1, cols1 = rep('black', length(lats1)),
               sizes1 = rep(2, length(lats1)), labels1=labels1, data_lines, data_lines_bis)
print(a)
dev.off()

mino=mi
maxo=mx
quart1=(mx-mi)/4 +mi
quart=3*(mx-mi)/4 +mi
mid=(mx-mi)/2 +mi
pdf(paste('legend_temperature_june_noaa.pdf', sep=''))
ggplot(df, aes(x, y, fill = z)) + geom_raster() +
  scale_fill_gradientn(colours=col_vec,na.value = "transparent",
                       breaks=c(0,0.25,0.5,0.75,1),labels=c(round(mino, 3), round(quart1, 2),round(mid, 1),
                                                            round(quart, 1),round(maxo, 0)),
                       limits=c(0,1))
dev.off()


