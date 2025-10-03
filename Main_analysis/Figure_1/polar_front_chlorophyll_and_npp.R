setwd("/env/cns/scratch_TaraOcean/BioAdvection_II/MetaT_4/Groups_metaT/Polar_front/")
#setwd("~/docs_from_inti/Polar_front")
library('ncdf4')
#library("rgdal")
library("raster")
library("ggplot2")
library('readxl')
library('gridExtra')
library('tidyr')
#library('DT')
library('dplyr')
library('RColorBrewer')
library('matlab')
library('maptools')
#library('rgeos')
library('mapproj')
library('ggrepel')
#library('SDMTools')
library('viridis')
library('pals')

type=commandArgs(trailingOnly = T)[1]
data=readRDS(paste('data_',type,'_stanford.rds', sep=''))
data=data[data$la>50,]

proj <- "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"

if (type=='chl') {
  col_vec=parula(99)
} else if (type=='npp'){
  col_vec=viridis(99)
}


data("wrld_simpl", package = "maptools")
#rgeos::set_RGEOS_CheckValidity(2L)
wm_ggplot <- crop(wrld_simpl, extent(-180, 180, 50, 90))
# Defines the x axes required
x_lines <- seq(-120,180, by = 60)

arctic_map <- function(longs, lats, cols, sizes,  colo, lbs, lims, brks,na, longs1, lats1, cols1,sizes1, labels1,bw, datal, datal_bis){
  v <- ggplot() +
    #aes(x = longs, y = lats, z = cols) +
    #geom_tile(aes(fill=cols)) +
    #stat_contour(aes(x = longs, y = lats, z = cols,fill=..level..), geom="polygon", binwidth=1) +

    #stat_contour(aes(x = longs, y = lats, z = cols), binwidth=1, colour='black') +
    geom_tile(aes(x = longs, y = lats,fill=cols))+
    scale_fill_gradientn(colours=colo,labels=lbs, limits=lims, breaks=brks,name=na, guide="none")+
    geom_contour(color="grey", alpha=0.5) +
    #stat_contour(aes(x = longs, y = lats, z = cols), binwidth=bw, colour='black') +
    # geom_contour(geom="polygon", aes(x=longs, y=lats, z=cols, fill= ..level..))+

    #geom_point(aes(x=longs, y=lats), shape=15, fill=cols, color=cols,
    #          size=sizes) +
    # geom_line(data = data,aes(x=lo, y=la, group=grp, color=dists))+
    #scale_colour_gradientn(colors=colo,
    #                       labels=lbs, limits=lims, breaks=brks,
    #                       name=na) +
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
      geom_segment(aes(y = 50, yend = 90, x = x_lines, xend = x_lines), linetype = "dashed", size=0.3) +
      # Change theme to remove axes and ticks
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks=element_blank())  +
      geom_point(aes(x=longs1, y=lats1), shape=16, fill=cols1, color=cols1,
               size=sizes1)
    #geom_label_repel(aes(x=longs1, y=lats1,label=labels1), size=5)
  #                        high = "red",name='Tmin\n(years)',
  #                        labels=c(0,0.5,1,1.5),breaks=c(0,0.5,1,1.5),
  #                        limits=c(0,1.5)) +

  # geom_point(aes(x=, y=), shape=19, fill=cols[1], color=cols[1],
  #            sizes=c(min(sizes),min(sizes)+(max(sizes)-min(sizes))/2, max(sizes)))
  return(v)
}

val=log10(data$val)
mino=10^min(val, na.rm=T)
maxo=10^max(val, na.rm=T)
mid=10^( (max(val, na.rm = T)-min(val, na.rm = T))/2 +min(val, na.rm = T) )
quart= 10^( 3*(max(val, na.rm = T)-min(val, na.rm = T))/4 +min(val, na.rm = T) )
quart1=10^( (max(val, na.rm = T)-min(val, na.rm = T))/4 +min(val, na.rm = T) )

df=data.frame(x=0,y=0, z=1)
pdf(paste('legend_', type, '_stanford.pdf', sep=''))
ggplot(df, aes(x, y, fill = z)) + geom_raster() +
  scale_fill_gradientn(colours=col_vec,na.value = "transparent",
                       breaks=c(0,0.25,0.5,0.75,1),labels=c(round(mino, 3), round(quart1, 2),round(mid, 1),
                                                            round(quart, 1),round(maxo, 0)),
                       limits=c(0,1))
dev.off()


env <- read.table('../env_arctic_3.txt')
env <- env[env$Station %in% c(155:168),]
lats1 <- env$Latitude
lons1 <- env$Longitude
labels1 <- paste(env$Station, '', sep='')

#pdf(paste('polar_fronti_hist_',type,'_stanford.pdf', sep=''))
#hist(val)
#dev.off()

data_lines=data.frame(x=rep(seq(-180, 180, 360/999), 4), y=c(rep(55, 1000), rep(65, 1000), rep(75, 1000), rep(85, 1000)),
                      gp=c(rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000)))

data_lines_bis=data.frame(x=seq(-180, 180, 360/999), y=50, gp=1)

mini=1
maxi=100
lo0=data$lo
lt0=data$la
dat0=data$colo
bw=1

vals_plot=10^val
plotted_dat=data.frame('lat'=lt0, 'lon'=lo0, 'val'=vals_plot)
saveRDS(plotted_dat, paste('plotted_data_', type, '_stanford.rds', sep=''))

pdf(paste('polar_front_',type,'_stanford.pdf', sep=''))
a <-arctic_map(longs = lo0, lats = lt0, cols = dat0,
               sizes = rep(0.01, length(lt0)),
               colo = col_vec,
               lbs = round(c(mini,  maxi),1), lims = c(mini, maxi),
               brks = c(mini, maxi),
               na = type,
               longs1 = lons1, lats1 =lats1, cols1 =rep('black' , length(lats1)),
               sizes1 = rep(2, length(lats1)), labels1=labels1, bw=bw, data_lines, data_lines_bis)
print(a)
dev.off()


