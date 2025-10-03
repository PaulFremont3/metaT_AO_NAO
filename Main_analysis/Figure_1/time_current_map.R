library("raster")
library("ggplot2")
library('gridExtra')
library('tidyr')

library('dplyr')
library('RColorBrewer')
library('ggrepel')



proj <- "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"

data("wrld_simpl", package = "maptools")                                                                            
wm_ggplot <- crop(wrld_simpl, extent(-180, 180, 25, 90)) 
# Defines the x axes required
x_lines <- seq(-120,180, by = 60)

arctic_map <- function(longs, lats, cols, sizes, data, colo, lbs, lims, brks,na, datal, datal_bis){
  v <- ggplot() +
    geom_polygon(data = wm_ggplot, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey", alpha = 1) +
    geom_line(data=datal, aes(x=x, y=y, group=gp),size = 0.1, linetype = 'dashed')+
    geom_line(data=datal_bis, aes(x=x, y=y, group=gp),size = 1) +
    # Convert to polar coordinates
    coord_map("ortho", orientation = c(90, 0, 0)) +
    scale_y_continuous(breaks = seq(25, 90, by = 5), labels = NULL) +
    
    # Removes Axes and labels
    scale_x_continuous(breaks = NULL) +
    xlab("") + 
    ylab("") +
    
    # Adds labels
    geom_text(aes(x = 180, y = seq(35, 85, by = 20), hjust = -0.2, label = paste0(seq(35, 85, by = 20), "°N"))) +
    geom_text(aes(x = x_lines, y = 39, label = c("120°W", "60°W", "0°", "60°E", "120°E", "180°W"))) +
    
    # Adds axes
    #geom_hline(aes(yintercept = 25), size = 1)  +
    geom_segment(aes(y = 25, yend = 90, x = x_lines, xend = x_lines),size=0.3, linetype = "dashed") +
    
    # Change theme to remove axes and ticks
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks=element_blank()) +
    geom_point(aes(x=longs, y=lats), shape=19, fill=cols, color=cols, 
               size=sizes) +
    geom_line(data = data,aes(x=lo, y=la, group=grp, color=dists))+
    scale_colour_gradientn(colors=colorRampPalette(colors = colo)(100),
                           labels=lbs, limits=lims, breaks=brks,
                           name=na))
  return(v)
}

Tmin1 <- readr::read_tsv("../data/minAijji_tarrive_min_surface_1000.csv") %>% 
  rename(Station_1 = ...1) %>% 
  dplyr::select(Station_1, ends_with("_SUR")) %>% 
  filter(grepl("_SUR", Station_1)) %>%     
  gather(-Station_1, key = "Station_2", value = "Tmin") %>% 
  mutate(type = "Tmin")
Tmin1$Station_1 <- as.character(sapply(Tmin1$Station_1, 
                          FUN = function(x){paste(strsplit(x, split = '_SUR')[[1]],
                                                  collapse = '')}))
Tmin1$Station_2 <- as.character(sapply(Tmin1$Station_2, 
                                       FUN = function(x){paste(strsplit(x, split = '_SUR')[[1]],
                                                               collapse = '')}))


env_arctic <- read.table('../data/env_arctic_3.txt')
env_arctic$Station <- paste(env_arctic$Station, '', sep='')
stations <- env_arctic$Station
good_stat <- stations[stations>'142' & stations <'201']
latis <- env_arctic$Latitude
lonis <- env_arctic$Longitude



data_lines=data.frame(x=rep(seq(-180, 180, 360/999), 3), y=c(rep(35, 1000),rep(55, 1000), rep(75, 1000)),
                      gp=c(rep(1, 1000), rep(2, 1000), rep(3, 1000)))

data_lines_bis=data.frame(x=seq(-180, 180, 360/999), y=25, gp=1)

stations <- c("158"  ,"173", "175", "178","180", "188", "189",
                     "193", "194", "196","143", "144", "145", "146", "147", "148", "149", "150" ,
                       "151", "152",'155','163', '168', '191')
stations<- sort(stations)
stations = stations[!(stations %in% c('143','191', '146', '149'))]

Tminar <- Tmin1[Tmin1$Station_1 %in% stations & Tmin1$Station_2 %in% stations ,]
init_stat <- '144'
subt <-  Tminar[Tminar$Station_1==init_stat & Tminar$Station_2!=init_stat,]
init_data <-subt[which.min(subt$Tmin),]
passed <- c(init_stat)
while (length(unique(passed))!= length(stations)){
  n_stat <- init_data$Station_2[dim(init_data)[1]]
  n_subt <-  Tminar[Tminar$Station_1==n_stat & Tminar$Station_2!=n_stat ,]
  n_dat <- n_subt[which.min(n_subt$Tmin),]
  passed <- append(passed, n_dat$Station_1)
  if ( !(n_dat$Station_2 %in% passed) ){
    init_data <- rbind(init_data, n_dat)
  } else{
    for (i in 2:dim(n_subt)[1]){
      n_dat <-  n_subt[order(n_subt$Tmin, decreasing = F)[i],]
      if ( !(n_dat$Station_2 %in% passed) ){
        init_data <- rbind(init_data,n_dat)
        break
      }
    }
  }
  
}
init_data0 <- c('158', '168' )
init_data$la1 <- env_arctic$Latitude[match(init_data$Station_1, env_arctic$Station)]
init_data$lo1 <- env_arctic$Longitude[match(init_data$Station_1, env_arctic$Station)]
init_data$la2 <- env_arctic$Latitude[match(init_data$Station_2, env_arctic$Station)]
init_data$lo2 <- env_arctic$Longitude[match(init_data$Station_2, env_arctic$Station)]

colfunc<-rev(colorRampPalette(c("red","yellow","springgreen","royalblue"))(101))
init_data$col= colfunc[round(100*(init_data$Tmin-min(init_data$Tmin, na.rm=T))/(max(init_data$Tmin,na.rm=T)-min(init_data$Tmin,na.rm=T))+1)]

arctic_map_bis <- function(longs, lats, cols, sizes, data, colo, lbs, lims, brks,na, labels, Ocean, labs, datal, datal_bis){
  v <- ggplot() +
  geom_polygon(data = wm_ggplot, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey", alpha = 1) +
  geom_line(data=datal, aes(x=x, y=y, group=gp),size = 0.1, linetype = 'dashed')+
  geom_line(data=datal_bis, aes(x=x, y=y, group=gp),size = 1) +
  # Convert to polar coordinates
  coord_map("ortho", orientation = c(90, 0, 0)) +
  scale_y_continuous(breaks = seq(25, 90, by = 5), labels = NULL) +
  
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  
  # Adds labels
  geom_text(aes(x = 180, y = seq(35, 85, by = 20), hjust = -0.2, label = paste0(seq(35, 85, by = 20), "°N"))) +
  geom_text(aes(x = x_lines, y = 39, label = c("120°W", "60°W", "0°", "60°E", "120°E", "180°W"))) +
  
  # Adds axes
  #geom_hline(aes(yintercept = 25), size = 1)  +
  geom_segment(aes(y = 25, yend = 90, x = x_lines, xend = x_lines),size=0.3, linetype = "dashed") +
  
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks=element_blank()) +
  geom_point(aes(x=longs, y=lats), shape=19, fill=cols, color=cols, 
             size=sizes) +
  geom_segment(data=data,
    aes(x = lo1, y = la1, xend = lo2, yend =la2),
    arrow = arrow(
      length = unit(0.02, "npc"), 
      type="open"# Describes arrow head (open or closed)
      
    ),
    colour=data$col,
    size = 1.5
  )  +
  geom_label_repel(aes(x=longs, y=lats,label=labels, colour=Ocean), size=3) + 
  scale_colour_manual(values = labs )
  return(v)
}
Ocean <- rep('', length(stations))
Ocean[stations<='155'] <- 'Atlantic Ocean'
Ocean[stations>'155'] <- 'Arctic Ocean'
labs <- c('Atlantic Ocean'='darkorange', 'Arctic Ocean'='darkblue')
statos = env_arctic$Station
lonis1=lonis[!(statos %in% c('143','191', '146', '149'))]
latis1=latis[!(statos %in% c('143','191', '146', '149'))]
statos = statos[!(statos %in% c('143','191', '146', '149'))]
pdf('time_current_map.pdf')
arctic_map_bis(longs = lonis1, lats = latis1, cols = rep('black', length(latis1)),
               sizes = rep(1, length(latis1)),data = init_data ,
               colo = colfunc,
               lbs = c(min(init_data$Tmin, na.rm=T), max(init_data$Tmin, na.rm=T)), lims = c(min(init_data$Tmin, na.rm=T), max(init_data$Tmin, na.rm=T)), 
               brks = c(min(init_data$Tmin, na.rm=T), max(init_data$Tmin, na.rm=T)),
               na = 'Tmin\n(years)', labels = statos, Ocean=Ocean, labs=labs, datal=data_lines, datal_bis=data_lines_bis)
mi <- round(min(init_data$Tmin, na.rm=T)*12, 2)
ma <- round(max(init_data$Tmin, na.rm=T)*12, 2)
pnt <- cbind(x =c(0,2,2,0), y =c(0,50,50,0))
plot(0,0, col='white', xlim=c(0,17), ylim=c(0, 60), axes=FALSE, frame.plot=F, xlab = '', ylab='')
SDMTools::legend.gradient(pnt, scales::alpha(colfunc, 1), limits=c(mi,ma), title ='Tmin (months)' , cex=1)
dev.off()



u=ncdf4::nc_open('../data/sst.day.mean.2013.nc')

sst=ncdf4::ncvar_get(u, 'sst')

march1=31+28+1
april31=31+28+31+30

sst_sub=sst[,,march1:april31]

sst_sub_mean=apply(sst_sub, c(1,2), mean)

sst_sub_mean_arc=sst_sub_mean[,561:720]

lat=ncdf4::ncvar_get(u, 'lat')
lon=ncdf4::ncvar_get(u, 'lon')
lon[lon>180]=lon[lon>180]-360
temp_dat=sst_sub_mean
temp_dat[,lat<25] <- NA
temp_dat_plot <- as.vector(temp_dat)

temp_dat=sst_sub_mean
temp_dat[,lat<25] <- NA
temp_dat_plot <- as.vector(temp_dat)



grid_coords <- expand.grid(lon=lon, lat=lat)
lons <- grid_coords$lon
lats <- grid_coords$lat
print(min(temp_dat_plot, na.rm = T))
print(max(temp_dat_plot, na.rm = T))

arctic_map_front <- function(longs, lats, cols, sizes, data, colo, lbs, lims, brks,na, labels, Ocean, labs, datal, datal_bis, cols_f, lo, la){
  v <- ggplot() +
  geom_polygon(data = wm_ggplot, aes(x = long, y = lat, group = group), fill = "grey", colour = "grey", alpha = 1) +
  geom_line(data=datal, aes(x=x, y=y, group=gp),size = 0.1, linetype = 'dashed')+
  geom_line(data=datal_bis, aes(x=x, y=y, group=gp),size = 1) +
  stat_contour(aes(x = lo, y = la, z = cols_f), breaks=c(0,3), colour='black') +
  # Convert to polar coordinates
  coord_map("ortho", orientation = c(90, 0, 0)) +
  scale_y_continuous(breaks = seq(25, 90, by = 5), labels = NULL) +

  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  xlab("") +
  ylab("") +

  # Adds labels
  geom_text(aes(x = 180, y = seq(35, 85, by = 20), hjust = -0.2, label = paste0(seq(35, 85, by = 20), "°N"))) +
  geom_text(aes(x = x_lines, y = 39, label = c("120°W", "60°W", "0°", "60°E", "120°E", "180°W"))) +

  # Adds axes
  #geom_hline(aes(yintercept = 25), size = 1)  +
  geom_segment(aes(y = 25, yend = 90, x = x_lines, xend = x_lines),size=0.3, linetype = "dashed") +

  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks=element_blank()) +
  geom_point(aes(x=longs, y=lats), shape=19, fill=cols, color=cols,
             size=sizes) +
  geom_segment(data=data,
    aes(x = lo1, y = la1, xend = lo2, yend =la2),
    arrow = arrow(
      length = unit(0.02, "npc"),
      type="open"# Describes arrow head (open or closed)

    ),
    colour=data$col,
    size = 1.5
  )  +
  geom_label_repel(aes(x=longs, y=lats,label=labels, colour=Ocean), size=3) +
  scale_colour_manual(values = labs )
  return(v)
}

sel <- lats>25
sel1=!is.infinite(temp_dat_plot) & !is.na(temp_dat_plot)
cols_t=temp_dat_plot[sel & sel1]
pdf('time_current_map_and_front.pdf')
arctic_map_front(longs = lonis1, lats = latis1, cols = rep('black', length(latis1)),
               sizes = rep(1, length(latis1)),data = init_data ,
               colo = colfunc,
               lbs = c(min(init_data$Tmin, na.rm=T), max(init_data$Tmin, na.rm=T)), lims = c(min(init_data$Tmin, na.rm=T), max(init_data$Tmin, na.rm=T)),
               brks = c(min(init_data$Tmin, na.rm=T), max(init_data$Tmin, na.rm=T)),
               na = 'Tmin\n(years)', labels = statos, Ocean=Ocean, labs=labs, datal=data_lines, datal_bis=data_lines_bis, cols_f=cols_t, la=lats[sel & sel1],
               lo=lons[sel & sel1])
mi <- round(min(init_data$Tmin, na.rm=T)*12, 2)
ma <- round(max(init_data$Tmin, na.rm=T)*12, 2)
pnt <- cbind(x =c(0,2,2,0), y =c(0,50,50,0))
plot(0,0, col='white', xlim=c(0,17), ylim=c(0, 60), axes=FALSE, frame.plot=F, xlab = '', ylab='')
SDMTools::legend.gradient(pnt, scales::alpha(colfunc, 1), limits=c(mi,ma), title ='Tmin (months)' , cex=1)
dev.off()

