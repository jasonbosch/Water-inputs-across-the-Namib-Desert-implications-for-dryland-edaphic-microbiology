#Namib iButton Data Analysis
#This script analyses iButton data collected from the Namib Desert as well as data from SASSCAL and NASA 
#to make inferences about soil moisture and microbial activity.
#Written by: Jason Bosch
#Last update: 04/11/2021

###############
####SET UP#####
###############

##The locations for the results folder and the raw data that needs to be imported will have to be adjusted to suit your setup.
##The data files can be downloaded from Zenodo: https://doi.org/10.5281/zenodo.5577945

#Load the required libraries
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(lubridate)
library(stringr)
library(vegan)
library(ggrepel)
library(ggpubfigs)
library(sf)
library(rgdal)
library(broom)
library(RColorBrewer)
library(ggspatial)
library(rnaturalearth)

#Set working directory
setwd("~/PostDoc/02_Namib_10yr_project/iButton_data/analysis/")

###Import the cleaned data###

#Bring in the iButton data
iButton_complete <- read.csv("iButton_Complete_dataset.csv",stringsAsFactors = F,colClasses = c(Date="POSIXct",Year="character"))

#Makes sure everything stays in the right order instead of alphabetical
iButton_complete$Group <- factor(iButton_complete$Group, levels = unique(iButton_complete$Group))
iButton_complete$Zone <- factor(iButton_complete$Zone, levels = unique(iButton_complete$Zone))

#Need the hour as a separate column where all the dd/MM/YYYY info is gone, in this case all converted to the day of the command
#Has to be done on the day that the script is run or the limits of some graphs will break.
Times <- strftime(iButton_complete$Date, format="%H:%M:%S")
Times <- as.POSIXct(Times, format="%H:%M:%S")
iButton_complete[,"Time"] <- Times

#Bring in NASA data
NASA_data <- read.csv("NASA_iButton_merged_data.csv",colClasses = c(YYYYMMDD="POSIXct"))

#Factor it to keep the right order
NASA_data$Zone <- factor(NASA_data$Zone, levels = unique(NASA_data$Zone), ordered = TRUE)
NASA_data$MM <- factor(NASA_data$MM, levels = sort(unique(NASA_data$MM)), ordered = TRUE)
NASA_data$Site <- factor(NASA_data$Site, levels = sort(unique(NASA_data$Site)), ordered = TRUE)

#Set the seasons (see PCA below for determination)
iButton_complete[month(iButton_complete[,"Date"])%in%c(1,2,3,4,10,11,12),"Season"] <- "Summer"
iButton_complete[month(iButton_complete[,"Date"])%in%c(5,6,7,8,9),"Season"] <- "Winter"
NASA_data[NASA_data[,"MM"]%in%c(1,2,3,4,10,11,12),"Season"] <- "Summer"
NASA_data[NASA_data[,"MM"]%in%c(5,6,7,8,9),"Season"] <- "Winter"

#Bring in weather station data

#Bring in data we need
weather_station_complete <- read.csv("weather_station_complete.csv",stringsAsFactors = F,colClasses = c(Date="POSIXct"))
weather_station_site <- read.csv("weather_station_by_site.csv",stringsAsFactors = F,colClasses = c(Date="POSIXct"))

#Set the season (Determined via PCA plot of relevant values)
weather_station_complete[month(weather_station_complete[,"Date"])%in%c(1,2,3,4,10,11,12),"Season"] <- "Summer"
weather_station_complete[month(weather_station_complete[,"Date"])%in%c(5,6,7,8,9),"Season"] <- "Winter"

#Set the zones as factors for correct ordering in graphs
weather_station_complete$Zone <- factor(weather_station_complete$Zone,levels = c("Fog","Rain"), ordered = T)
weather_station_site$Zone <- factor(weather_station_site$Zone,levels = c("Fog","Rain"), ordered = T)

#Combine all the data sets
Namib_Weather_all_data <- merge(NASA_data,weather_station_site, by.x = c("YYYYMMDD","Site"), by.y = c("Date","Site"))
#Clean up one of the zone names to use in analyses
colnames(Namib_Weather_all_data)[25] <- "Zone"

#Bring in the Shallow Fog site

#Import data
iButton_shallow <-read.csv("iButton_shallow.csv",stringsAsFactors = F,colClasses = c(Date="POSIXct"))

#Set times
Times <- strftime(iButton_shallow[,"Date"], format="%H:%M:%S")
Times <- as.POSIXct(Times, format="%H:%M:%S")
iButton_shallow[,"Time"] <- Times

#Seasons
iButton_shallow[month(iButton_shallow[,"Date"])%in%c(1,2,3,4,10,11,12),"Season"] <- "Summer"
iButton_shallow[month(iButton_shallow[,"Date"])%in%c(5,6,7,8,9),"Season"] <- "Winter"

###################
#####ANALYSIS######
###################

#Set seed
set.seed(123456)

#Set constant colours
#colourbrewer2 used as colour selection aid: https://colorbrewer2.org

rain_cols <- rev(colorRampPalette( rev(brewer.pal(9, "Blues")) )(14))
rain_cols <- rain_cols[-1]
names(rain_cols) <- c("0-50","50-100","100-150","150-200","200-250","250-300","300-350","350-400","400-450","450-500","500-550","550-600","more than 600")

feature_cols <- friendly_pal("muted_nine")[7:9] #c("#1b9e77","#d95f02","#7570b3")
names(feature_cols) <- c("Sample Site","City","Weather Station")

process_cols <- friendly_pal("muted_nine")[1:4]
names(process_cols) <- c("Microbial community growth","Photosynthesis","Amino acid metabolism","Basal metabolism")

year_cols <- c("#0571b0","#92c5de","#ca0020","#f4a582")
names(year_cols) <- c("RH 2018-2019","RH 2019-2020","Temp 2018-2019","Temp 2019-2020")

zone_cols <- friendly_pal("contrast_three")[1:2]
names(zone_cols) <- c("Fog","Rain")

site_cols <- c("#756bb1","#bcbddc","#efedf5","#feedde","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#8c2d04")
names(site_cols) <- sort(unique(Namib_Weather_all_data$Site))

season_cols <- c("#ef8a62","#67a9cf")
names(season_cols) <- c("Summer","Winter")

Variable_cols <- c("#e66101","#5e3c99","#b2abd2","#fdb863","#b2df8a","#1f78b4","#a6cee3","black","red")
names(Variable_cols) <- c("iButton soil temp.","NASA air temp.","SASSCAL air temp.","SASSCAL soil temp.","iButton soil RH (%)","NASA air RH (%)","SASSCAL air RH (%)","SASSCAL rain (mm)","NASA rain (mm)")

##########

#Study Site Map

#Load shapefile with the rainfall. Data comes from: https://www.uni-koeln.de/sfb389/e/e1/download/atlas_namibia/index_e.htm
shapefile <- readOGR(dsn = "../../Map/Namib_rainfall/", layer = "Average annual rainfall")
spdf_fortified <- tidy(shapefile, region = "CLASS")
spdf_fortified$id <- factor(spdf_fortified$id,levels = c("0-50","50-100","100-150","150-200","200-250","250-300","300-350","350-400","400-450","450-500","500-550","550-600","more than 600"))
colnames(spdf_fortified)[7] <- "Rainfall"

#Map Features
map_objects <- as.data.frame(matrix(data = c("Windhoek",-22.57,17.083611,"City",
                                             "Walvis Bay",-22.956111,14.508056,"City",
                                             "Swakopmond",-22.683333,14.533333,"City",
                                             "C14-02",-23.00107,14.67158,"Sample Site",
                                             "C14-04",-23.01800,14.85978,"Sample Site",
                                             "C14-06",-23.06642,15.03977,"Sample Site",
                                             "C14-08",-23.14323,15.20855,"Sample Site",
                                             "C14-10",-23.24627,15.36033,"Sample Site",
                                             "C14-12",-23.31098,15.53390,"Sample Site",
                                             "C14-14",-23.32545,15.71455,"Sample Site",
                                             "C14-16",-23.32073,15.86228,"Sample Site",
                                             "C14-18",-23.34538,16.01060,"Sample Site",
                                             "C14-20",-23.24475,16.14272,"Sample Site",
                                             "Coastal Met",-23.05631,14.625947,"Weather Station",
                                             "Kleinberg",-22.989279,14.72793,"Weather Station",
                                             "Sophies Hoogte",-23.006815,14.890866,"Weather Station",
                                             "Vogelfederberg",-23.097969,15.029032,"Weather Station",
                                             "Garnet Koppie",-23.115385,15.305036,"Weather Station",
                                             "Rooisand",-23.294528,16.114667,"Weather Station",
                                             "Gobabeb",-23.560433333333332,15.041,"Weather Station"),
                                    ncol = 4,dimnames = list(NULL,c("Name","Lat","Lon","Feature")),byrow = TRUE))

NASA_blocks <- as.data.frame(matrix(data = c(14.5,-23.5,1,"NASA 1",
                                             15,-23.5,2,"NASA 1",
                                             15,-23,3,"NASA 1",
                                             14.5,-23,4,"NASA 1",
                                             14.5,-23.5,5,"NASA 1",
                                             15,-23.5,1,"NASA 2",
                                             15.5,-23.5,2,"NASA 2",
                                             15.5,-23,3,"NASA 2",
                                             15,-23,4,"NASA 2",
                                             15,-23.5,5,"NASA 2",
                                             15.5,-23.5,1,"NASA 3",
                                             16,-23.5,2,"NASA 3",
                                             16,-23,3,"NASA 3",
                                             15,-23,4,"NASA 3",
                                             15,-23.5,5,"NASA 3",
                                             16,-23.5,1,"NASA 4",
                                             16.5,-23.5,2,"NASA 4",
                                             16.5,-23,3,"NASA 4",
                                             16,-23,4,"NASA 4",
                                             16,-23.5,5,"NASA 4"),
                                    ncol = 4,dimnames = list(NULL,c("long","lat","piece","group")),byrow = TRUE))
NASA_blocks[,1] <- as.numeric(NASA_blocks[,1])
NASA_blocks[,2] <- as.numeric(NASA_blocks[,2])
NASA_blocks[,3] <- as.numeric(NASA_blocks[,3])

Zoom_box <- as.data.frame(matrix(data = c(14.2,-23.6,1,"Zoom",
                                          16.8,-23.6,2,"Zoom",
                                          16.8,-22.8,3,"Zoom",
                                          14.2,-22.8,4,"Zoom",
                                          14.2,-23.6,5,"Zoom"),
                                 ncol = 4,dimnames = list(NULL,c("long","lat","piece","group")),byrow = TRUE))
Zoom_box[,1] <- as.numeric(Zoom_box[,1])
Zoom_box[,2] <- as.numeric(Zoom_box[,2])
Zoom_box[,3] <- as.numeric(Zoom_box[,3])

#Add major roads
roads <- ne_download(scale = "large",category = "cultural",type = "roads",returnclass = "sf")
roads <- st_crop(roads, xmin = 11, xmax = 25, ymin = -30, ymax = -16)

#Plot maps

#large
large_map <- ggplot(spdf_fortified, aes(x = long, y = lat)) + 
  coord_equal() +
  geom_spatial_polygon(aes(fill = Rainfall, group=group), crs = 4326) +
  geom_sf(data = roads,colour = "light grey",inherit.aes = FALSE) + 
  geom_spatial_polygon(data = Zoom_box,aes(x = long, y = lat),colour = "red",alpha = 0, crs = 4326) +
  geom_spatial_point(data = map_objects, mapping = aes(x=as.numeric(Lon), y= as.numeric(Lat), colour=Feature), crs = 4326) +
  coord_sf(xlim = c(11.50, 26), ylim = c(-30.00, -16.00), expand = FALSE, default_crs = 4326, crs = 4326) +
  scale_fill_manual(values = rain_cols) +
  scale_color_manual(values = feature_cols) +
  geom_spatial_text_repel(data = map_objects[map_objects$Feature%in%c("City"),], mapping = aes(x=as.numeric(Lon), y= as.numeric(Lat), label = Name),max.overlaps = 15,direction = "y", crs = 4326) +
  theme(panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25), panel.background = element_rect(fill = "light grey")) +
  labs( x = "Longitude", y = "Latitude") +
  theme(legend.position = "none") +
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "br", width_hint = 0.5)

#Zoomed
small_map <- ggplot(spdf_fortified, aes(x = long, y = lat)) + 
  coord_equal() +
  geom_spatial_polygon(aes(fill = Rainfall, group=group), crs = 4326) +
  geom_sf(data = roads,colour = "light grey",inherit.aes = FALSE) + 
  geom_spatial_polygon(data = NASA_blocks,aes(x = long, y = lat),colour = "dark gray",alpha = 0, crs = 4326) +
  geom_spatial_point(data = map_objects, mapping = aes(x=as.numeric(Lon), y= as.numeric(Lat), colour=Feature), crs = 4326) +
  geom_spatial_text_repel(data = map_objects[!map_objects$Name%in%c("Swakopmond","Windhoek"),], mapping = aes(x=as.numeric(Lon), y= as.numeric(Lat), label = Name),max.overlaps = 21,direction = "y", crs = 4326, force = 1.5, min.segment.length = 0.5, force_pull = 0.75) +
  coord_sf(xlim = c(14.2, 16.80), ylim = c(-23.6, -22.80), expand = FALSE, default_crs = 4326, crs = 4326) +
  scale_fill_manual(values = rain_cols, name="Annual Rainfall (mm)") +
  scale_color_manual(values = feature_cols) +
  theme(panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25), panel.background = element_rect(fill = "light grey")) +
  labs( x = "Longitude", y = "Latitude") +
  theme(legend.position = "bottom") +
  guides(colour=guide_legend(nrow=3,byrow=TRUE)) +
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "br", width_hint = 0.5)

#Save figure
combined_map <- grid.arrange(large_map + labs(tag = "A"),small_map + labs(tag = "B"), nrow=1,layout_matrix = matrix(c(1,2,2),ncol = 3, byrow=TRUE))
ggsave(filename = paste("Figure 1.svg",sep=""), combined_map, width = 16, height = 5)

##########

#Graph and save the raw data per year for all sites

raw_2019 <- ggplot(iButton_complete[iButton_complete$Year==2019,], aes(x = Date)) + 
  geom_line(aes(y = RH), color = "blue", size = .2) + 
  geom_line(aes(y = Temp), color = "red", size = .2) + 
  scale_x_datetime(breaks = "1 month", minor_breaks = "1 week", date_labels = "%b-%y") + 
  labs(title = "Raw iButton Readings (2018-2019)",  y = NULL) + 
  scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  facet_wrap(~Site,ncol = 2,dir = "v") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1))

raw_2020 <- ggplot(iButton_complete[iButton_complete$Year==2020,], aes(x = Date)) + 
  geom_line(aes(y = RH), color = "blue", size = .2) + 
  geom_line(aes(y = Temp), color = "red", size = .2) + 
  scale_x_datetime(breaks = "1 month", minor_breaks = "1 week", date_labels = "%b-%y") + 
  labs(title = "Raw iButton Readings (2019-2020)",  y = NULL) + 
  scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  facet_wrap(~Site,ncol = 2,dir = "v") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1))

ggsave(filename = "Figure S1.svg", raw_2019, width = 16, height = 9)
ggsave(filename = "Figure S2.svg", raw_2020, width = 16, height = 9)

##########

#Weather Graph/Table
weather_table <- as.data.frame(matrix(nrow = 40,ncol = 11,dimnames = list(NULL,c("Site","Year","Variable","Label","Min","QOne","Median","Mean","QThree","Max","sd"))))
weather_table$Site <- rep(unique(iButton_complete$Site),4)
weather_table$Year <- c(rep(unique(iButton_complete$Year)[1],20),rep(unique(iButton_complete$Year)[2],20))
weather_table$Variable <- rep(c(rep("RH",10),rep("Temp",10)),2)
weather_table$Label <- paste(weather_table$Variable,weather_table$Year,sep = " ")
weather_table$Label <- gsub("2019","2018-2019",weather_table$Label)
weather_table$Label <- gsub("2020","2019-2020",weather_table$Label)

for (n in 1:nrow(weather_table)) {
  weather_table[n,5:11] <- c(summary(iButton_complete[iButton_complete$Site==weather_table$Site[n] & iButton_complete$Year==weather_table$Year[n],weather_table$Variable[n]])[1:6], sd(iButton_complete[iButton_complete$Site==weather_table$Site[n] & iButton_complete$Year==weather_table$Year[n],weather_table$Variable[n]],na.rm = T))
}

#Make cleaner name for plotting and remove the empty Site 14 rows and 
weather_table$Variable <- rep(c(rep("Relative Humidity",10),rep("Temperature",10)),2)
weather_table <- weather_table[c(-27,-37),]

#Graph the table
weather_summary_graph <- ggplot(weather_table,aes(x = Site)) +
  annotate("rect", xmin=2, xmax=6, ymin=0, ymax=100, alpha=0.2, fill="light blue") +
  geom_line(aes(x = Site, y = Min, colour = Label),linetype = "dotted") +
  geom_line(aes(x = Site, y = Max, colour = Label),linetype = "dotted") +
  geom_line(aes(x = Site, y = QOne, colour = Label),linetype = "dashed") +
  geom_line(aes(x = Site, y = QThree, colour = Label),linetype = "dashed") +
  geom_line(aes(x = Site, y = Mean, colour = Label)) +
  scale_y_continuous(limits = c(0,100), n.breaks = 10, 
                     sec.axis = sec_axis(trans = ~ .x,
                                         name = "Temperature (°C)",
                                         breaks = c(0,10,20,30,40,50,60,70,80,90,100))) +
  scale_x_continuous(limits = c(2,20), breaks = c(2,4,6,8,10,12,14,16,18,20)) +
  facet_wrap(~Variable) +
  labs(title = "Annual Weather Conditions", y = "Relative Humidity (%)", colour="") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "bottom") +
  scale_colour_manual(values = year_cols)

ggsave(filename = "Figure 2.svg", weather_summary_graph, width = 8, height = 4.5)

#Save the table
write.table(weather_table,file="weather_table.csv",row.names = F,col.names = T,quote = T,sep = ",")

##########

#Using PCA to find/confirm the data groupings

#The years need to be done separately because of the huge difference in terms of rainfall.
#Due to the fact that Rooisand stopped transmitting humidity data in October 2018, humidity data was excluded from the PCA.

#Prep the dataframes
PCA_df <- subset(Namib_Weather_all_data,select = -Humidity)
PCA_df <- na.omit(PCA_df)
PCA_df_2019 <- PCA_df[PCA_df$YYYYMMDD >= "2018-04-08" & PCA_df$YYYYMMDD <= "2019-03-18",]
PCA_df_2020 <- PCA_df[PCA_df$YYYYMMDD >= "2019-04-30" & PCA_df$YYYYMMDD <= "2020-04-07",]

#Run the PCAs
#2019
PCA_df_2019_sub <- subset(PCA_df_2019, select = c(PRECTOTCORR,QV2M,RH2M,T2M,T2M_MAX,T2M_MIN,iTemp,iRH,Air.temp.avg.,Air.temp.min.,Air.temp.max.,Precipitation,Soil.Temp.avg.))
pca_res_2019 <- prcomp(PCA_df_2019_sub, scale. = TRUE)
dtp_2019 <- data.frame('Site' = as.character(PCA_df_2019$Site),'Zone' = as.character(PCA_df_2019$Zone), 'Month' = as.character(PCA_df_2019$MM),'Season' = as.character(PCA_df_2019$Season), pca_res_2019$x[,1:2]) 
PCAloadings_2019 <- data.frame(Variables = rownames(pca_res_2019$rotation), pca_res_2019$rotation)
#2020
PCA_df_2020_sub <- subset(PCA_df_2020, select = c(PRECTOTCORR,QV2M,RH2M,T2M,T2M_MAX,T2M_MIN,iTemp,iRH,Air.temp.avg.,Air.temp.min.,Air.temp.max.,Precipitation,Soil.Temp.avg.))
pca_res_2020 <- prcomp(PCA_df_2020_sub, scale. = TRUE)
dtp_2020 <- data.frame('Site' = as.character(PCA_df_2020$Site),'Zone' = as.character(PCA_df_2020$Zone), 'Month' = as.character(PCA_df_2020$MM),'Season' = as.character(PCA_df_2020$Season), pca_res_2020$x[,1:2]) 
PCAloadings_2020 <- data.frame(Variables = rownames(pca_res_2020$rotation), pca_res_2020$rotation)

#Keep the ordering
dtp_2019$Zone <- factor(dtp_2019$Zone, levels = unique(dtp_2019$Zone), ordered = TRUE)
dtp_2019$Month <- factor(dtp_2019$Month, levels = sort(as.numeric(unique(dtp_2019$Month))), ordered = TRUE)
dtp_2019$Site <- factor(dtp_2019$Site, levels = sort(as.numeric(unique(dtp_2019$Site))), ordered = TRUE)
dtp_2020$Zone <- factor(dtp_2020$Zone, levels = unique(dtp_2020$Zone), ordered = TRUE)
dtp_2020$Month <- factor(dtp_2020$Month, levels = sort(as.numeric(unique(dtp_2020$Month))), ordered = TRUE)
dtp_2020$Site <- factor(dtp_2020$Site, levels = sort(as.numeric(unique(dtp_2020$Site))), ordered = TRUE)

#Plot the graphs

#Site
#2019
site_graph_2019 <- ggplot(data = dtp_2019) + 
  geom_point(aes(x = PC1, y = PC2, col = Site)) + 
  theme_minimal() + 
  scale_y_continuous(n.breaks = 7) + 
  scale_x_continuous(n.breaks = 6) + 
  scale_color_manual(name = "Site", values = site_cols) +
  labs(title ="Sites (2018-2019)", x = paste("PC1 (",round(summary(pca_res_2019)$importance[2]*100,digits = 2),"%)",sep = ""), y = paste("PC2 (",round(summary(pca_res_2019)$importance[5]*100,digits = 2),"%)",sep = "")) +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
#2020
site_graph_2020 <- ggplot(data = dtp_2020) + 
  geom_point(aes(x = PC1, y = PC2, col = Site)) + 
  theme_minimal() + 
  scale_y_continuous(n.breaks = 6) + 
  scale_x_continuous(n.breaks = 6) + 
  scale_color_manual(name = "Site", values = site_cols) +
  theme(legend.position = "bottom") +
  labs(title ="Sites (2019-2020)", x = paste("PC1 (",round(summary(pca_res_2020)$importance[2]*100,digits = 2),"%)",sep = ""), y = paste("PC2 (",round(summary(pca_res_2020)$importance[5]*100,digits = 2),"%)",sep = "")) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5))

#Month & Season
#2019
month_graph_2019 <- ggplot(data = dtp_2019) + 
  geom_point(aes(x = PC1, y = PC2, col = Month)) + 
  theme_minimal() +
  scale_color_brewer(type = "qual", palette = "Set3") + 
  labs(title ="Months")
#Summer: 1, 2, 3, 4, 10, 11, 12
#Winter: 5, 6, 7, 8, 9
season_graph_2019 <- ggplot(data = dtp_2019) + 
  geom_point(aes(x = PC1, y = PC2, col = Season)) + 
  theme_minimal() + 
  scale_y_continuous(n.breaks = 7) + 
  scale_x_continuous(n.breaks = 6) + 
  scale_colour_manual(values=c("Red","Blue")) + 
  labs(title ="Determined Seasons (2018-2019)", x = paste("PC1 (",round(summary(pca_res_2019)$importance[2]*100,digits = 2),"%)",sep = ""), y = paste("PC2 (",round(summary(pca_res_2019)$importance[5]*100,digits = 2),"%)",sep = "")) +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
#2020
month_graph_2020 <- ggplot(data = dtp_2020) + 
  geom_point(aes(x = PC1, y = PC2, col = Month)) + 
  theme_minimal() +
  scale_color_brewer(type = "qual", palette = "Set3") + 
  labs(title ="Months")
season_graph_2020 <- ggplot(data = dtp_2020) + 
  geom_point(aes(x = PC1, y = PC2, col = Season)) + 
  theme_minimal() + 
  scale_y_continuous(n.breaks = 6) + 
  scale_x_continuous(n.breaks = 6) + 
  scale_colour_manual(values=c("Red","Blue")) + 
  theme(legend.position = "bottom") +
  labs(title ="Determined Seasons (2019-2020)", x = paste("PC1 (",round(summary(pca_res_2020)$importance[2]*100,digits = 2),"%)",sep = ""), y = paste("PC2 (",round(summary(pca_res_2020)$importance[5]*100,digits = 2),"%)",sep = "")) +
  guides(colour=guide_legend(nrow=1,byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5))

#Loading
#2019
loading_graph_2019 <- ggplot(data = PCAloadings_2019) + 
  geom_segment(aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")), color = "black") + 
  geom_text_repel(aes(label=Variables, x = PCAloadings_2019[,2], y=PCAloadings_2019[,3]),max.overlaps = 13) +
  theme_minimal() + 
  scale_y_continuous(n.breaks = 6) + 
  scale_x_continuous(n.breaks = 6) + 
  labs(title ="Loading (2018-2019)", x = "PC1", y = "PC2") +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
#2020
loading_graph_2020 <- ggplot(data = PCAloadings_2020) + 
  geom_segment(aes(x = 0, y = 0, xend = (PC1), yend = (PC2)), arrow = arrow(length = unit(1/2, "picas")), color = "black") + 
  geom_text_repel(aes(label=Variables, x = PCAloadings_2020[,2], y=PCAloadings_2020[,3]),max.overlaps = 13) +
  theme_minimal() +
  scale_y_continuous(n.breaks = 6) + 
  scale_x_continuous(n.breaks = 6) +
  labs(title ="Loading (2019-2020)", x = "PC1", y = "PC2") +
  theme(plot.title = element_text(hjust = 0.5))
  
#Create a combined image
combined_PCA <- grid.arrange(site_graph_2019 + labs(tag = "A"),season_graph_2019 + labs(tag = "B"), loading_graph_2019 + labs(tag = "C"),
                             site_graph_2020 + labs(tag = "D"),season_graph_2020 + labs(tag = "E"), loading_graph_2020 + labs(tag = "F"), nrow=2)
ggsave(filename = paste("Figure 3.svg",sep=""), combined_PCA, width = 16, height = 9)

#Cluster according to euclidean distance
#2019
PCAcoords_2019 <- dtp_2019[,c("PC1","PC2")]
PCAdist_2019 <- dist(PCAcoords_2019)
#2020
PCAcoords_2020 <- dtp_2020[,c("PC1","PC2")]
PCAdist_2020 <- dist(PCAcoords_2020)

#Check for significant differences between zones
#2019
adonis2(PCAdist_2019 ~ Zone, data = dtp_2019, method='eu')
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = PCAdist_2019 ~ Zone, data = dtp_2019, method = "eu")
#            Df SumOfSqs      R2      F Pr(>F)    
# Zone        1   3794.5 0.12121 471.31  0.001 ***
# Residual 3417  27510.1 0.87879                  
# Total    3418  31304.6 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
beta_2019_site <- betadisper(PCAdist_2019, dtp_2019$Zone)
permutest(beta_2019_site)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups       1   28.6 28.5516 15.964    999  0.001 ***
# Residuals 3417 6111.1  1.7884                         
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(beta_2019_site)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#               diff        lwr        upr    p adj
# Fog-Rain -0.199453 -0.2973266 -0.1015795 6.59e-05
#2020
adonis2(PCAdist_2020 ~ Zone, data = dtp_2020, method='eu')
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = PCAdist_2020 ~ Zone, data = dtp_2020, method = "eu")
#             Df SumOfSqs      R2      F Pr(>F)    
# Zone        1     4254 0.14883 528.94  0.001 ***
# Residual 3025    24329 0.85117                  
# Total    3026    28583 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
beta_2020_site <- betadisper(PCAdist_2020, dtp_2020$Zone)
permutest(beta_2020_site)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups       1   20.8 20.8357 11.264    999  0.001 ***
# Residuals 3025 5595.3  1.8497                         
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(beta_2020_site)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#                diff        lwr         upr     p adj
# Rain-Fog -0.1753974 -0.2778659 -0.07292886 0.0007998

#Check for significant differences between seasons
#2019
adonis2(PCAdist_2019 ~ Season, data = dtp_2019, method='eu')
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = PCAdist_2019 ~ Season, data = dtp_2019, method = "eu")
#            Df SumOfSqs      R2      F Pr(>F)    
# Season      1    10695 0.34165 1773.3  0.001 ***
# Residual 3417    20609 0.65835                  
# Total    3418    31305 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
beta_2019_season <- betadisper(PCAdist_2019, dtp_2019$Season)
permutest(beta_2019_season)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups       1   64.7  64.707 41.721    999  0.001 ***
# Residuals 3417 5299.6   1.551                         
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(beta_2019_season)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#                     diff        lwr        upr p adj
# Winter-Summer -0.2766886 -0.3606761 -0.1927011     0
#2020
adonis2(PCAdist_2020 ~ Season, data = dtp_2020, method='eu')
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = PCAdist_2020 ~ Season, data = dtp_2020, method = "eu")
#            Df SumOfSqs      R2      F Pr(>F)    
# Season      1   7347.7 0.25707 1046.7  0.001 ***
# Residual 3025  21235.0 0.74293                  
# Total    3026  28582.6 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
beta_2020_season <- betadisper(PCAdist_2020, dtp_2020$Season)
permutest(beta_2020_season)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups       1    0.0 0.00437 0.0023    999  0.962
# Residuals 3025 5850.3 1.93399             
TukeyHSD(beta_2020_season)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#                       diff        lwr        upr     p adj
# Winter-Summer -0.002412966 -0.1019841 0.09715816 0.9621051

##########

#Compare how well the different data sources give the same information

#Split the data between the two years
comp_data_2019 <- Namib_Weather_all_data[Namib_Weather_all_data$YYYYMMDD >= "2018-04-08" & Namib_Weather_all_data$YYYYMMDD <= "2019-03-18",]
comp_data_2020 <- Namib_Weather_all_data[Namib_Weather_all_data$YYYYMMDD >= "2019-04-30" & Namib_Weather_all_data$YYYYMMDD <= "2020-04-07",]

#Compare iButton RH with NASA RH2M and PRECTOTCORR and weather station
comp_RH <- list()
for (x in c("comp_data_2019","comp_data_2020")) {
  if (str_extract(string = x,pattern = "[0-9]{4}") == "2019") {
    y <- "2018-2019"
  }
  if (str_extract(string = x,pattern = "[0-9]{4}") == "2020") {
    y <- "2019-2020"
  }
  for (n in levels(get(x)$Zone)) {
    graph <- ggplot(get(x)[get(x)$Zone==n,], aes(YYYYMMDD)) + 
      stat_summary(aes(y = iRH, col="iButton soil RH (%)"), fun=mean, geom="line") + 
      stat_summary(aes(y = RH2M, col="NASA air RH (%)"), fun=mean,fun.args = list(na.rm = T), geom="line") + 
      stat_summary(aes(y = Humidity, col="SASSCAL air RH (%)"), fun=mean,fun.args = list(na.rm = T), geom="line") + 
      stat_summary(aes(y = Precipitation * 5, col="SASSCAL rain (mm)"), fun=mean,fun.args = list(na.rm = T), geom="line",alpha = 0.5) + 
      stat_summary(aes(y = PRECTOTCORR * 5, col="NASA rain (mm)"), fun=mean, geom="line",alpha = 0.5) + 
      scale_x_datetime(breaks = "2 month", date_labels = "%b %y") + 
      labs(title = paste("Water availability - ",n," (",y,")",sep = ""), y = "Relative Humidity (%)", x = "Date") + 
      scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                         sec.axis = sec_axis(trans = ~ .x / 5,
                                             name = "Precipitation (mm)",
                                             breaks = c(0,2,4,6,8,10,12,14,16,18,20))) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(name = "Data", values = Variable_cols[5:9])
    comp_RH[paste(n,"-",y,sep = "")] <- list(graph)
    rm(graph)
  }
}

#Compare iButton Temp with NASA T2M, T2M_Max and T2M_Min
comp_Temp <- list()
for (x in c("comp_data_2019","comp_data_2020")) {
  if (str_extract(string = x,pattern = "[0-9]{4}") == "2019") {
    y <- "2018-2019"
  }
  if (str_extract(string = x,pattern = "[0-9]{4}") == "2020") {
    y <- "2019-2020"
  }
  for (n in levels(get(x)$Zone)) {
    graph <- ggplot(get(x)[get(x)$Zone==n,], aes(YYYYMMDD)) + 
      stat_summary(aes(y = iTemp, col="iButton soil temp."), fun=mean,fun.args = list(na.rm = T), geom="line") + 
      stat_summary(aes(y = T2M, col="NASA air temp."), fun=mean,fun.args = list(na.rm = T), geom="line") +
      stat_summary(aes(y = Air.temp.avg., col="SASSCAL air temp."), fun=mean,fun.args = list(na.rm = T), geom="line") +
      stat_summary(aes(y = Soil.Temp.avg., col="SASSCAL soil temp."), fun=mean,fun.args = list(na.rm = T), geom="line") +
      scale_x_datetime(breaks = "2 month", date_labels = "%b %y") + 
      labs(title = paste("Temperature - ",n," (",y,")",sep = ""), y = "Temperature (°C)", x = "Date") + 
      scale_y_continuous(limits = c(0,45), breaks = c(0,5,10,15,20,25,30,35,40,45)) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(name = "Data", values = Variable_cols[1:4])
    comp_Temp[paste(n,"-",y,sep = "")] <- list(graph)
    rm(graph)
  }
}

#Nice comparison figure
NASA_iButton_comparison <- grid.arrange(comp_RH[[1]] + labs(tag = "A"),comp_RH[[3]] + labs(tag = "B"),comp_RH[[2]] + labs(tag = "C"),comp_RH[[4]] + labs(tag = "D"),comp_Temp[[1]] + labs(tag = "E"),comp_Temp[[3]] + labs(tag = "F"),comp_Temp[[2]] + labs(tag = "G"),comp_Temp[[4]] + labs(tag = "H"),ncol = 2)
ggsave(filename = "Figure 4.svg", NASA_iButton_comparison, width = 12, height = 9)

#Correlation table
corr_table <- as.data.frame(matrix(nrow = 10, ncol = 6,dimnames = list(NULL,c("Comparison","Overall (r2)","Fog 2018-2019 (r2)","Rain 2018-2019 (r2)","Fog 2019-2020 (r2)","Rain 2019-2020 (r2)"))))
corr_table[,"Comparison"] <- c("iRH vs RH2M","iRH vs Humidity","RH2M vs Humidity","iTemp vs T2M","iTemp vs Air.temp.avg.","iTemp vs Soil.Temp.avg.","T2M vs Air.temp.avg.","iRH vs Precipitation","iRH vs PRECTOTCORR","Precipitation vs PRECTOTCORR")
#overall
corr_table[1,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$iRH~Namib_Weather_all_data$RH2M))$r.squared
corr_table[2,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$iRH~Namib_Weather_all_data$Humidity))$r.squared
corr_table[3,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$RH2M~Namib_Weather_all_data$Humidity))$r.squared
corr_table[4,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$iTemp~Namib_Weather_all_data$T2M))$r.squared
corr_table[5,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$iTemp~Namib_Weather_all_data$Air.temp.avg.))$r.squared
corr_table[6,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$iTemp~Namib_Weather_all_data$Soil.Temp.avg.))$r.squared
corr_table[7,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$T2M~Namib_Weather_all_data$Air.temp.avg.))$r.squared
corr_table[8,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$iRH~Namib_Weather_all_data$Precipitation))$r.squared
corr_table[9,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$iRH~Namib_Weather_all_data$PRECTOTCORR))$r.squared
corr_table[10,"Overall (r2)"] <- summary(lm(Namib_Weather_all_data$Precipitation~Namib_Weather_all_data$PRECTOTCORR))$r.squared
#zones 2018-2019
for (x in unique(comp_data_2019$Zone)) {
  if (x == "Fog") {
    y <- "Fog 2018-2019 (r2)"
  }
  if (x == "Rain") {
    y <- "Rain 2018-2019 (r2)"
  }
  corr_table[1,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"iRH"]~comp_data_2019[comp_data_2019$Zone==x,"RH2M"]))$r.squared
  corr_table[2,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"iRH"]~comp_data_2019[comp_data_2019$Zone==x,"Humidity"]))$r.squared
  corr_table[3,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"RH2M"]~comp_data_2019[comp_data_2019$Zone==x,"Humidity"]))$r.squared
  corr_table[4,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"iTemp"]~comp_data_2019[comp_data_2019$Zone==x,"T2M"]))$r.squared
  corr_table[5,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"iTemp"]~comp_data_2019[comp_data_2019$Zone==x,"Air.temp.avg."]))$r.squared
  corr_table[6,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"iTemp"]~comp_data_2019[comp_data_2019$Zone==x,"Soil.Temp.avg."]))$r.squared
  corr_table[7,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"T2M"]~comp_data_2019[comp_data_2019$Zone==x,"Air.temp.avg."]))$r.squared
  corr_table[8,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"iRH"]~comp_data_2019[comp_data_2019$Zone==x,"Precipitation"]))$r.squared
  corr_table[9,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"iRH"]~comp_data_2019[comp_data_2019$Zone==x,"PRECTOTCORR"]))$r.squared
  corr_table[10,y] <- summary(lm(comp_data_2019[comp_data_2019$Zone==x,"Precipitation"]~comp_data_2019[comp_data_2019$Zone==x,"PRECTOTCORR"]))$r.squared
}
#zones 2019-2020
for (x in unique(comp_data_2020$Zone)) {
  if (x == "Fog") {
    y <- "Fog 2019-2020 (r2)"
  }
  if (x == "Rain") {
    y <- "Rain 2019-2020 (r2)"
  }
  corr_table[1,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"iRH"]~comp_data_2020[comp_data_2020$Zone==x,"RH2M"]))$r.squared
  corr_table[2,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"iRH"]~comp_data_2020[comp_data_2020$Zone==x,"Humidity"]))$r.squared
  corr_table[3,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"RH2M"]~comp_data_2020[comp_data_2020$Zone==x,"Humidity"]))$r.squared
  corr_table[4,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"iTemp"]~comp_data_2020[comp_data_2020$Zone==x,"T2M"]))$r.squared
  corr_table[5,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"iTemp"]~comp_data_2020[comp_data_2020$Zone==x,"Air.temp.avg."]))$r.squared
  corr_table[6,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"iTemp"]~comp_data_2020[comp_data_2020$Zone==x,"Soil.Temp.avg."]))$r.squared
  corr_table[7,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"T2M"]~comp_data_2020[comp_data_2020$Zone==x,"Air.temp.avg."]))$r.squared
  corr_table[8,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"iRH"]~comp_data_2020[comp_data_2020$Zone==x,"Precipitation"]))$r.squared
  corr_table[9,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"iRH"]~comp_data_2020[comp_data_2020$Zone==x,"PRECTOTCORR"]))$r.squared
  corr_table[10,y] <- summary(lm(comp_data_2020[comp_data_2020$Zone==x,"Precipitation"]~comp_data_2020[comp_data_2020$Zone==x,"PRECTOTCORR"]))$r.squared
}

write.table(corr_table,"Table S1.csv",row.names = F,col.names = T,quote = T,sep = ",")

##########

#Rainfall effect

#Can reduce the table to just the precipitation, iRH and iTemp
rain_effect <- Namib_Weather_all_data[,c("YYYYMMDD","Zone","Precipitation","iRH","iTemp")]
#Then average for each day for each zone
rain_effect_day <- data.frame()
for (x in unique(rain_effect$YYYYMMDD)) {
  for (y in unique(rain_effect$Zone)) {
    current <- data.frame()
    current[1,"YYYYMMDD"] <- as.POSIXct(x,origin = "1970-01-01")
    current[1,"Zone"] <- y
    current[1,"Precipitation"] <- mean(rain_effect[rain_effect$YYYYMMDD==x & rain_effect$Zone==y,"Precipitation"],na.rm = TRUE)
    current[1,"iRH"] <- mean(rain_effect[rain_effect$YYYYMMDD==x & rain_effect$Zone==y,"iRH"],na.rm = TRUE)
    current[1,"iTemp"] <- mean(rain_effect[rain_effect$YYYYMMDD==x & rain_effect$Zone==y,"iTemp"],na.rm = TRUE)
    rain_effect_day <- rbind(rain_effect_day,current)
  }
}
rain_effect_day$Zone <- factor(rain_effect_day$Zone, levels = c("Fog","Rain"), ordered = TRUE)

#Split the data between the two years
rain_effect_day_2019 <- rain_effect_day[rain_effect_day$YYYYMMDD >= "2018-04-08" & rain_effect_day$YYYYMMDD <= "2019-03-18",]
rain_effect_day_2020 <- rain_effect_day[rain_effect_day$YYYYMMDD >= "2019-04-30" & rain_effect_day$YYYYMMDD <= "2020-04-07",]

#Function for finding peaks in data
#https://github.com/stas-g/findPeaks 
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#Create a function that can find the rain events and avoid doubling up code
find_rain_events <- function (INPUT) {
  OUTPUT_FIG <- list()
  OUTPUT_TABLE <- data.frame()
  
  for (n in unique(INPUT$Zone)) {
    #rain peak position
    peak_pos <- find_peaks(INPUT[INPUT$Zone==n,"Precipitation"],m=4)
    
    #rain start positions
    start_pos <- vector()
    i <- 1
    for (x in peak_pos) {
      step <- 1
      done <- FALSE
      while (done == FALSE) {
        if (INPUT[INPUT$Zone==n,"Precipitation"][x - step] > 0) {
          step <- step + 1
        }
        else {
          done <- TRUE
        }
      }
      start_pos[i] <- x - step + 1
      i <- i+1
    }
    
    #
    
    #rain end positions
    end_pos <- vector()
    i <- 1
    for (x in peak_pos) {
      step <- 1
      done <- FALSE
      while (done == FALSE) {
        if (INPUT[INPUT$Zone==n,"Precipitation"][x + step] > 0) {
          step <- step + 1
        }
        else {
          done <- TRUE
        }
      }
      end_pos[i] <- x + step - 1
      i <- i+1
    }
    
    #amount of rain between these points
    precip_tot <- vector()
    i <- 1
    for (x in start_pos) {
      precip_tot[i] <- sum(INPUT[INPUT$Zone==n,"Precipitation"][x:end_pos[i]])
      i <- i + 1
    }
    
    #Mean iRH the preceding week
    mean_iRH <- vector()
    i <- 1
    for (x in start_pos) {
      if (x-7 > 0) {
        mean_iRH[i] <- mean(INPUT[INPUT$Zone==n,"iRH"][(x-7):x],na.rm = TRUE)
      }
      else {
        mean_iRH[i] <- mean(INPUT[INPUT$Zone==n,"iRH"][0:x],na.rm = TRUE)
      }
      i <- i+1
    }
    
    #Duration iRH is over preceding mean
    eff_end_pos <- vector()
    i <- 1
    for (x in end_pos) {
      step <- 1
      done <- FALSE
      while (done == FALSE) {
        if (INPUT[INPUT$Zone==n,"iRH"][x + step] > mean_iRH[i] & (x + step)<length(INPUT[INPUT$Zone==n,"iRH"])) {
          step <- step + 1
        }
        else {
          done <- TRUE
        }
      }
      eff_end_pos[i] <- x + step - 1
      i <- i+1
    }
    
    #Max iRH in that period
    max_iRH <- vector()
    i <- 1
    for (x in end_pos) {
      max_iRH[i] <- max(INPUT[INPUT$Zone==n,"iRH"][start_pos[i]:eff_end_pos[i]])
      i <- i + 1
    }
    
    #Construct a table
    current_rain_table <- data.frame(matrix(nrow = length(start_pos),dimnames = list(NULL,"Zone")))
    current_rain_table[,"Zone"] <- rep(n,length(start_pos))
    current_rain_table[,"Rain start"] <- INPUT[INPUT$Zone==n,"YYYYMMDD"][start_pos]
    current_rain_table[,"Rain end"] <- INPUT[INPUT$Zone==n,"YYYYMMDD"][end_pos]
    current_rain_table[,"Rain duration (days)"] <- end_pos - start_pos
    current_rain_table[,"Total precipitation (mm)"] <- precip_tot
    current_rain_table[,"Effect end"] <- INPUT[INPUT$Zone==n,"YYYYMMDD"][eff_end_pos]
    current_rain_table[,"Preceeding iRH (%)"] <- mean_iRH
    current_rain_table[,"Rain max iRH (%)"] <- max_iRH
    current_rain_table[,"Rain effect duration (days)"] <- eff_end_pos - end_pos
    current_rain_table[,"Rain effect magnitude (%)"] <- max_iRH - mean_iRH
    
    #Add drying gradient
    for (z in 1:length(start_pos)) {
      valx <- (start_pos[z]:eff_end_pos[z])
      valy <- (INPUT[INPUT$Zone==n,"iRH"][start_pos[z]:eff_end_pos[z]])
      current_rain_table[z,"Drying rate (%/day)"] <- coef(lm(valy~valx))[2]
    }
    
    #Merge the tables
    OUTPUT_TABLE <- rbind(OUTPUT_TABLE,current_rain_table)
    
    #One of the rain events ends up duplicated. Needs to be trimmed out.
    OUTPUT_TABLE <- OUTPUT_TABLE[!duplicated(OUTPUT_TABLE),]
    
    #Plot the rain events
    current_rain_plot <- ggplot(INPUT[INPUT$Zone==n,], aes(YYYYMMDD)) + 
      geom_line(aes(y = iRH), colour = "blue") + 
      geom_line(aes(y = Precipitation * 5)) + 
      scale_x_datetime(breaks = "2 month", date_labels = "%b %y") + 
      labs(title = paste("Rain events - ",n," Zone (",paste(unique(year(INPUT$YYYYMMDD))[1],unique(year(INPUT$YYYYMMDD))[2],sep="-"),")",sep = ""), y = "Relative Humidity (%)", x = "Date") + 
      scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100),
                         sec.axis = sec_axis(trans = ~ .x / 5,
                                             name = "Precipitation (mm)",
                                             breaks = c(0,2,4,6,8,10,12,14,16,18,20))) +
      geom_vline(xintercept = INPUT[INPUT$Zone==n,"YYYYMMDD"][peak_pos],colour="blue", alpha = 0.6) +
      geom_vline(xintercept = INPUT[INPUT$Zone==n,"YYYYMMDD"][eff_end_pos],colour="red", alpha = 0.6) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
    
    OUTPUT_FIG[n] <- list(current_rain_plot)
    
    
  }
  
  #Rain effect duration vs precipitation
  OUTPUT_PRECIPDURA <- ggplot(OUTPUT_TABLE, aes_(x = as.name("Total precipitation (mm)"),y = as.name("Rain effect duration (days)"))) +
    geom_point(aes(col=Zone)) +
    geom_smooth(method = "lm", alpha = 0.3, colour = "black") +
    scale_x_continuous(limits = c(0,15), breaks = c(0,5,10,15)) +
    scale_y_continuous(limits = c(0,50), breaks = c(0,10,20,30,40,50)) +
    labs(title = paste("Total precipitation versus effect duration (",paste(unique(year(INPUT$YYYYMMDD))[1],unique(year(INPUT$YYYYMMDD))[2],sep="-"),")",sep = ""), y = "Effect duration (days)", x = "Total precipitation (mm)") + 
    scale_colour_manual(values = zone_cols) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  #Drying rate graph
  OUTPUT_DRYINGRATE <- ggplot(OUTPUT_TABLE, aes_(y = as.name("Drying rate (%/day)"), x = as.name("Zone"))) +
    geom_boxplot(aes(fill=Zone)) +
    labs(title = "Drying rates in the different zones", y = "Drying rate (%RH/day)") + 
    scale_fill_manual(values = zone_cols) +
    geom_signif(comparisons = list(c("Fog","Rain")), map_signif_level=TRUE, colour="black") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  #Result
  result <- list(OUTPUT_TABLE,OUTPUT_FIG,OUTPUT_PRECIPDURA,OUTPUT_DRYINGRATE)
  return(result)
}

#2019 Rain events
rain_res_2019 <- find_rain_events(rain_effect_day_2019)
rain_res_2020 <- find_rain_events(rain_effect_day_2020)

#Save our rain data
rain_fig <- grid.arrange(rain_res_2019[[2]][[2]] + labs(tag = "A"),rain_res_2019[[2]][[1]] + labs(tag = "B"),rain_res_2019[[3]] + labs(tag = "C"),rain_res_2020[[2]][[2]] + labs(tag = "D"),rain_res_2020[[2]][[1]] + labs(tag = "E"),rain_res_2020[[3]] + labs(tag = "F"),nrow = 2)
ggsave(filename = "Figure 5.svg", rain_fig, width = 16, height = 9)

rain_event_table <- rbind(rain_res_2019[[1]],rain_res_2020[[1]])
rain_event_table$Year <- c(rep("2018-2019",nrow(rain_res_2019[[1]])),rep("2019-2020",nrow(rain_res_2020[[1]])))
write.table(rain_event_table,"Table S2.csv",row.names = F,col.names = T,quote = T,sep = ",")

##########

#Water activity table

#Set up the table and list the sites, processes and their water activity limits
water_activity_table <- as.data.frame(matrix(nrow = 1, ncol = 4, dimnames = list(NULL,c("Site","Process","Hours","Year"))))
sites <- c(2,4,6,8,10,12,14,16,18,20)
processes <- c("Microbial community growth","Photosynthesis","Amino acid metabolism","Basal metabolism")
limits <- c(90,80,30,20)

#Loop through the variables to build the table
#Count how many readings are above the threshold per site and multiply by 4 to get the total number of hours
for (SITE in sites) {
  for (YEAR in unique(iButton_complete$Year)) {
    for (n in 1:length(processes)) {
      water_activity_table <- rbind(water_activity_table,
            c(SITE,
              processes[n],
              round((nrow(iButton_complete[iButton_complete$Site == SITE & iButton_complete$Year == YEAR & iButton_complete$RH >= limits[n],]))*4),
              YEAR))
    }
  }
}

#Format the table
water_activity_table <- water_activity_table[-1,]
water_activity_table$Site <- as.numeric(water_activity_table$Site)
water_activity_table$Process <- factor(water_activity_table$Process,levels = c("Microbial community growth","Photosynthesis","Amino acid metabolism","Basal metabolism"))
water_activity_table$Hours <- as.numeric(water_activity_table$Hours)
water_activity_table$Year <- as.numeric(water_activity_table$Year)

#Remove the missing, data 2020 site 14, which calls incorrect row numbers
water_activity_table <- subset(water_activity_table, !(water_activity_table$Site==14 & water_activity_table$Year==2020))

#Will need an inflation parameter to compensate for not having readings for every single day of the year.
#Inflation parameter is total hours in an average year, divided by the hours measured in the year
total_time_year <- 365.25*24
total_time_2019 <- as.numeric(difftime(max(iButton_complete[iButton_complete$Year==2019,"Date"]),min(iButton_complete[iButton_complete$Year==2019,"Date"]),units = "hours"))
total_time_2020 <- as.numeric(difftime(max(iButton_complete[iButton_complete$Year==2020,"Date"]),min(iButton_complete[iButton_complete$Year==2020,"Date"]),units = "hours"))
water_activity_table[which(water_activity_table$Year == 2019),"Hours"] <- water_activity_table[which(water_activity_table$Year == 2019),"Hours"] * (total_time_year/total_time_2019)
water_activity_table[which(water_activity_table$Year == 2020),"Hours"] <- water_activity_table[which(water_activity_table$Year == 2020),"Hours"] * (total_time_year/total_time_2020)
water_activity_table[,"%Year"] <- water_activity_table[,"Hours"]/8766 * 100

#Graph the results
metabolic_windows <- ggplot(water_activity_table,aes(x = Site, y = Hours)) +
  annotate("rect", xmin=2, xmax=6, ymin=0, ymax=8766, alpha=0.2, fill="light blue") +
  geom_line(data = water_activity_table[water_activity_table$Year==2019,],aes(x = Site, y = Hours, colour = Process)) +
  geom_line(data = water_activity_table[water_activity_table$Year==2020,],aes(x = Site, y = Hours, colour = Process),linetype = "dashed") +
  scale_y_continuous(limits = c(0,8766), breaks = c(0,877,1753,2630,3506,4383,5260,6136,7013,7889,8766), 
                     sec.axis = sec_axis(trans = ~ .x / 8766 * 100,
                                         name = "Percent of year",
                                         breaks = c(0,10,20,30,40,50,60,70,80,90,100))) +
  scale_x_continuous(limits = c(2,20), breaks = c(2,4,6,8,10,12,14,16,18,20)) +
  labs(title = "Annual Metabolic Windows (1-2 cm deep)",colour="") + 
  theme_minimal() +
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values = process_cols) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE))

ggsave(filename = "Figure 6.svg", metabolic_windows, width = 8, height = 4.5)

#Export the table
write.table(water_activity_table,"water_activity_table.csv",row.names = F,col.names = T,quote = T,sep = ",")

