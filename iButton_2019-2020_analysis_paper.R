#Namib iButton Data Analysis
#This script analyses iButton data collected from the Namib Desert as well as data from SASSCAL and NASA to make inferences
#about soil moisture and microbial activity for the 2019-2020 period. 
#While all the analysis for the 2018-2019 period is possible, the only part used in the publication was the water activity 
#table and total rainfall. Other analyses that were performed for the 2018-2019 data but not published for the 2019-2020 data
#are commented out.
#Written by: Jason Bosch
#Last update: 15/06/2021

###############
####SET UP#####
###############

##The locations for the results folder and the raw data that needs to be imported will have to be adjusted to suit your setup.##

#Load the required libraries
library("ggplot2")
library("gridExtra")
library("ggsignif")
library("lubridate")
library("nasapower")
library("stringr")
library("vegan")

#Set working directory
setwd("~/PostDoc/02_Namib_10yr_project/iButton_data/2020_analysis/")

###Original Data Preparation###

#This step only needed to be run once to bring all the data into a single file.
#It is unnecessary for the published iButton data and is included for reference purposes only.
#It is necessary for retrieving/preparing the NASA data and for preparing the SASSCAL data which must be downloaded
#from their website.
#Each section can be uncommented as necessary.

#iButton preparation

# #Bring in data
# sites <- c("2","4","6","8","10","12","14","16","18","20")
# 
# #Have to use try() for everything because some data is missing and this will stop everything from grinding to a halt.
# 
# #2020
# for (x in sites) {
#   try(import_rh <- read.csv(paste("../C14_Transect_iButtons_2020/TR-",x,"_RH_Apr_2021.csv",sep = ""),stringsAsFactors = F,skip = 19))
#   try(import_temp <- read.csv(paste("../C14_Transect_iButtons_2020/TR-",x,"_T_Apr_2021.csv",sep = ""),stringsAsFactors = F,skip = 19))
#   try(import_file <- merge.data.frame(import_rh, import_temp, by = "Date.Time"))
#   try(for (n in 1:nrow(import_file)) {
#     import_file[n,6] <- as.POSIXct(strptime(import_file[n,1],format = "%d/%m/%y %H:%M:%S"))
#   })
#   if (nchar(x) == 1) {
#     try(assign(paste("iButton_2020_0",x,sep = ""),import_file))
#   }
#   else {
#     try(assign(paste("iButton_2020_",x,sep = ""),import_file))
#   }
#   try(rm(import_file))
#   try(rm(import_temp))
#   try(rm(import_rh))
# }
# 
# #Create an empty blank file for samples with missing data to fill in graphing gap from missing site 14
# 
# iButton_2020_14 <- iButton_2020_02
# iButton_2020_14[,3] <- NaN
# iButton_2020_14[,5] <- NaN
# 
# #Create a single data frame with all the data
# sites <- c("02","04","06","08","10","12","14","16","18","20")
# iButton_complete <- data.frame()
# for (x in 1:length(sites)) {
#   current <- get(paste("iButton_2020_",sites[x],sep = ""))
#   current[,7] <- sites[x]
#   current[,8] <- "2020"
#   current[,9] <- paste("2020_",sites[x],sep="")
#   iButton_complete <- rbind(iButton_complete,current)
# }
# 
# #Include the rainfall zones
# iButton_complete[iButton_complete[,7]%in%c("02","04","06"),10] <- "Fog"
# iButton_complete[iButton_complete[,7]%in%c("08","10","12","14"),10] <- "Low rain"
# iButton_complete[iButton_complete[,7]%in%c("16","18","20"),10] <- "High rain"
# 
# #Label all the columns
# colnames(iButton_complete) <- c("date.time","unit.RH","RH","unit.Temp","Temp","Date","Site","Year","Group","Zone")
# 
# #Write the data to simplify everything else that will be done later
# write.csv(iButton_complete,"iButton_Complete_dataset",row.names = F)

#NASA preparation

# #Bring in the data for all the sites and for measurements of interest
# 
# # PRECTOT (Precipitation, mm day-1)
# # QV2M (Specific Humidity at 2 Meters, kg kg-1)
# # RH2M (Relative Humidity at 2 Meters, %)
# # T2M (Temperature at 2 Meters, C)
# # T2M_MAX (Maximum Temperature at 2 Meters, C)
# # T2M_MIN (Minimum Temperature at 2 Meters, C)
# # WS2M (Wind Speed at 2 Meters, m/s)
# 
# #This copies co-ords for the sites but they are too close together to get unique data.
# #Since we're averaging per zone, it's probably not a disaster
# #Need to increase the first site's latitude by 0.01 as it is on the border of block. This keeps it with the second site.
# coords <- as.data.frame(matrix(nrow = 10,ncol = 3,dimnames = list(NULL,c("Site","Lat","Lon"))))
# coords$Site <- c("2","4","6","8","10","12","14","16","18","20")
# coords$Lat <- c("-23.01","-23.02","-23.07","-23.14","-23.25","-23.31","-23.33", "-23.32","-23.35","-23.24")
# coords$Lon <- c("14.67","14.86","15.04","15.21","15.36","15.53","15.71","15.86","16.01","16.14")
# 
# NASA_data <- data.frame()
# for (x in 1:nrow(coords)) {
# current_data <- get_power(community = "AG",
#                           pars = c("PRECTOT","QV2M","RH2M","T2M","T2M_MAX","T2M_MIN","WS2M"),
#                           temporal_average = "DAILY",
#                           lonlat = c(as.numeric(coords[x,"Lon"]),as.numeric(coords[x,"Lat"])),
#                           dates = c("2019-05-01","2020-04-06"))
#   NASA_data <- rbind(NASA_data,current_data)
# }
# 
# write.csv(NASA_data,"NASA_data",row.names = F)
# 
# #To remove any weird formatting of the data frame
# NASA_data <- read.csv("NASA_data")
# 
# #Since the coords got slightly off we need to round them first
# NASA_data$LON <- signif(NASA_data$LON,4)
# NASA_data$LAT <- signif(NASA_data$LAT,4)
# 
# #Need to make sure NASA data mentions the sites correctly
# for (x in 1:nrow(coords)) {
#   NASA_data[NASA_data$LAT == as.numeric(coords$Lat[x]), "Site"] <- coords$Site[x]
# }
# 
# #Find everything in iButton data with the same data and average the iButton data for each day
# #Slow
# 
# for (x in 1:nrow(NASA_data)) {
#   NASA_data[x,c("iRH","iTemp")] <- colMeans(iButton_complete[iButton_complete$Site == NASA_data[x,"Site"] & str_extract(iButton_complete$Date, "^.{10}") == NASA_data[x,"YYYYMMDD"],c("RH","Temp")])
# }
# 
# #Add the zone markers
# NASA_data[NASA_data[,"Site"]%in%c("2","4","6"),"Zone"] <- "Fog"
# NASA_data[NASA_data[,"Site"]%in%c("8","10","12","14"),"Zone"] <- "Low rain"
# NASA_data[NASA_data[,"Site"]%in%c("16","18","20"),"Zone"] <- "High rain"
# 
# write.csv(NASA_data,"NASA_iButton_merged_data",row.names = F)

#SASSCAL preparation

# #Bring in the data files
# #Marble Koppie, Gobabeb Met and Narais-Duruchaus will be excluded as they are located far from the sites
# CoastalMet <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_CoastalMet_E7631_CRNB.csv",sep=";",comment.char = "",header = T)
# Kleinberg <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_Kleinberg-FN_E7630_CRNB.csv",sep=";",comment.char = "",header = T)
# SophiesHoogte <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_SophiesHoogte_E7629_CRNB.csv",sep=";",comment.char = "",header = T)
# Vogelfederberg <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_Vogelfederberg_E7626_CRNB.csv",sep=";",comment.char = "",header = T)
# GarnetKoppie <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_GarnetKoppie_E7628_CRNB.csv",sep=";",comment.char = "",header = T)
# Rooisand <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_Rooisand_103_MBBW.csv",sep=";",comment.char = "",header = T)
# 
# #Keep only columns of interest which are also present in all datasets
# CoastalMet <- subset(CoastalMet,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Wind.speed..vc.avg.,Wind.speed..max.,Soil.Temp.1..avg.))
# Kleinberg <- subset(Kleinberg,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Wind.speed..vc.avg.,Wind.speed..max.,Soil.Temp.1..avg.))
# SophiesHoogte <- subset(SophiesHoogte,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Wind.speed..vc.avg.,Wind.speed..max.,Soil.Temp.1..avg.))
# Vogelfederberg <- subset(Vogelfederberg,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Wind.speed..vc.avg.,Wind.speed..max.,Soil.Temp.1..avg.))
# GarnetKoppie <- subset(GarnetKoppie,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Wind.speed..vc.avg.,Wind.speed..max.,Soil.Temp.1..avg.))
# Rooisand <- subset(Rooisand,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Wind.speed..vc.avg.,Wind.speed..max.,Soil.Temp.1..avg.))
# 
# #Remove the top line which has the units
# CoastalMet <- CoastalMet[-1,]
# Kleinberg <- Kleinberg[-1,]
# SophiesHoogte <- SophiesHoogte[-1,]
# Vogelfederberg <- Vogelfederberg[-1,]
# GarnetKoppie <- GarnetKoppie[-1,]
# Rooisand <- Rooisand[-1,]
# 
# #For some reason Rooisand has data for points which are missing in all other data sets. This causes downstream problems.
# #Simplest solution is to remove those data points so that everything matches.
# Rooisand <- subset(Rooisand,Rooisand$Date%in%CoastalMet$Date)
# 
# #Combine all the data
# weather_station_complete <- Reduce(rbind,list(CoastalMet,Kleinberg,SophiesHoogte,Vogelfederberg,GarnetKoppie,Rooisand))
# 
# #Some stations have missing values which need to be removed
# for (x in 1:ncol(weather_station_complete)) {
#   weather_station_complete[which(weather_station_complete[,x]==" - "),x] <- NA
# }
# 
# #Replace with nicer column names
# colnames(weather_station_complete) <- c("Station","Date","Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.")
# 
# #Include the rainfall zones
# weather_station_complete[weather_station_complete[,"Station"]%in%c("Coastal Met","Kleinberg-FN","Sophies Hoogte","Vogelfederberg"),"Zone"] <- "Fog"
# weather_station_complete[weather_station_complete[,"Station"]%in%c("Garnet Koppie"),"Zone"] <- "Low rain"
# weather_station_complete[weather_station_complete[,"Station"]%in%c("Rooisand"),"Zone"] <- "High rain"
# 
# #Create a new POSIX compliant date column.
# for (n in 1:nrow(weather_station_complete)) {
#   weather_station_complete[n,"Date2"] <- as.POSIXct(weather_station_complete[n,"Date"],format = "%Y-%m-%d")
# }
# weather_station_complete["Date"] <- weather_station_complete["Date2"]
# weather_station_complete <- subset(weather_station_complete, select = -Date2)
# 
# #Save the new file for later use
# write.csv(weather_station_complete,"weather_station_complete",row.names = F)
# 
# #Bring in data we need (also changes the column types from character)
# weather_station_complete <- read.csv("weather_station_complete",stringsAsFactors = F,colClasses = c(Date="POSIXct"))
# 
# #The weather station data will be assigned/averaged for the closest site in its zone.
# #Site 2 will get the average of coastal Met and Kleinberg
# #Site 4 will get Sophies Hoogte
# #Site 6 will get Vogelfederberg
# #All low rain zone will get Garnet Koppie
# #All high rain zone will get Rooisand
# 
# #To merge it with the other data, we will need to summarise per zone.
# weather_station_averaged <- data.frame()
# for (x in 1:length(unique(weather_station_complete$Date))) {
#     current_res <- data.frame()
#     current_res[1:10,"Date"] <- unique(weather_station_complete$Date)[x]
#     current_res[,"Site"] <- c("2","4","6","8","10","12","14","16","18","20")
#     current_res[1,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.","Zone")] <-
#       c(colMeans(weather_station_complete[weather_station_complete$Station%in%c("Coastal Met","Kleinberg-FN") & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.")],na.rm = T),"Fog")
#     current_res[2,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.","Zone")] <-
#       c(weather_station_complete[weather_station_complete$Station=="Sophies Hoogte" & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.")],"Fog")
#     current_res[3,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.","Zone")] <-
#       c(weather_station_complete[weather_station_complete$Station=="Vogelfederberg" & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.")],"Fog")
#     current_res[4:7,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.","Zone")] <-
#       c(weather_station_complete[weather_station_complete$Station=="Garnet Koppie" & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.")],"Low rain")
#     current_res[8:10,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.","Zone")] <-
#       c(weather_station_complete[weather_station_complete$Station=="Rooisand" & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Wind.speed.avg.","Wind.speed.max.","Soil.Temp.avg.")],"High rain")
#     weather_station_averaged <- rbind(weather_station_averaged,current_res)
# }
# 
# #Write the averaged data
# write.csv(weather_station_averaged,"weather_station_averaged",row.names = F)

###Import the cleaned data###

#Bring in the data
iButton_complete <- read.csv("iButton_Complete_dataset",stringsAsFactors = F,colClasses = c(Date="POSIXct",Year="character"))

#Makes sure everything stays in the right order instead of alphabetical
iButton_complete[,9] <- factor(iButton_complete[,9], levels = unique(iButton_complete[,9]))
iButton_complete[,10] <- factor(iButton_complete[,10], levels = unique(iButton_complete[,10]))

#Remove the RH values which are out of bounds and do not make sense
iButton_complete[which(iButton_complete[,3] > 100),3] <- 100
iButton_complete[which(iButton_complete[,3] < 0),3] <- NaN

#Need the hour as a separate column where all the dd/MM/YYYY info is gone, in this case all converted to the day of the command
#Has to be done on the day that the script is run or the limits of some graphs will break.
Times <- strftime(iButton_complete[,6], format="%H:%M:%S")
Times <- as.POSIXct(Times, format="%H:%M:%S")
iButton_complete[,"Time"] <- Times

#Bring in NASA data
NASA_data <- read.csv("NASA_iButton_merged_data",colClasses = c(YYYYMMDD="POSIXct"))

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
weather_station_complete <- read.csv("weather_station_complete",stringsAsFactors = F,colClasses = c(Date="POSIXct"))
weather_station_averaged <- read.csv("weather_station_averaged",stringsAsFactors = F,colClasses = c(Date="POSIXct"))

#Set the season (Determined via PCA plot of relevant values)
weather_station_complete[month(weather_station_complete[,"Date"])%in%c(1,2,3,4,10,11,12),"Season"] <- "Summer"
weather_station_complete[month(weather_station_complete[,"Date"])%in%c(5,6,7,8,9),"Season"] <- "Winter"

#Set the zones as factors for correct ordering in graphs
weather_station_complete$Zone <- factor(weather_station_complete$Zone,levels = c("Fog","Low rain", "High rain"), ordered = T)
weather_station_averaged$Zone <- factor(weather_station_averaged$Zone,levels = c("Fog","Low rain", "High rain"), ordered = T)

#Combine all the data sets
Namib_Weather_all_data <- merge(NASA_data,weather_station_averaged, by.x = c("YYYYMMDD","Site"), by.y = c("Date","Site"))
#Clean up one of the zone names to use in analyses
colnames(Namib_Weather_all_data)[28] <- "Zone"

###################
#####ANALYSIS######
###################

#Setting constant colours
#colourbrewer2 for initial colour selection: https://colorbrewer2.org
#Canva's colour wheel for Triadic colour: https://www.canva.com/colors/color-wheel/ 

Zone_cols <- c("#31a354","#756bb1","#e6550d","#5431A3","#B1756B","#0DE655","#A35431","#6BB175","#550DE6")
names(Zone_cols) <- c("Fog","Low rain","High rain","NASA Fog","NASA Low rain","NASA High rain", "Station Fog", "Station Low rain", "Station High rain")

Variable_cols <- c("#e66101","#5e3c99","#b2abd2","#fdb863","#b2df8a","#1f78b4","#a6cee3","red","black")
names(Variable_cols) <- c("iButton soil temp.","NASA air temp.","SASSCAL air temp.","SASSCAL soil temp.","iButton soil RH (%)","NASA air RH (%)","SASSCAL air RH (%)","SASSCAL rain (mm)","NASA rain (mm)")

##########

# #Graph the data per year for all sites
# 
# sites <- c("02","04","06","08","10","12","14","16","18","20")
# raw_graphs <- list()
# for (n in 1:length(sites)) {
#   graph <- ggplot(iButton_complete[iButton_complete$Group==paste("2020_",sites[n],sep=""),], aes(Date)) + 
#     geom_line(aes(y = RH), color = "blue", size = .4) + geom_line(aes(y = Temp), color = "red", size = .4) + 
#     scale_x_datetime(breaks = "1 month", minor_breaks = "1 week", date_labels = "%b-%y") + 
#     labs(title = paste("C14_",sites[n],"_2020",sep=""),  y = NULL, x = NULL) + 
#     scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
#     theme(axis.text.x=element_text(angle=15, hjust=1)) + 
#     theme_minimal()
#   raw_graphs[paste("graph_2020_",sites[n],sep = "")] <- list(graph)
#   rm(graph)
# }
# 
# raw_image <- grid.arrange(raw_graphs[[1]] + labs(tag = "A"),raw_graphs[[2]] + labs(tag = "B"),raw_graphs[[3]] + labs(tag = "C"),raw_graphs[[4]] + labs(tag = "D"),raw_graphs[[5]] + labs(tag = "E"),raw_graphs[[6]] + labs(tag = "F"),raw_graphs[[7]] + labs(tag = "G"),raw_graphs[[8]] + labs(tag = "H"),raw_graphs[[9]] + labs(tag = "I"),raw_graphs[[10]] + labs(tag = "J"),nrow = 5)
# ggsave(filename = "Raw_readings.svg", raw_image, width = 16, height = 9, dpi = 600)
# ggsave(filename = "Raw_readings.jpg", raw_image, width = 16, height = 9, dpi = 100)

##########

# #Using PCA to find/confirm the data groupings
# 
# #PCA (Can not use humidity because it was missing from Rooisand data)
# PCA_test <- subset(Namib_Weather_all_data, select = c(Site,Zone,MM,Season,PRECTOT,QV2M,RH2M,T2M,T2M_MAX,T2M_MIN,iTemp,iRH,Air.temp.avg.,Air.temp.min.,Air.temp.max.,Precipitation,Soil.Temp.avg.))
# PCA_test <- na.omit(PCA_test)
# test <- subset(PCA_test, select = c(PRECTOT,QV2M,RH2M,T2M,T2M_MAX,T2M_MIN,iTemp,iRH,Air.temp.avg.,Air.temp.min.,Air.temp.max.,Precipitation,Soil.Temp.avg.))
# pca_res <- prcomp(test, scale. = TRUE)
# dtp <- data.frame('Site' = as.character(PCA_test$Site),'Zone' = as.character(PCA_test$Zone), 'Month' = as.character(PCA_test$MM),'Season' = as.character(PCA_test$Season), pca_res$x[,1:2]) 
# PCAloadings <- data.frame(Variables = rownames(pca_res$rotation), pca_res$rotation)
# 
# #Keep the ordering
# dtp$Zone <- factor(dtp$Zone, levels = unique(dtp$Zone), ordered = TRUE)
# dtp$Month <- factor(dtp$Month, levels = sort(as.numeric(unique(dtp$Month))), ordered = TRUE)
# dtp$Site <- factor(dtp$Site, levels = sort(as.numeric(unique(dtp$Site))), ordered = TRUE)
# 
# #Plot sites
# site_graph <- ggplot(data = dtp) + 
#   geom_point(aes(x = PC1, y = PC2, col = Site)) + 
#   theme_minimal() + 
#   scale_colour_manual(values=c("#74c476","#31a354","#006d2c","#cbc9e2","#9e9ac8","#756bb1","#54278f","#fd8d3c","#e6550d","#a63603")) + 
#   labs(title ="Sites")
# site_graph_loadings <- site_graph + 
#   geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), arrow = arrow(length = unit(1/2, "picas")), color = "black") + 
#   annotate("text", x = (PCAloadings$PC1*8), y = (PCAloadings$PC2*8), label = PCAloadings$Variables)
# 
# #plot zones
# zone_graph <- ggplot(data = dtp) + 
#   geom_point(aes(x = PC1, y = PC2, col = Zone)) + 
#   theme_minimal() + 
#   scale_colour_manual(values = Zone_cols) + 
#   labs(title ="Zones") 
# zone_graph_loadings <- zone_graph + 
#   geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), arrow = arrow(length = unit(1/2, "picas")), color = "black") + 
#   annotate("text", x = (PCAloadings$PC1*8), y = (PCAloadings$PC2*8), label = PCAloadings$Variables)
# 
# #Plot months
# month_graph <- ggplot(data = dtp) + 
#   geom_point(aes(x = PC1, y = PC2, col = Month)) + theme_minimal() +
#   #scale_color_manual(values=c("red","red","red","red","blue","blue","blue","blue","blue","red","red","red")) +
#   scale_color_brewer(type = "qual", palette = "Set3") + 
#   labs(title ="Months")
# month_graph_loadings <- month_graph + 
#   geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), arrow = arrow(length = unit(1/2, "picas")), color = "black") + 
#   annotate("text", x = (PCAloadings$PC1*8), y = (PCAloadings$PC2*8), label = PCAloadings$Variables)
# 
# #Putative seasons based on how the data is sorted in a PCA
# #Summer: 1, 2, 3, 4, 10, 11, 12
# #Winter: 5, 6, 7, 8, 9
# 
# #Plot seasons
# season_graph <- ggplot(data = dtp) + 
#   geom_point(aes(x = PC1, y = PC2, col = Season)) + 
#   theme_minimal() + 
#   scale_colour_manual(values=c("Red","Blue")) + 
#   labs(title ="Determined Seasons")
# season_graph_loadings <- season_graph + 
#   geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), arrow = arrow(length = unit(1/2, "picas")), color = "black") + 
#   annotate("text", x = (PCAloadings$PC1*8), y = (PCAloadings$PC2*8), label = PCAloadings$Variables)
# 
# #Create a combined image
# data_grouping_image <- grid.arrange(site_graph,site_graph_loadings,zone_graph,zone_graph_loadings,month_graph,month_graph_loadings,season_graph,season_graph_loadings,ncol=2)
# ggsave(filename = paste("data_grouping_image.png",sep=""), data_grouping_image, width = 9, height = 16, dpi = 600)
# ggsave(filename = paste("data_grouping_image.jpg",sep=""), data_grouping_image, width = 9, height = 16, dpi = 100)
# #SVG version has too many points so is very slow to load.
# #ggsave(filename = paste("data_grouping_image.svg",sep=""), data_grouping_image, width = 9, height = 16, dpi = 600)
# 
# #Create a combined image
# smaller_group <- grid.arrange(zone_graph + labs(tag = "A"),zone_graph_loadings + labs(tag = "B"),season_graph + labs(tag = "C"),season_graph_loadings + labs(tag = "D"),ncol=2)
# ggsave(filename = paste("data_grouping_image_small.png",sep=""), smaller_group, width = 16, height = 9, dpi = 600)
# ggsave(filename = paste("data_grouping_image_small.jpg",sep=""), smaller_group, width = 16, height = 9, dpi = 200)
# 
# #Cluster according to euclidean distance
# PCAcoords <- dtp[,c("PC1","PC2")]
# PCAdist <- dist(PCAcoords)
# 
# #Adonis shows different centroids for the zones
# adonis(PCAdist ~ Zone, data = dtp, method='eu')
# # Call:
# # adonis(formula = PCAdist ~ Zone, data = dtp, method = "eu") 
# # 
# # Permutation: free
# # Number of permutations: 999
# # 
# # Terms added sequentially (first to last)
# # 
# #             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# # Zone         2    4170.9 2085.44  254.66 0.14282  0.001 ***
# # Residuals 3057   25033.8    8.19         0.85718           
# # Total     3059   29204.7                 1.00000           
# # ---
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# #However we do see differences in beta dispersal
# beta <- betadisper(PCAdist, dtp$Zone)
# permutest(beta)
# # Permutation test for homogeneity of multivariate dispersions
# # Permutation: free
# # Number of permutations: 999
# # 
# # Response: Distances
# #             Df Sum Sq Mean Sq      F N.Perm Pr(>F)   
# # Groups       2   30.1 15.0340 7.6607    999  0.002 **
# # Residuals 3057 5999.3  1.9625                        
# # ---
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# TukeyHSD(beta)
# # Tukey multiple comparisons of means
# # 95% family-wise confidence level
# # 
# # Fit: aov(formula = distances ~ group, data = df)
# # 
# # $group
# #                           diff         lwr        upr     p adj
# # Low rain-Fog       -0.07523578 -0.22047817 0.07000661 0.4445513
# # High rain-Fog       0.16257770  0.01701339 0.30814202 0.0240661
# # High rain-Low rain  0.23781348  0.09224917 0.38337780 0.0003829
# adonis(PCAdist ~ Season, data = dtp, method='eu')
# # Call:
# # adonis(formula = PCAdist ~ Season, data = dtp, method = "eu") 
# # 
# # Permutation: free
# # Number of permutations: 999
# # 
# # Terms added sequentially (first to last)
# # 
# #             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# # Season       1    8162.5  8162.5  1186.2 0.27949  0.001 ***
# # Residuals 3058   21042.1     6.9         0.72051           
# # Total     3059   29204.7                 1.00000           
# # ---
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##########

# #Daily graphs
# 
# #Smoothed curve uses the standard setting: geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")
# #Non-smoothed follows the smoothed line very well for all except fog RH which is too jerky. Chose to use smoothed for all
# 
# #To avoid weird artifacts in the plotting, the data is replicated over three days but zoomed in to only show one.
# 
# #Preparing the data
# Daily_df1 <- iButton_complete[,c("Time","RH","Temp","Zone","Season")]
# Daily_df2 <- iButton_complete[,c("Time","RH","Temp","Zone","Season")]
# day(Daily_df2$Time) <- day(Daily_df2$Time) - 1
# Daily_df3 <- iButton_complete[,c("Time","RH","Temp","Zone","Season")]
# day(Daily_df3$Time) <- day(Daily_df3$Time) + 1
# 
# daily_full <- rbind.data.frame(Daily_df1,Daily_df2,Daily_df3)
# 
# #plotting
# 
# #Combined
# daily_graphs <- list()
# for (n in 1:length(levels(daily_full$Zone))) {
#   graph <- ggplot(daily_full[daily_full[,"Zone"]==levels(daily_full$Zone)[n],], aes(Time)) + 
#     geom_point(aes(y = RH), color = "light blue", alpha = 0.25, size = .4, position = "jitter") + 
#     geom_smooth(aes(y = RH), color = "blue",) + 
#     geom_point(aes(y = Temp), color = "pink", alpha = 0.25,size = .4, position = "jitter") + 
#     geom_smooth(aes(y = Temp), color = "red",) +
#     scale_x_datetime(breaks = "2 hour", minor_breaks = "1 hour", date_labels = "%H:%M") +
#     coord_cartesian(xlim = (c(as.POSIXct("00:00:00", format="%H:%M:%S"),as.POSIXct("23:59:59", format="%H:%M:%S")))) +
#     labs(title = paste(levels(daily_full$Zone)[n]," zone - Annual",sep = ""), y = NULL, x = "Time") + 
#     scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
#     theme_minimal()
#   daily_graphs[paste("Daily_",gsub(" ","_",levels(daily_full$Zone)[n]),"_Zone",sep = "")] <- list(graph)
#   rm(graph)
# }
# 
# #Summer
# summer_graphs <- list()
# for (n in 1:length(levels(daily_full$Zone))) {
#   graph <- ggplot(daily_full[daily_full[,"Zone"]==levels(daily_full$Zone)[n] & daily_full$Season=="Summer",], aes(Time)) + 
#     geom_point(aes(y = RH), color = "light blue", alpha = 0.25, size = .4, position = "jitter") + 
#     geom_smooth(aes(y = RH), color = "blue",) + 
#     geom_point(aes(y = Temp), color = "pink", alpha = 0.25,size = .4, position = "jitter") + 
#     geom_smooth(aes(y = Temp), color = "red",) +
#     scale_x_datetime(breaks = "2 hour", minor_breaks = "1 hour", date_labels = "%H:%M") +
#     coord_cartesian(xlim = (c(as.POSIXct("00:00:00", format="%H:%M:%S"),as.POSIXct("23:59:59", format="%H:%M:%S")))) +
#     labs(title = paste(levels(daily_full$Zone)[n]," zone - Summer",sep = ""), y = NULL, x = "Time") + 
#     scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
#     theme_minimal()
#   summer_graphs[paste(gsub(" ","_",levels(daily_full$Zone)[n])," Zone - Summer",sep = "")] <- list(graph)
#   rm(graph)
# }
# 
# #Winter
# winter_graphs <- list()
# for (n in 1:length(levels(daily_full$Zone))) {
#   graph <- ggplot(daily_full[daily_full[,"Zone"]==levels(daily_full$Zone)[n] & daily_full$Season=="Winter",], aes(Time)) + 
#     geom_point(aes(y = RH), color = "light blue", alpha = 0.25, size = .4, position = "jitter") + 
#     geom_smooth(aes(y = RH), color = "blue",) + 
#     geom_point(aes(y = Temp), color = "pink", alpha = 0.25,size = .4, position = "jitter") + 
#     geom_smooth(aes(y = Temp), color = "red",) +
#     scale_x_datetime(breaks = "2 hour", minor_breaks = "1 hour", date_labels = "%H:%M") +
#     coord_cartesian(xlim = (c(as.POSIXct("00:00:00", format="%H:%M:%S"),as.POSIXct("23:59:59", format="%H:%M:%S")))) +
#     labs(title = paste(levels(daily_full$Zone)[n]," zone - Winter",sep = ""), y = NULL, x = "Time") + 
#     scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
#     theme_minimal()
#   winter_graphs[paste("Daily_",gsub(" ","_",levels(daily_full$Zone)[n]),"_Zone_Winter",sep = "")] <- list(graph)
#   rm(graph)
# }
# 
# #Save images
# daily_image <- grid.arrange(daily_graphs[[1]],daily_graphs[[2]],daily_graphs[[3]],nrow = 3)
# summer_images <- grid.arrange(summer_graphs[[1]],summer_graphs[[2]],summer_graphs[[3]],nrow = 3)
# winter_images <- grid.arrange(winter_graphs[[1]],winter_graphs[[2]],winter_graphs[[3]],nrow = 3)
# all_images <- grid.arrange(daily_graphs[[1]] + labs(tag = "A"),daily_graphs[[2]]+ labs(tag = "B"),daily_graphs[[3]]+ labs(tag = "C"),summer_graphs[[1]]+ labs(tag = "D"),summer_graphs[[2]]+ labs(tag = "E"),summer_graphs[[3]]+ labs(tag = "F"),winter_graphs[[1]]+ labs(tag = "G"),winter_graphs[[2]]+ labs(tag = "H"),winter_graphs[[3]]+ labs(tag = "I"),nrow = 3)
# #SVG files are huge here because of all the individual points that are saved!
# #ggsave(filename = "Daily_Graphs.svg", all_images, width = 16, height = 9, dpi = 600)
# ggsave(filename = "Daily_Graphs.png", all_images, width = 16, height = 9, dpi = 600)
# ggsave(filename = "Daily_Graphs.jpg", all_images, width = 16, height = 9, dpi = 100)

##########

# #Comparison of iButton data
# iButton_RH <- ggplot(NASA_data, aes(YYYYMMDD)) +
#   stat_summary(aes(y = iRH, group = Zone, col = Zone), fun=mean, geom="line") +
#   stat_summary(aes(y = iRH, group = Zone, col = Zone), fun.data =mean_sdl,fun.args = list(mult = 1), geom="errorbar", alpha = 0.3) +
#   scale_x_datetime(breaks = "1 month", date_labels = "%b %Y") + 
#   labs(title = "iButton Relative Humidity", y = NULL, x = "Time") + 
#   scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
#   theme(axis.text.x=element_text(angle=45, hjust=1)) +
#   theme_minimal() +
#   scale_color_manual(values = Zone_cols)
# 
# iButton_Temp <- ggplot(NASA_data, aes(YYYYMMDD)) +
#   stat_summary(aes(y = iTemp, group = Zone, col = Zone), fun=mean, geom="line") +
#   stat_summary(aes(y = iTemp, group = Zone, col = Zone), fun.data =mean_sdl,fun.args = list(mult = 1), geom="errorbar", alpha = 0.3) +
#   scale_x_datetime(breaks = "1 month", date_labels = "%b %Y") + 
#   labs(title = "iButton Temperature", y = "Temperature (°C)", x = "Time") + 
#   scale_y_continuous(limits = c(0,45), breaks = c(0,5,10,15,20,25,30,35,40,45)) +
#   theme(axis.text.x=element_text(angle=45, hjust=1)) +
#   theme_minimal() +
#   scale_color_manual(values = Zone_cols)
# 
# ggsave(filename = paste("iButton_RH.svg",sep=""), iButton_RH, width = 16, height = 9, dpi = 600)
# ggsave(filename = paste("iButton_Temp.svg",sep=""), iButton_Temp, width = 16, height = 9, dpi = 600)
# ggsave(filename = paste("iButton_RH.jpg",sep=""), iButton_RH, width = 16, height = 9, dpi = 100)
# ggsave(filename = paste("iButton_Temp.jpg",sep=""), iButton_Temp, width = 16, height = 9, dpi = 100)
# 
# #Compare iButton RH with NASA RH2M and PRECTOT and weather station
# comp_RH <- list()
# for (n in levels(Namib_Weather_all_data$Zone)) {
#   graph <- ggplot(Namib_Weather_all_data[Namib_Weather_all_data$Zone==n,], aes(YYYYMMDD)) + 
#     stat_summary(aes(y = iRH, col="iButton soil RH (%)"), fun=mean, geom="line") + 
#     stat_summary(aes(y = RH2M, col="NASA air RH (%)"), fun=mean,fun.args = list(na.rm = T), geom="line") + 
#     stat_summary(aes(y = Humidity, col="SASSCAL air RH (%)"), fun=mean,fun.args = list(na.rm = T), geom="line") + 
#     stat_summary(aes(y = Precipitation, col="SASSCAL rain (mm)"), fun=mean,fun.args = list(na.rm = T), geom="line",alpha = 0.5) + 
#     stat_summary(aes(y = PRECTOT, col="NASA rain (mm)"), fun=mean, geom="line",alpha = 0.5) + 
#     scale_x_datetime(breaks = "1 month", date_labels = "%b %y") + 
#     labs(title = paste("Water availability - ",n,sep = ""), y = NULL, x = "Date") + 
#     scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
#     theme_minimal() +
#     scale_color_manual(name = "Group", values = Variable_cols)
#   comp_RH[n] <- list(graph)
#   rm(graph)
# }
# 
# #Compare iButton Temp with NASA T2M, T2M_Max and T2M_Min
# comp_Temp <- list()
# for (n in levels(Namib_Weather_all_data$Zone)) {
#   graph <- ggplot(Namib_Weather_all_data[Namib_Weather_all_data$Zone==n,], aes(YYYYMMDD)) + 
#     stat_summary(aes(y = iTemp, col="iButton soil temp."), fun=mean,fun.args = list(na.rm = T), geom="line") + 
#     stat_summary(aes(y = T2M, col="NASA air temp."), fun=mean,fun.args = list(na.rm = T), geom="line") +
#     stat_summary(aes(y = Air.temp.avg., col="SASSCAL air temp."), fun=mean,fun.args = list(na.rm = T), geom="line") +
#     stat_summary(aes(y = Soil.Temp.avg., col="SASSCAL soil temp."), fun=mean,fun.args = list(na.rm = T), geom="line") +
#     scale_x_datetime(breaks = "1 month", date_labels = "%b %y") + 
#     labs(title = paste("Temperature - ",n,sep = ""), y = "Temperature (°C)", x = "Date") + 
#     scale_y_continuous(limits = c(0,45), breaks = c(0,5,10,15,20,25,30,35,40,45)) +
#     theme_minimal() +
#     scale_color_manual(name = "Group", values = Variable_cols)
#   comp_Temp[n] <- list(graph)
#   rm(graph)
# }
# 
# #Nice comparison figure
# NASA_iButton_comparison <- grid.arrange(comp_RH[[1]] + labs(tag = "A"),comp_RH[[2]] + labs(tag = "B"),comp_RH[[3]] + labs(tag = "C"),comp_Temp[[1]] + labs(tag = "D"),comp_Temp[[2]] + labs(tag = "E"),comp_Temp[[3]] + labs(tag = "F"),ncol = 3)
# ggsave(filename = "NASA_iButton_comparison.svg", NASA_iButton_comparison, width = 21, height = 9, dpi = 600)
# ggsave(filename = "NASA_iButton_comparison.jpg", NASA_iButton_comparison, width = 21, height = 9, dpi = 100)
# 
# #Correlation table
# corr_table <- as.data.frame(matrix(nrow = 9, ncol = 5,dimnames = list(NULL,c("Comparison","Overall","Fog","Low rain","High rain"))))
# corr_table[,"Comparison"] <- c("iRH vs RH2M","iRH vs Humidity","RH2M vs Humidity","iTemp vs T2M","iTemp vs Air.temp.avg.","iTemp vs Soil.Temp.avg.","T2M vs Air.temp.avg.","iRH vs Precipitation","iRH vs PRECTOT")
# #overall
# corr_table[1,"Overall"] <- summary(lm(Namib_Weather_all_data$iRH~Namib_Weather_all_data$RH2M))$r.squared
# corr_table[2,"Overall"] <- summary(lm(Namib_Weather_all_data$iRH~Namib_Weather_all_data$Humidity))$r.squared
# corr_table[3,"Overall"] <- summary(lm(Namib_Weather_all_data$RH2M~Namib_Weather_all_data$Humidity))$r.squared
# corr_table[4,"Overall"] <- summary(lm(Namib_Weather_all_data$iTemp~Namib_Weather_all_data$T2M))$r.squared
# corr_table[5,"Overall"] <- summary(lm(Namib_Weather_all_data$iTemp~Namib_Weather_all_data$Air.temp.avg.))$r.squared
# corr_table[6,"Overall"] <- summary(lm(Namib_Weather_all_data$iTemp~Namib_Weather_all_data$Soil.Temp.avg.))$r.squared
# corr_table[7,"Overall"] <- summary(lm(Namib_Weather_all_data$T2M~Namib_Weather_all_data$Air.temp.avg.))$r.squared
# corr_table[8,"Overall"] <- summary(lm(Namib_Weather_all_data$iRH~Namib_Weather_all_data$Precipitation))$r.squared
# corr_table[9,"Overall"] <- summary(lm(Namib_Weather_all_data$iRH~Namib_Weather_all_data$PRECTOT))$r.squared
# 
# #zones
# for (x in unique(Namib_Weather_all_data$Zone)) {
#   corr_table[1,x] <- try(summary(lm(Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"iRH"]~Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"RH2M"]))$r.squared)
#   corr_table[2,x] <- try(summary(lm(Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"iRH"]~Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"Humidity"]))$r.squared)
#   corr_table[3,x] <- try(summary(lm(Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"RH2M"]~Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"Humidity"]))$r.squared)
#   corr_table[4,x] <- try(summary(lm(Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"iTemp"]~Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"T2M"]))$r.squared)
#   corr_table[5,x] <- try(summary(lm(Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"iTemp"]~Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"Air.temp.avg."]))$r.squared)
#   corr_table[6,x] <- try(summary(lm(Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"iTemp"]~Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"Soil.Temp.avg."]))$r.squared)
#   corr_table[7,x] <- try(summary(lm(Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"T2M"]~Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"Air.temp.avg."]))$r.squared)
#   corr_table[8,x] <- try(summary(lm(Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"iRH"]~Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"Precipitation"]))$r.squared)
#   corr_table[9,x] <- try(summary(lm(Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"iRH"]~Namib_Weather_all_data[Namib_Weather_all_data$Zone==x,"PRECTOT"]))$r.squared)
# }
# 
# write.table(corr_table,"measure_correlations.csv",row.names = F,col.names = T,quote = T,sep = ",")
# write.table(corr_table,"TableS5.csv",row.names = F,col.names = T,quote = T,sep = ",")


##########

# #Box plot of the iButton data
# Box_RH <- ggplot(iButton_complete, aes(x = Zone, y = RH)) + 
#   geom_boxplot(aes(fill=Zone)) + 
#   scale_fill_manual(values=Zone_cols) + 
#   theme_minimal() +
#   scale_y_continuous(n.breaks = 7) +
#   ylab(label = "Relative humidity (%)") +
#   theme(legend.position = "none") +
#   geom_signif(comparisons = list(c("Fog","Low rain"),c("Low rain","High rain"),c("Fog","High rain")), y_position = c(max(iButton_complete$RH,na.rm = TRUE)*1.05,(max(iButton_complete$RH,na.rm = TRUE))*1.1,(max(iButton_complete$RH,na.rm = TRUE))*1.15), map_signif_level=TRUE, colour="black",na.rm = TRUE)
# 
# Box_Temp <- ggplot(iButton_complete, aes(x = Zone, y = Temp)) + 
#   geom_boxplot(aes(fill=Zone)) + 
#   scale_fill_manual(values=Zone_cols) + 
#   theme_minimal() +
#   scale_y_continuous(n.breaks = 7) +
#   ylab(label = "Temperature (°C)") +
#   theme(legend.position = "none") +
#   geom_signif(comparisons = list(c("Fog","Low rain"),c("Low rain","High rain"),c("Fog","High rain")), y_position = c(max(iButton_complete$Temp,na.rm = TRUE)*1.05,(max(iButton_complete$Temp,na.rm = TRUE))*1.1,(max(iButton_complete$Temp,na.rm = TRUE))*1.15), map_signif_level=TRUE, colour="black",na.rm = TRUE)
# 
# supp_boxes <- grid.arrange(Box_RH + labs(tag = "A"),Box_Temp + labs(tag = "B"),nrow = 2)
# ggsave(filename = "FigS2.svg", supp_boxes, width = 10, height = 10, dpi = 600)

##########

# #Weather table
# weather_tables <- list()
# VARIABLES <- c("RH","Temp")
# for (VAR in VARIABLES) {
#   Weathertable <- data.frame()
#   #Transect
#   Weathertable[1,1:7] <- c(summary(iButton_complete[,VAR])[1:6], sd(iButton_complete[,VAR],na.rm = T))
#   Weathertable[1,8:14] <- c(summary(iButton_complete[iButton_complete$Season == "Summer",VAR])[1:6], sd(iButton_complete[iButton_complete$Season == "Summer",VAR],na.rm = T))
#   Weathertable[1,15:21] <- c(summary(iButton_complete[iButton_complete$Season == "Winter",VAR])[1:6], sd(iButton_complete[iButton_complete$Season == "Winter",VAR],na.rm = T))
#   #Fog / Low rain / High rain
#   for (n in 1:length(levels(iButton_complete$Zone))) {
#     Weathertable[1+n,1:7] <- c(summary(iButton_complete[iButton_complete$Zone==levels(iButton_complete$Zone)[n],VAR])[1:6], sd(iButton_complete[iButton_complete$Zone==levels(iButton_complete$Zone)[n],VAR],na.rm = T))
#     Weathertable[1+n,8:14] <- c(summary(iButton_complete[iButton_complete$Zone==levels(iButton_complete$Zone)[n] & iButton_complete$Season == "Summer",VAR])[1:6], sd(iButton_complete[iButton_complete$Zone==levels(iButton_complete$Zone)[n] & iButton_complete$Season == "Summer",VAR],na.rm = T))
#     Weathertable[1+n,15:21] <- c(summary(iButton_complete[iButton_complete$Zone==levels(iButton_complete$Zone)[n] & iButton_complete$Season == "Winter",VAR])[1:6], sd(iButton_complete[iButton_complete$Zone==levels(iButton_complete$Zone)[n] & iButton_complete$Season == "Winter",VAR],na.rm = T))
#   }
#   weather_tables[paste(VAR,sep = "")] <- list(Weathertable)
#   rm(Weathertable)
# }
# #Set up the main table columns
# Weathertable <- data.frame()
# Weathertable[1,1:23] <- c("","",c(rep("Annual",7),rep("Summer",7),rep("Winter",7)))
# Weathertable[2,1:23] <- c("","",rep(c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","sd"),3))
# #Bring in all the previously generated tables
# for (n in 1:length(VARIABLES)) {
#   Weathertable[((3*n)+(n-1)),1:23] <- c(VARIABLES[n],"Transect",weather_tables[[VARIABLES[n]]][1,])
#   Weathertable[(3*n+(n-1)+1),1:23] <- c(VARIABLES[n],"Fog",weather_tables[[VARIABLES[n]]][2,])
#   Weathertable[(3*n+(n-1)+2),1:23] <- c(VARIABLES[n],"Low Rain",weather_tables[[VARIABLES[n]]][3,])
#   Weathertable[(3*n+(n-1)+3),1:23] <- c(VARIABLES[n],"High Rain",weather_tables[[VARIABLES[n]]][4,])
# }
# #Save the table
# write.table(Weathertable,file="Weather_Table.csv",row.names = F,col.names = F,quote = T,sep = ",")

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
rain_effect_day$Zone <- factor(rain_effect_day$Zone, levels = c("Fog","Low rain","High rain"), ordered = TRUE)

#Remove the days when some zones have missing data.
rain_effect_day <- na.omit(rain_effect_day)

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

rain_figs <- list()
rain_table <- data.frame()
for (n in unique(rain_effect_day$Zone)) {
  #rain peak position
  peak_pos <- find_peaks(rain_effect_day[rain_effect_day$Zone==n,"Precipitation"],m=4)
  
  #rain start positions
  start_pos <- vector()
  i <- 1
  for (x in peak_pos) {
    step <- 1
    done <- FALSE
    while (done == FALSE) {
      if (rain_effect_day[rain_effect_day$Zone==n,"Precipitation"][x - step] > 0) {
        step <- step + 1
      }
      else {
        done <- TRUE
      }
    }
    start_pos[i] <- x - step + 1
    i <- i+1
  }
  
  #rain end positions
  end_pos <- vector()
  i <- 1
  for (x in peak_pos) {
    step <- 1
    done <- FALSE
    while (done == FALSE) {
      if (rain_effect_day[rain_effect_day$Zone==n,"Precipitation"][x + step] > 0) {
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
    precip_tot[i] <- sum(rain_effect_day[rain_effect_day$Zone==n,"Precipitation"][x:end_pos[i]])
    i <- i + 1
  }
  
  #Mean iRH the preceding week
  mean_iRH <- vector()
  i <- 1
  for (x in start_pos) {
    if (x-7 > 0) {
      mean_iRH[i] <- mean(rain_effect_day[rain_effect_day$Zone==n,"iRH"][(x-7):x],na.rm = TRUE)
    }
    else {
      mean_iRH[i] <- mean(rain_effect_day[rain_effect_day$Zone==n,"iRH"][0:x],na.rm = TRUE)
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
      if ((x + step + 1) > length(rain_effect_day[rain_effect_day$Zone==n,"iRH"])) {
        done <- TRUE
      }
      if (rain_effect_day[rain_effect_day$Zone==n,"iRH"][x + step] > mean_iRH[i]) {
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
    max_iRH[i] <- max(rain_effect_day[rain_effect_day$Zone==n,"iRH"][start_pos[i]:eff_end_pos[i]])
    i <- i + 1
  }
  
  #Construct a table
  current_rain_table <- data.frame(matrix(nrow = length(start_pos),dimnames = list(NULL,"Zone")))
  current_rain_table[,"Zone"] <- rep(n,length(start_pos))
  current_rain_table[,"Rain_start"] <- rain_effect_day[rain_effect_day$Zone==n,"YYYYMMDD"][start_pos]
  current_rain_table[,"Rain_end"] <- rain_effect_day[rain_effect_day$Zone==n,"YYYYMMDD"][end_pos]
  current_rain_table[,"Rain_duration"] <- end_pos - start_pos
  current_rain_table[,"Total_precipitation"] <- precip_tot
  current_rain_table[,"Effect_end"] <- rain_effect_day[rain_effect_day$Zone==n,"YYYYMMDD"][eff_end_pos]
  current_rain_table[,"Preceeding_iRH"] <- mean_iRH
  current_rain_table[,"Rain_max_iRH"] <- max_iRH
  current_rain_table[,"Rain_effect_duration"] <- eff_end_pos - start_pos
  current_rain_table[,"Rain_effect_magnitude"] <- max_iRH - mean_iRH
  
  #Add drying gradient
  for (z in 1:length(start_pos)) {
    valx <- (start_pos[z]:eff_end_pos[z])
    valy <- (rain_effect_day[rain_effect_day$Zone==n,"iRH"][start_pos[z]:eff_end_pos[z]])
    current_rain_table[z,"Drying_rate"] <- coef(lm(valy~valx))[2]
  }
  
  #Merge the tables
  rain_table <- rbind(rain_table,current_rain_table)
  
  #Plot the rain events
  current_rain_plot <- ggplot(rain_effect_day[rain_effect_day$Zone==n,], aes(YYYYMMDD)) + 
    geom_line(aes(y = iRH), colour = "blue") + 
    geom_line(aes(y = Precipitation)) + 
    scale_x_datetime(breaks = "1 month", date_labels = "%b %y") + 
    labs(title = paste("Rain events - ",n,sep = ""), y = NULL, x = "Date") + 
    scale_y_continuous(limits = c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
    geom_vline(xintercept = rain_effect_day[rain_effect_day$Zone==n,"YYYYMMDD"][peak_pos],colour="blue", alpha = 0.6) +
    geom_vline(xintercept = rain_effect_day[rain_effect_day$Zone==n,"YYYYMMDD"][eff_end_pos],colour="red", alpha = 0.6) +
    theme_minimal() +
    theme(legend.position = "none")
  
  rain_figs[n] <- list(current_rain_plot)
}

#Order zones
rain_table$Zone <- factor(rain_table$Zone, levels = c("Fog","Low rain","High rain"), ordered = TRUE)

#Rain effect duration vs precipitation
precip_dura <- ggplot(rain_table, aes(x = Total_precipitation,y = Rain_effect_duration)) +
  geom_point(aes(col=Zone)) +
  geom_smooth(method = "lm", alpha = 0.3, colour = "black") +
  scale_x_continuous(limits = c(0,18), breaks = c(0,5,9,14,18)) +
  scale_y_continuous(limits = c(0,24), breaks = c(0,6,12,18,24)) +
  labs(title = "Total precipitation per event versus the duration of effect", y = "Effect duration (days)", x = "Total precipitation (mm)") + 
  scale_colour_manual(values = Zone_cols) +
  theme_minimal()

#Drying rate graph
drying_rate_graph <- ggplot(rain_table, aes(y = Drying_rate, x = Zone)) +
  geom_boxplot(aes(col=Zone)) +
  labs(title = "Drying rate in the different zones", y = "Drying rate (%RH/day)") + 
  scale_colour_manual(values = Zone_cols) +
  geom_signif(comparisons = list(c("Fog","Low rain"),c("High rain","Low rain"),c("Fog","High rain")), y_position = c(1,1.5,2), map_signif_level=TRUE, colour="black") +
  theme_minimal()


#Save our rain data
rain_fig <- grid.arrange(rain_figs[[1]] + labs(tag = "A"),rain_figs[[2]] + labs(tag = "B"),rain_figs[[3]] + labs(tag = "C"),precip_dura + labs(tag = "D"),ncol = 2)
ggsave(filename = "Rain_graph.jpg", rain_fig, width = 16, height = 9, dpi = 200)
ggsave(filename = "Rain_graph.svg", rain_fig, width = 16, height = 9, dpi = 600)

ggsave(filename = "Drying_rate.jpg", drying_rate_graph, width = 16, height = 9, dpi = 200)
ggsave(filename = "Drying_rate.svg", drying_rate_graph, width = 16, height = 9, dpi = 600)

write.table(rain_table,"rain_table.csv",row.names = F,col.names = T,quote = T,sep = ",")

##########

#Water availability and activity table

water_activity_table <- as.data.frame(matrix(nrow = 5, ncol = 8, dimnames = list(NULL,c("Process","Water Activity Limit","Fog Zone1","Fog Zone2","Low rainfall zone1","Low rainfall zone2","High rainfall zone1","High rainfall zone2"))))
water_activity_table[,"Process"] <- c("Cell division","Photosynthesis","Extracellular protease","Amino acid metabolism","Basal metabolism (Maintenance energy)")
water_activity_table[,"Water Activity Limit"] <- c("0.9","0.8","0.4","0.3","0.2")
total_time <- as.numeric(difftime(max(iButton_complete[,"Date"]),min(iButton_complete[,"Date"]),units = "hours"))
#We get how many readings are above the threshold then divide by 3 because there are three sites per zone with data.
#Each reading covers a period of 4 hours so then we multiply by 4 to get the total hours.
#Fog
water_activity_table[1,"Fog Zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "Fog" & iButton_complete$RH >= 90,]))/3*4)
water_activity_table[2,"Fog Zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "Fog" & iButton_complete$RH >= 80,]))/3*4)
water_activity_table[3,"Fog Zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "Fog" & iButton_complete$RH >= 40,]))/3*4)
water_activity_table[4,"Fog Zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "Fog" & iButton_complete$RH >= 30,]))/3*4)
water_activity_table[5,"Fog Zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "Fog" & iButton_complete$RH >= 20,]))/3*4)
water_activity_table[1,"Fog Zone2"] <- round(water_activity_table[1,"Fog Zone1"]/total_time*100)
water_activity_table[2,"Fog Zone2"] <- round(water_activity_table[2,"Fog Zone1"]/total_time*100)
water_activity_table[3,"Fog Zone2"] <- round(water_activity_table[3,"Fog Zone1"]/total_time*100)
water_activity_table[4,"Fog Zone2"] <- round(water_activity_table[4,"Fog Zone1"]/total_time*100)
water_activity_table[5,"Fog Zone2"] <- round(water_activity_table[5,"Fog Zone1"]/total_time*100)
#Low rain
water_activity_table[1,"Low rainfall zone1"] <- round((nrow(na.omit(iButton_complete[iButton_complete$Zone == "Low rain" & iButton_complete$RH >= 90,])))/3*4)
water_activity_table[2,"Low rainfall zone1"] <- round((nrow(na.omit(iButton_complete[iButton_complete$Zone == "Low rain" & iButton_complete$RH >= 80,])))/3*4)
water_activity_table[3,"Low rainfall zone1"] <- round((nrow(na.omit(iButton_complete[iButton_complete$Zone == "Low rain" & iButton_complete$RH >= 40,])))/3*4)
water_activity_table[4,"Low rainfall zone1"] <- round((nrow(na.omit(iButton_complete[iButton_complete$Zone == "Low rain" & iButton_complete$RH >= 30,])))/3*4)
water_activity_table[5,"Low rainfall zone1"] <- round((nrow(na.omit(iButton_complete[iButton_complete$Zone == "Low rain" & iButton_complete$RH >= 20,])))/3*4)
water_activity_table[1,"Low rainfall zone2"] <- round(water_activity_table[1,"Low rainfall zone1"]/total_time*100)
water_activity_table[2,"Low rainfall zone2"] <- round(water_activity_table[2,"Low rainfall zone1"]/total_time*100)
water_activity_table[3,"Low rainfall zone2"] <- round(water_activity_table[3,"Low rainfall zone1"]/total_time*100)
water_activity_table[4,"Low rainfall zone2"] <- round(water_activity_table[4,"Low rainfall zone1"]/total_time*100)
water_activity_table[5,"Low rainfall zone2"] <- round(water_activity_table[5,"Low rainfall zone1"]/total_time*100)
#High rain
water_activity_table[1,"High rainfall zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "High rain" & iButton_complete$RH >= 90,]))/3*4)
water_activity_table[2,"High rainfall zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "High rain" & iButton_complete$RH >= 80,]))/3*4)
water_activity_table[3,"High rainfall zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "High rain" & iButton_complete$RH >= 40,]))/3*4)
water_activity_table[4,"High rainfall zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "High rain" & iButton_complete$RH >= 30,]))/3*4)
water_activity_table[5,"High rainfall zone1"] <- round((nrow(iButton_complete[iButton_complete$Zone == "High rain" & iButton_complete$RH >= 20,]))/3*4)
water_activity_table[1,"High rainfall zone2"] <- round(water_activity_table[1,"High rainfall zone1"]/total_time*100)
water_activity_table[2,"High rainfall zone2"] <- round(water_activity_table[2,"High rainfall zone1"]/total_time*100)
water_activity_table[3,"High rainfall zone2"] <- round(water_activity_table[3,"High rainfall zone1"]/total_time*100)
water_activity_table[4,"High rainfall zone2"] <- round(water_activity_table[4,"High rainfall zone1"]/total_time*100)
water_activity_table[5,"High rainfall zone2"] <- round(water_activity_table[5,"High rainfall zone1"]/total_time*100)
#export
write.table(water_activity_table,"water_activity_table.csv",row.names = F,col.names = T,quote = T,sep = ",")
#Table will be completed in Calc to fine tune the display.