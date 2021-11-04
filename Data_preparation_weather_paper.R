#Namib Weather Data Preparation
#This script sets up the data files which were analysed. The results of this script are made available for download and this
#script solely serves as a reference for how the raw data files were joined together in preparation for the analysis.
#Written by: Jason Bosch
#Last update: 19/10/2021

###############

#Load the required libraries
library(nasapower)
library(stringr)

#Set working directory
setwd("~/PostDoc/02_Namib_10yr_project/iButton_data/analysis/")

###############

##iButton Preparation

#This uses the original files retrieved from the iButtons and combined them into a single data set which can be used
#for subsequent analyses.

#Bring in data
sites <- c("2","4","6","8","10","12","14","16","18","20")
#2019
for (x in sites) {
  try(import_rh <- read.csv(paste("../C14 Transect_iButtons_2019/C14 transect_site ",x,"_RH_2019.csv",sep = ""),stringsAsFactors = F,skip = 19))
  try(import_temp <- read.csv(paste("../C14 Transect_iButtons_2019/C14 transect_site ",x,"_temp_2019.csv",sep = ""),stringsAsFactors = F,skip = 19))
  try(import_file <- merge.data.frame(import_rh, import_temp, by = "Date.Time"))
  try(for (n in 1:nrow(import_file)) {
    import_file[n,6] <- as.POSIXct(strptime(import_file[n,1],format = "%Y/%m/%d %I:%M:%S %p"))
  })
  if (nchar(x) == 1) {
    try(assign(paste("iButton_2019_0",x,sep = ""),import_file))
  }
  else {
    try(assign(paste("iButton_2019_",x,sep = ""),import_file))
  }
  try(rm(import_file))
  try(rm(import_temp))
  try(rm(import_rh))
}
#Because the C14-14 iButton was missed and collected on a subsequent trip, the date format is different.
#This section will fix the date format and perform what should have been done on import.
for (n in 1:nrow(iButton_2019_14)) {
  iButton_2019_14[n,6] <- as.POSIXct(strptime(iButton_2019_14[n,1],format = "%d/%m/%y %H:%M:%S"))
}
#2020
for (x in sites) {
  try(import_rh <- read.csv(paste("../C14_Transect_iButtons_2020/TR-",x,"_RH_Apr_2021.csv",sep = ""),stringsAsFactors = F,skip = 19))
  try(import_temp <- read.csv(paste("../C14_Transect_iButtons_2020/TR-",x,"_T_Apr_2021.csv",sep = ""),stringsAsFactors = F,skip = 19))
  try(import_file <- merge.data.frame(import_rh, import_temp, by = "Date.Time"))
  try(for (n in 1:nrow(import_file)) {
    import_file[n,6] <- as.POSIXct(strptime(import_file[n,1],format = "%d/%m/%y %H:%M:%S"))
  })
  if (nchar(x) == 1) {
    try(assign(paste("iButton_2020_0",x,sep = ""),import_file))
  }
  else {
    try(assign(paste("iButton_2020_",x,sep = ""),import_file))
  }
  try(rm(import_file))
  try(rm(import_temp))
  try(rm(import_rh))
}
#Site 14's data is missing so we create a blank file to fill in the gaps
iButton_2020_14 <- iButton_2020_02
iButton_2020_14[,3] <- NaN
iButton_2020_14[,5] <- NaN

#Create a single data frame with all the data
sites <- c("02","04","06","08","10","12","14","16","18","20")
iButton_complete <- data.frame()
#2019
for (x in 1:length(sites)) {
  current <- get(paste("iButton_2019_",sites[x],sep = ""))
  current[,7] <- sites[x]
  current[,8] <- "2019"
  current[,9] <- paste("2019_",sites[x],sep="")
  iButton_complete <- rbind(iButton_complete,current)
}
#2020
for (x in 1:length(sites)) {
  current <- get(paste("iButton_2020_",sites[x],sep = ""))
  current[,7] <- sites[x]
  current[,8] <- "2020"
  current[,9] <- paste("2020_",sites[x],sep="")
  iButton_complete <- rbind(iButton_complete,current)
}

#Specify if it falls in the fog zone.
iButton_complete[iButton_complete[,7]%in%c("02","04","06"),10] <- "Fog"
iButton_complete[iButton_complete[,7]%in%c("08","10","12","14","16","18","20"),10] <- "Rain"

#Label all the columns
colnames(iButton_complete) <- c("date.time","unit.RH","RH","unit.Temp","Temp","Date","Site","Year","Group","Zone")

#Remove the RH values which are out of bounds and do not make sense
iButton_complete[which(iButton_complete$RH > 100),"RH"] <- 100
iButton_complete[which(iButton_complete$RH < 0),"RH"] <- NaN

#Write the data to a table
write.csv(iButton_complete,"iButton_Complete_dataset.csv",row.names = F)

###############

##NASA preparation

#This pulls in data from NASA Power satellite results for all sites and measurements of interest.

# PRECTOTCORR (Precipitation Corrected, mm)
# QV2M (Specific Humidity at 2 Meters, g/kg)
# RH2M (Relative Humidity at 2 Meters, %)
# T2M (Temperature at 2 Meters, C)
# T2M_MAX (Maximum Temperature at 2 Meters, C)
# T2M_MIN (Minimum Temperature at 2 Meters, C)

#Site co-ordinates
#Due to the resolution, many sites will unfortunately share the same data.
#Need to increase the first site's latitude by 0.01 as it is on the border of block.
coords <- as.data.frame(matrix(nrow = 10,ncol = 3,dimnames = list(NULL,c("Site","Lat","Lon"))))
coords$Site <- c("2","4","6","8","10","12","14","16","18","20")
coords$Lat <- c("-23.01","-23.02","-23.07","-23.14","-23.25","-23.31","-23.33", "-23.32","-23.35","-23.24")
coords$Lon <- c("14.67","14.86","15.04","15.21","15.36","15.53","15.71","15.86","16.01","16.14")

NASA_data <- data.frame()
#2019
for (x in 1:nrow(coords)) {
current_data <- get_power(community = "ag",
                          pars = c("PRECTOTCORR","QV2M","RH2M","T2M","T2M_MAX","T2M_MIN"),
                          temporal_api = "DAILY",
                          lonlat = c(as.numeric(coords[x,"Lon"]),as.numeric(coords[x,"Lat"])),
                          dates = c("2018-04-09","2019-03-17"))
  NASA_data <- rbind(NASA_data,current_data)
}
#2020
for (x in 1:nrow(coords)) {
current_data <- get_power(community = "ag",
                          pars = c("PRECTOTCORR","QV2M","RH2M","T2M","T2M_MAX","T2M_MIN"),
                          temporal_api = "DAILY",
                          lonlat = c(as.numeric(coords[x,"Lon"]),as.numeric(coords[x,"Lat"])),
                          dates = c("2019-05-01","2020-04-06"))
  NASA_data <- rbind(NASA_data,current_data)
}

write.csv(NASA_data,"NASA_data.csv",row.names = F)

#To remove any weird formatting of the data frame
NASA_data <- read.csv("NASA_data.csv")

#Since the coords got slightly off we need to round them first
NASA_data$LON <- signif(NASA_data$LON,4)
NASA_data$LAT <- signif(NASA_data$LAT,4)

#Need to make sure NASA data mentions the sites correctly
for (x in 1:nrow(coords)) {
  NASA_data[NASA_data$LAT == as.numeric(coords$Lat[x]), "Site"] <- coords$Site[x]
}

#Find everything in iButton data with the same data and average the iButton data for each day (slow)
for (x in 1:nrow(NASA_data)) {
  NASA_data[x,c("iRH","iTemp")] <- colMeans(iButton_complete[iButton_complete$Site == NASA_data[x,"Site"] & str_extract(iButton_complete$Date, "^.{10}") == NASA_data[x,"YYYYMMDD"],c("RH","Temp")])
}

#Add the zone markers
NASA_data[NASA_data[,"Site"]%in%c("2","4","6"),"Zone"] <- "Fog"
NASA_data[NASA_data[,"Site"]%in%c("8","10","12","14","16","18","20"),"Zone"] <- "Rain"

write.csv(NASA_data,"NASA_iButton_merged_data.csv",row.names = F)

###############

#SASSCAL preparation

#The SASSCAL data was originally downloaded from http://www.sasscalweathernet.org/

#Bring in the data files
#2019
CoastalMet_2019 <- read.table("../../Gobabeb_weather_data/SASSCAL WeatherNet Data for the 9 Stations Avail requested/SASSCAL_WeatherNet_Data April2018_2019 Daily values/2018-04-01_D_CoastalMet_E7631_TPAC.csv",sep=";",comment.char = "",header = T)
Kleinberg_2019 <- read.table("../../Gobabeb_weather_data/SASSCAL WeatherNet Data for the 9 Stations Avail requested/SASSCAL_WeatherNet_Data April2018_2019 Daily values/2018-04-01_D_Kleinberg-FN_E7630_TPAC.csv",sep=";",comment.char = "",header = T)
SophiesHoogte_2019 <- read.table("../../Gobabeb_weather_data/SASSCAL WeatherNet Data for the 9 Stations Avail requested/SASSCAL_WeatherNet_Data April2018_2019 Daily values/2018-04-01_D_SophiesHoogte_E7629_TPAC.csv",sep=";",comment.char = "",header = T)
Vogelfederberg_2019 <- read.table("../../Gobabeb_weather_data/SASSCAL WeatherNet Data for the 9 Stations Avail requested/SASSCAL_WeatherNet_Data April2018_2019 Daily values/2018-04-01_D_Vogelfederberg_E7626_TPAC.csv",sep=";",comment.char = "",header = T)
GarnetKoppie_2019 <- read.table("../../Gobabeb_weather_data/SASSCAL WeatherNet Data for the 9 Stations Avail requested/SASSCAL_WeatherNet_Data April2018_2019 Daily values/2018-04-01_D_GarnetKoppie_E7628_TPAC.csv",sep=";",comment.char = "",header = T)
Rooisand_2019 <- read.table("../../Gobabeb_weather_data/SASSCAL WeatherNet Data for the 9 Stations Avail requested/SASSCAL_WeatherNet_Data April2018_2019 Daily values/2018-04-01_D_Rooisand_103_ONRA.csv",sep=";",comment.char = "",header = T)
#2020
CoastalMet_2020 <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_CoastalMet_E7631_CRNB.csv",sep=";",comment.char = "",header = T)
Kleinberg_2020 <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_Kleinberg-FN_E7630_CRNB.csv",sep=";",comment.char = "",header = T)
SophiesHoogte_2020 <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_SophiesHoogte_E7629_CRNB.csv",sep=";",comment.char = "",header = T)
Vogelfederberg_2020 <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_Vogelfederberg_E7626_CRNB.csv",sep=";",comment.char = "",header = T)
GarnetKoppie_2020 <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_GarnetKoppie_E7628_CRNB.csv",sep=";",comment.char = "",header = T)
Rooisand_2020 <- read.table("../../Gobabeb_weather_data/2020_data/2019-05-01_D_Rooisand_103_MBBW.csv",sep=";",comment.char = "",header = T)

#Keep only columns of interest which are also present in all datasets
#2019
CoastalMet_2019 <- subset(CoastalMet_2019,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
Kleinberg_2019 <- subset(Kleinberg_2019,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
SophiesHoogte_2019 <- subset(SophiesHoogte_2019,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
Vogelfederberg_2019 <- subset(Vogelfederberg_2019,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
GarnetKoppie_2019 <- subset(GarnetKoppie_2019,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
Rooisand_2019 <- subset(Rooisand_2019,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
#2020
CoastalMet_2020 <- subset(CoastalMet_2020,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
Kleinberg_2020 <- subset(Kleinberg_2020,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
SophiesHoogte_2020 <- subset(SophiesHoogte_2020,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
Vogelfederberg_2020 <- subset(Vogelfederberg_2020,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
GarnetKoppie_2020 <- subset(GarnetKoppie_2020,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))
Rooisand_2020 <- subset(Rooisand_2020,select = c(Station.Name,Date,Air.temp...avg.,Air.temp...min.,Air.temp...max.,Precip...total.,Humidity,Soil.Temp.1..avg.))

#Remove the top line which has the units
#2019
CoastalMet_2019 <- CoastalMet_2019[-1,]
Kleinberg_2019 <- Kleinberg_2019[-1,]
SophiesHoogte_2019 <- SophiesHoogte_2019[-1,]
Vogelfederberg_2019 <- Vogelfederberg_2019[-1,]
GarnetKoppie_2019 <- GarnetKoppie_2019[-1,]
Rooisand_2019 <- Rooisand_2019[-1,]
#2020
CoastalMet_2020 <- CoastalMet_2020[-1,]
Kleinberg_2020 <- Kleinberg_2020[-1,]
SophiesHoogte_2020 <- SophiesHoogte_2020[-1,]
Vogelfederberg_2020 <- Vogelfederberg_2020[-1,]
GarnetKoppie_2020 <- GarnetKoppie_2020[-1,]
Rooisand_2020 <- Rooisand_2020[-1,]

#Rooisand_2020 has datapoints which are missing in all other 2020 data sets. This causes downstream problems if not removed.
Rooisand_2020 <- subset(Rooisand_2020,Rooisand_2020$Date%in%CoastalMet_2020$Date)

#Combine all the data
weather_station_complete <- Reduce(rbind,list(CoastalMet_2019,Kleinberg_2019,SophiesHoogte_2019,Vogelfederberg_2019,GarnetKoppie_2019,Rooisand_2019,
                                              CoastalMet_2020,Kleinberg_2020,SophiesHoogte_2020,Vogelfederberg_2020,GarnetKoppie_2020,Rooisand_2020))

#Some stations have missing values which need to be removed
for (x in 1:ncol(weather_station_complete)) {
  weather_station_complete[which(weather_station_complete[,x]==" - "),x] <- NA
}

#Replace with nicer column names
colnames(weather_station_complete) <- c("Station","Date","Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.")

#Fog or Rainfall zone?
weather_station_complete[weather_station_complete[,"Station"]%in%c("Coastal Met","Kleinberg-FN","Sophies Hoogte","Vogelfederberg"),"Zone"] <- "Fog"
weather_station_complete[weather_station_complete[,"Station"]%in%c("Garnet Koppie","Rooisand"),"Zone"] <- "Rain"

#Create a new POSIX compliant date column.
for (n in 1:nrow(weather_station_complete)) {
  weather_station_complete[n,"Date2"] <- as.POSIXct(weather_station_complete[n,"Date"],format = "%Y-%m-%d")
}
weather_station_complete["Date"] <- weather_station_complete["Date2"]
weather_station_complete <- subset(weather_station_complete, select = -Date2)

#Save the new file for later use
write.csv(weather_station_complete,"weather_station_complete.csv",row.names = F)

#Bring in data we need (also changes the column types from character)
weather_station_complete <- read.csv("weather_station_complete.csv",stringsAsFactors = F,colClasses = c(Date="POSIXct"))

#Each site will be assigned/averaged the data from the closest weather station
#Site 2 will get the average of coastal Met and Kleinberg
#Site 4 will get Sophies Hoogte
#Site 6 will get Vogelfederberg
#Sites 8, 10 and 12 will get Garnet Koppie
#Sites 14, 16, 18 and 20 will get Rooisand

#To merge it with the other data, we will need to summarise per site
weather_station_averaged <- data.frame()
for (x in 1:length(unique(weather_station_complete$Date))) {
  current_res <- data.frame()
  current_res[1:10,"Date"] <- unique(weather_station_complete$Date)[x]
  current_res[,"Site"] <- c("2","4","6","8","10","12","14","16","18","20")
  current_res[1,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.","Zone")] <-
    c(colMeans(weather_station_complete[weather_station_complete$Station%in%c("Coastal Met","Kleinberg-FN") & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.")],na.rm = T),"Fog")
  current_res[2,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.","Zone")] <-
    c(weather_station_complete[weather_station_complete$Station=="Sophies Hoogte" & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.")],"Fog")
  current_res[3,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.","Zone")] <-
    c(weather_station_complete[weather_station_complete$Station=="Vogelfederberg" & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.")],"Fog")
  current_res[4:6,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.","Zone")] <-
    c(weather_station_complete[weather_station_complete$Station=="Garnet Koppie" & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.")],"Rain")
  current_res[7:10,c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.","Zone")] <-
    c(weather_station_complete[weather_station_complete$Station=="Rooisand" & weather_station_complete$Date == unique(weather_station_complete$Date)[x],c("Air.temp.avg.","Air.temp.min.","Air.temp.max.","Precipitation","Humidity","Soil.Temp.avg.")],"Rain")
  weather_station_averaged <- rbind(weather_station_averaged,current_res)
}

#Write the weather station site data
write.csv(weather_station_averaged,"weather_station_by_site.csv",row.names = F)