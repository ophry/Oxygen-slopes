# Oxygen-slopes
install.packages("devtools")
library(devtools)
install_github('Echinophoria/container')
install.packages("reshape")

library(container)
read_pHlog
################################################################################################################
#### This functions structures the txt file from firesting. Dont runn the last part. need to modify to split channels. 

read_pHlog <- function (file) 
{ 
  file="data/2020-11-10_Sixth.txt"
  signe <- readLines(file, n = 44) ## definerar objekt, børja læsa på linje 45
  log_name <- grep("#Log Nam", signe, value = TRUE)# matchar loggnamn med Signe 
  log_name <- substr(log_name, (gregexpr("\t",log_name)[[1]] + #nchar= number of characters in log name. vart string ska børja och vart det ska sluta. førsta delen behandlar bara metadata. 
                                  1), nchar(log_name))
  st <- grep("#settings", signe, value = FALSE)# skapar et objekt, settings= en string med text och itexten anvænder dom #. säger att pH-log værde är falskt= returns the string. depends on if you want the value or the position 
  lnm <- unlist(strsplit(signe[st], split = "\t"))[-1]# lnm=log namn? splitter pHlog och separerar med space= "\t"
  settings <- unlist(strsplit(signe[st + 1], split = "\t"))[-1]# nytt objekt med settings också hær delas vektorn upp och separeras med space. men vad betyder - eller +1? 
  
  settings <- as.numeric(settings)# säger att objekt settings är numeriskt. 
  names(settings) <- lnm # lognmanes enligt settings? 
  settings <- list(settings)# listar upp settings
  cal <- grep("#calibration", signe, value = FALSE)# object calibration, matchar mot "pHlog"
  lnm <- unlist(strsplit(signe[cal], split = "\t"))[-1] # simplify list of vector, splitt with space. 
  calibration <- unlist(strsplit(signe[cal + 1], split = "\t"))[-1]
  calibration <- as.numeric(calibration)# plotting data to check data distrubution? 
  names(calibration) <- lnm
  calibration <- list(calibration)
  hd <- grep("Date", signe, value = FALSE)# header  
  nms <- grep("Date", signe, value = TRUE)# return numeric names in a data frame 
  hdr <- read.table(textConnection(nms), sep = c("\t"), # hdr= highest density regions 
                    header = TRUE)
  hdr <- names(hdr)
  vars <- c()### name of variables 
  units <- c()## unit of eac varibale 
  for (i in 1:length(hdr)) {
    h <- hdr[i]
    h <- gsub(pattern = "\\.{2}", replacement = "\t", 
              h)## matches and remove untis for matching strings of data?
    h <- strsplit(h, "\t")
    vars[i] <- h[[1]][1]
    units[i] <- h[[1]][2]
  }
  names(units) <- vars
  units <- list(units)
  o2.log <- read.table(file, skip = hd, sep = "\t", header = FALSE)### At the moment i stopp here. 
  hdr <- read.table(textConnection(nms), sep = c("\t"), # hdr= highest density regions 
                    header = TRUE)
  
  hdr <- names(hdr)
  vars <- c()### name of variables 
  units <- c()## unit of eac varibale 
  for (i in 1:length(hdr)) {
    h <- hdr[i]
    h <- gsub(pattern = "\\.{2}", replacement = "\t", 
              h)## matches and remove untis for matching strings of data?
    h <- strsplit(h, "\t")
    vars[i] <- h[[1]][1]
    units[i] <- h[[1]][2]
  }
  names(units) <- vars
  units <- list(units)
  o2.log <- read.table("data/2020-11-10_Sixth.txt", skip = hd, sep = "\t", header = FALSE)
  names(o2.log) <- vars
  o2.log$Date <- strptime(o2.log$Date, format = "%Y-%m-%d %H:%M:%S", 
                          tz = "")## extract to columns of the data set. Time and pH 
  o2.log$Date <- as.POSIXct(o2.log$Date)#### manipulate objects of classes 
  attr(o2.log$Date, "tzone") <- "UTC"
  colnames(o2.log)
  names(o2.log)[4] <- "cha.1"
  names(o2.log)[22] <- "cha.2"
  names(o2.log)[3] <- "dt.1"
  signe <- list(o2.log = o2.log, settings = settings, calibration = calibration, 
                variable_units = units)
  
  return(signe)#. 
}

######## Tar ut listan o2.log ur Signe.  
o2.log = signe[[1]]

#### Convert dataset from wide to long format. Slår samman alla kanal-værden i en kolumn och lægger till en ny ID-kolumn. 
data2 <- o2.log 
library(tidyr)
data_long <- gather(data2, channel, measurement, cha.1, cha.2) ### Før få unika kolumner så måste precisera exakt vad den ska anvænda. 


#### Convert to Dissolved Oxygen mg/L. Antonio comment: Then DO is the column in O2.log with the saturation and temperature is one of the columns with sample temperature in the O2.log, pressure is also one of the columns and RH = 100.
#kjør functionen så att du ser att den dyker upp i environment. Kalla på funksjonen linje 104. 
conv_stO2 <- function(DO, Temperature, Pressure, Salinity, RH){
  sal_factor <- exp(-Salinity * (0.017674-10.754/(Temperature+273.15)^2))
  DO_st_press <- exp(-139.34411+(1.575701*10^5)/(273.15+Temperature)-(6.642308*10^7)/(273.15+Temperature)^2+(1.2438*10^10)/(273.15+Temperature)^3-(8.621949*10^11)/(273.15+Temperature)^4)
  VP_RH <- (10^(8.107131-1730.63/(235+Temperature)))*0.13332239 * RH / 100
  Theta <- 0.000975-(1.426*10^-5)*Temperature+(6.436*10^-8)*Temperature^2
  Press_factor <- ((Pressure/101.325)-(VP_RH/101.325))*(1-Theta*(Pressure/101.325))/(1-(VP_RH/101.325))*(1-Theta)
  DO_mg_l <- Press_factor*DO_st_press*sal_factor*DO/100
  return(DO_mg_l)
}

data_long$DO_mg_l_measurement <- conv_stO2(data_long$measurement, data_long$Sample.Temp, data_long$Pressure, 35, 100)

####Extrahera data och gør en ny df som plottas i ggplot. OBS gør det før en kanal åt gången. 
data3 <- data_long 
data3 = cbind(data_long$dt.1,data_long$DO_mg_l_measurement)### R fattar inte att det ær en dfmåste deffinera det
names(data3) <- c('time', 'DO_mg_l_measurement') ## C betyder "combine" 
data3 <- as.data.frame(data3)
names(data3) <- c('time', 'DO_mg_l_measurement')
######ggplot to check the data 
library(ggplot2) 
p<- ggplot()
p<- p + geom_point(data=data3, aes(x=time, y=DO_mg_l_measurement, color=DO_mg_l_measurement)) 
p<- p + theme_gray(base_size=20)
p<- p + labs(x="time", y="Oxygen_consumption")
p

#######Subtractinc control channel frome polychaete-channels
data_long %>% dplyr::mutate(control.extracted = "cha.1-cha.2") 
