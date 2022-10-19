#install + load appropriate packages, etc..
library(lubridate)
library(ggplot2)
library(nlme)
library(plotly)
library(dplyr)
library(tidyr)
library(magrittr)
library(cowplot)
library(emmeans)
library(ggpmisc)

###########################################################################
#2020 Farm + Urban Temp data from ARY for Overwintering temp guidelines

## Using Farm and Botanical Garden data for determining U vs R diff in winter temps
## First, wrangle data in separate tables, then combine into one set so you can visualize both U and R together

################################################################################
#read in data, check
setwd('C:/Users/prile/OneDrive/Documents/R/MartinLab_Projects/Acorn_Ant_OvW')

# Temperature over time - rural fall data set
WRtemp <- read.csv("Farm_16_MAR_2020 (1).csv")

#keep only temp and date cols
WRtemp <- WRtemp[,-c(3:4,6:7)]
head(WRtemp)
#change column names to better ones
colnames(WRtemp) <- c("Count", "Date", "RAir_Temp", "RSoil_Temp")

#split Date col into date and time cols with date including unique hour ID
WRtemp1 <- WRtemp %>%
  separate(col = Date,
           into = c("Date", "hour"), sep = " ")%>%
  filter(Count > 7594 & Count < 25738) #subset / filter to only include Dec 21 - Feb 21 temps

WRtemp1 <- WRtemp1 %>%
  separate(col = hour, into = c("hour", "m:s", sep = ":"))

#remove m:s and : columns
WRtemp1 <- WRtemp1[,-c(4:5)]

#change hour into number and Date into Date
WRtemp1$hour<- as.numeric(WRtemp1$hour)
class(WRtemp1$hour) #check
WRtemp1$Date <- mdy(WRtemp1$Date)


#check
WRtemp1[1:25,1:6]

#NEW goal: find mean of air and soil temp within unique hours 
WRtemp1 <- WRtemp1%>%
  group_by(Date, hour)%>%
  mutate(mean.air = mean(RAir_Temp))%>%
  mutate(mean.soil = mean(RSoil_Temp))

#check
WRtemp1[1:25,1:7]

### Part 3: Basic visualization: temp vs. date for rural
plot(mean.air ~ hour, data = WRtemp1)
plot(mean.soil ~ hour, data = WRtemp1)

ggplot(WRtemp1, aes(x = hour, y = mean.air))+
  geom_point(stat = 'summary', fun.y = 'mean') +
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0, fun.args = list(mult = 1))


###############
## Urban temp data - Cle Botanical Garden site
WUtemp <- read.csv("CBG_16_MAR_2020 (1).csv")
head(WUtemp)

#Part 2: Data Wrangling / Rearrangement
## Goal: adjust so that only the three month period of Dec 21 - Feb 21 is included
WUtemp <- WUtemp[,-c(3:4,6,8)]


WUtemp1$Date <- mdy(WUtemp1$Date)
class(WUtemp1$Date)

#change column names to better ones
colnames(WUtemp) <- c("Count", "Date", "UAir_Temp", "USoil_Temp")

#subset to only include Dec 21 - Feb 21 temps
WUtemp1 <- WUtemp %>%
  separate(col = Date,
           into = c("Date", "hour"), sep = " ")%>%
  filter(Count > 7594 & Count < 25738)

#separate Date into ind cols
WUtemp1 <- WUtemp1 %>%
  separate(col = hour, into = c("hour", "m:s", sep = ":"))

#remove extra : col and m:s col
WUtemp1 <- WUtemp1[,-c(4:5)]

#change to R date value
#change hour into factor to divide into groups for means
WUtemp1$hour<- as.numeric(WUtemp1$hour)
class(WUtemp1$hour) #check


#check date col
WUtemp1[1:25, 1:5]
summary(WUtemp1) #urban mean air temp = 2.23, mean soil = 2.5

#find mean of air and soil temp by hour - other method only worked for air temp
WUtemp1 <- WUtemp1%>%
  group_by(Date, hour)%>%
  mutate(mean.air = mean(UAir_Temp))%>%
  mutate(mean.soil = mean(USoil_Temp))

#check
WUtemp1[1:25,1:7]

### Part 3: Basic visualization: temp vs. date for rural
plot(mean.air ~ hour, data = WUtemp1)
plot(mean.soil ~ hour, data = WUtemp1)

ggplot(WUtemp1, aes(x = hour, y = mean.air))+
  geom_point(stat = 'summary', fun.y = 'mean') +
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0, fun.args = list(mult = 1))


########
##Combine data sets and then plot together to compare 
WinTemp <- merge(WUtemp1, WRtemp1, by = c("Count"))
head(WinTemp)

#change col names (change one urban and rural at a time for sep air and soil visualizations)
colnames(WinTemp) <- c("Count", "DateU", "hourU","Urban.Air", "Urban.Soil", "UrbanA",
                       "Urban", "DateR","hourR", "Rural.Air","Rural.Soil","RuralA","Rural")

### Data Wrangling Part 2: Pivot Longer
#GOAL: create one combined column of Air Temp (from  both sources)
#with new col of source (rural / urban), then new combined col of Soil Temp

WinTemp1 <- WinTemp %>%
  pivot_longer(c("Urban","Rural"), 
               names_to = "Source.pop", values_to = c("Mean.SoilTemp"))
WinTemp2<- WinTemp %>%  
  pivot_longer(c("Urban", "Rural"), 
               names_to = "Source.pop", values_to = c("Mean.AirTemp"))

WinTemp1$DateU <- mdy(WinTemp1$DateU)

#PROBLEM: need to fix duplicates: every count has 4 values instead of 2
#1. either remove every 3rd+ 4th row after 2nd + 3rd rows (eg. [,-c(2:3,6:7,10:11... )
#2 OR fix pivot longer above to place all source pops in one column and avoid duplication

#rown, rown+1, row+4, row+5, 

#remove extra date, time, and count cols
WinTemp1<- WinTemp1[,-c(4:10)]
#remove for air temps
WinTemp2 <- WinTemp2[,-c(4:11)]


#Part 3: visualize over full winter time period
## GOAL: visualize temp over time + monthly day:night avgs
plot(Urban.Soil ~ DateU, type = "l", col = "red", ylab = "Soil Temperature (C)", xlab = "Date",
     main = "Winter Urban and Rural Soil Temperatures", data = WinTemp1)
points(Rural.Soil ~ DateU, type = "l", col = "blue", data = WinTemp1)


WinTemp2$Source.pop<- as.factor(WinTemp1$Source.pop)
levels(WinTemp1$Source.pop)

class(WinTemp1$Source.pop)


#GOAL: visualize diurnal average temps
#plot(mean Temp C +/- SD ~ hour, colour = source pop)

library(ggplot2)
ggplot(WinTemp1, aes(x = hourU, y = Mean.SoilTemp, color = Source.pop))+
  geom_point(size = 4, stat = 'summary', fun.y = 'mean') +
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0, fun.args = list(mult = 1))+
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  labs(y = "Soil Temperature C)", x = "Hour", color = "Source Population")+
  ggtitle('Diurnal Hourly Winter Temperature')+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10))

ggplot(WinTemp2, aes(x = hourU, y = Mean.AirTemp, colour = Source.pop))+
  geom_point(size = 4, stat = 'summary', fun.y = 'mean') +
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0, fun.args = list(mult = 1))+ #Soil Temp
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  labs(y = "Air Temperature (C)", x = "Hour", colour = "Source Population")+
  ggtitle('Diurnal Hourly Winter Temperature')+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10))
    
##GOAL visualize monthly day:night averages
#first split date into month col
#then, check monthly avgs
class(WinTemp1$DateU)
WinTemp1 <- WinTemp1 %>%
  separate(col = DateU,
           into = c("Year", "Month", "Day"), sep = "-" )
WinTemp1 <- WinTemp1 %>%
  separate(col = DateR, 
           into = c("Year", "Month", "Day"), sep = "-" )
head(WinTemp1)

#change month to factor, then change names
WinTemp1$Month <- as.factor(WinTemp1$Month)
levels(WinTemp1$Month) <- c("November", "December", "January", "February", "March")

#draw out plot for monthly avgs
ggplot(WinTemp1, aes(x = Month, y = Mean.SoilTemp, color = Source.pop))+
  geom_point(size = 4, stat = 'summary', fun.y = 'mean') +
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0.05, fun.args = list(mult = 1))+
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  labs(y = "Soil Temperature C)", x = "Month", color = "Source Population")+
  ggtitle('Diurnal Hourly Winter Temperature')+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10))


##GOAL create day avg and night avg from each one
#need  mutate to create conditional column

#create hour number
WinTemp1$hourU <- as.numeric(WinTemp1$hourU)
library(purrr)
WinTemp1 <- WinTemp1 %>%  #create time col to make new factor to designate day / night
  mutate(time = ifelse(hourU >=8 & hourU <= 16, "day", "night"))

#now create a mean col for day and mean col for night
WinTemp1$time <- as.factor(WinTemp1$time)
class(WinTemp1$time)
WinTemp1 <- WinTemp1 %>%
  group_by(Month, Source.pop, time)%>%
  mutate(mean.time = mean(Mean.SoilTemp))

#plot data based on Temp ~ month + source pop (color) + time
ggplot(WinTemp1, aes(x = Month, y = mean.time, shape = time, colour = Source.pop))+
  geom_point(size = 4, stat = 'summary', fun.y = 'mean')+
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0.05, fun.args = list(mult = 1))+
  geom_line(aes(group = interaction(time, Source.pop)))+ #this allows lines to only go to either shape and source pop
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  labs(y = "Mean Soil Temperature (C)", x = "Month", color = "Source Population")+
  ggtitle('Monthly Day:Night Winter Temperature')+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10))


  

my_colorsTEMP <-c("cadetblue", "darkorange")
names(my_colorsTEMP) <- levels(WinTemp2$Source.pop)

my_colorsTEMPs <-c("cadetblue", "darkorange")
names(my_colorsTEMPs) <- levels(WinTemp1$Soil.Source.pop)



##############  Temp determination for OvW Fall stepdowns ###############3
################
#compute Fall temp means, medians, and quantiles
rur_temp_fall0 <- Fall_rural_temp$Temperature
rur_temp_fall1 <- Fall_rural_temp$Temperature1
rur_temp_fall2 <- Fall_rural_temp$Temperature2

temp_date_fall0 <- Fall_rural_temp$date
temp_date_fall1 <- Fall_rural_temp$date1
temp_date_fall2 <- Fall_rural_temp$date2

rur_temp_mean <- mean(Fall_rural_temp$Temperature, na.rm = T)
rur_temp_med <- median(Fall_rural_temp$Temperature, na.rm = T)

rur_temp_mean1 <- mean(Fall_rural_temp$Temperature1, na.rm = T)
rur_temp_med1 <- median(Fall_rural_temp$Temperature1, na.rm = T)

rur_temp_mean2 <- mean(Fall_rural_temp$Temperature2, na.rm = T)
rur_temp_med2 <- median(Fall_rural_temp$Temperature2, na.rm = T)
                        
rur_temp_quant <- quantile(Fall_rural_temp$Temperature, na.rm = T, 0.9)

#new Fall data frame w/ temp data
plot.type <- Fall_rural_temp$Plot.type
temp.summary <- data.frame(rur_temp_mean=rur_temp_mean, rur_temp_med=rur_temp_med, rur_temp_quant=rur_temp_quant, plot.type=plot.type)

#compute Winter temp means, medians, and quantiles
rur_temp_wint <- Winter_rural_temp$Temperature
rur_date_wint <- Winter_rural_temp$date

rur_tempw_mean <- mean(Winter_rural_temp$Temperature, na.rm = T)
rur_tempw_med <- median(Winter_rural_temp$Temperature, na.rm = T)
rur_tempw_quant <- quantile(Winter_rural_temp$Temperature, na.rm = T, 0.9)

#new Winter data frame w/ temp data
plot.typew <- Winter_rural_temp$Plot.type
temp.summary <- data.frame(rur_tempw_mean=rur_tempw_mean, rur_tempw_med=rur_tempw_med, rur_tempw_quant=rur_tempw_quant, plot.typew="Farm")


#plots of ibutton data over time
par(mfrow= c(2,1))
#plots of avg Fall temp data
boxplot(rur_temp_mean ~ plot.type, temp.summary, y_lim = c(-1, 5),xlab = "Farm Fall (mean)", ylab = "Fall Temperature(C)")
boxplot(rur_temp_med ~ plot.type, temp.summary, y_lim = c(25, 29),xlab = "Farm Fall (median)", ylab = "Temperature(C)")
boxplot(rur_temp_quant ~ plot.type, temp.summary, y_lim = c(32, 36), xlab = "Farm Fall (quantile)", ylab = "Temperature(C)")

#plots of avg Winter temp data
boxplot(rur_tempw_mean ~ plot.typew, temp.summary, y_lim = c(-1, 5),xlab = "Farm Winter (mean)", ylab = "Fall Temperature(C)")
boxplot(rur_tempw_med ~ plot.typew, temp.summary, y_lim = c(25, 29),xlab = "Farm Winter (median)", ylab = "Temperature(C)")
boxplot(rur_tempw_quant ~ plot.typew, temp.summary, y_lim = c(32, 36), xlab = "Farm Winter (quantile)", ylab = "Temperature(C)")

#boxplots of total temp
boxplot(rur_temp_fall ~ plot.type, y_lim = c(-1,6), xlab = "Fall Farm", ylab = "Temperature(C)")
boxplot(rur_temp_wint ~ plot.typew, y_lim = c(-1,6), xlab = "Winter Farm", ylab = "Temperature(C)")


#Plot over time
plot(Fall_rural_temp$Temperature, type = "l", col = "red", ylab = "Temperature (C)", xlab = "Time",
     main = "Fall + Winter Plot Temperatures")
  points(Winter_rural_temp$Temperature, type = "l", col = "blue")


####
#Urban site (Cleveland Botanical Garden - CBG)
# Temperature over time - urban fall data set
  Fall_urban_temp <- read.csv(file.choose())
  Winter_urban_temp <- read.csv(file.choose())
  
#compute urban Fall temp means, medians, and quantiles
  urb_temp_fall <- Fall_urban_temp$Temperature
  urb_date_fall <- Fall_urban_temp$date
  
  urb_temp_mean <- mean(Fall_urban_temp$Temperature, na.rm = T)
  urb_temp_med <- median(Fall_urban_temp$Temperature, na.rm = T)
  urb_temp_quant <- quantile(Fall_urban_temp$Temperature, na.rm = T, 0.9)
  
  #new Fall data frame w/ temp data
  uplot.type <- Fall_urban_temp$plot.type
  uftemp.summary <- data.frame(urb_temp_mean=urb_temp_mean, urb_temp_med=urb_temp_med, urb_temp_quant=urb_temp_quant, uplot.type="Urban")
  
  #compute urban Winter temp means, medians, and quantiles
  urb_temp_wint <- Winter_urban_temp$Temperature
  urb_date_wint <- Winter_urban_temp$date
  
  urb_tempw_mean <- mean(Winter_urban_temp$Temperature, na.rm = T)
  urb_tempw_med <- median(Winter_urban_temp$Temperature, na.rm = T)
  urb_tempw_quant <- quantile(Winter_urban_temp$Temperature, na.rm = T, 0.9)
  
  #new Winter data frame w/ temp data
 uplot.typew <- Winter_rural_temp$Plot.type
  uwtemp.summary <- data.frame(urb_tempw_mean=urb_tempw_mean, urb_tempw_med=urb_tempw_med, urb_tempw_quant=urb_tempw_quant, uplot.typew="Urban")
  
  #plots of ibutton data over time
  par(mfrow= c(2,1))
  #plots of avg Fall temp data
  boxplot(urb_temp_mean ~ uplot.type, uftemp.summary, y_lim = c(-1, 5),xlab = "Urban Fall (mean)", ylab = "Fall Temperature(C)")
  boxplot(urb_temp_med ~ plot.type, uftemp.summary, y_lim = c(25, 29),xlab = "Urban Fall (median)", ylab = "Temperature(C)")
  boxplot(urb_temp_quant ~ uplot.type, uftemp.summary, y_lim = c(32, 36), xlab = "Urban Fall (quantile)", ylab = "Temperature(C)")
  
  #plots of avg Winter temp data
  boxplot(urb_tempw_mean ~ uplot.typew, uwtemp.summary, y_lim = c(-1, 5),xlab = "Urban Winter (mean)", ylab = "Temperature(C)")
  boxplot(urb_tempw_med ~ uplot.typew, uwtemp.summary, y_lim = c(25, 29),xlab = "Urban Winter (median)", ylab = "Temperature(C)")
  boxplot(urb_tempw_quant ~ uplot.typew, uwtemp.summary, y_lim = c(32, 36), xlab = "Urban Winter (quantile)", ylab = "Temperature(C)")
  
  #boxplots of total temp
  boxplot(urb_temp_fall ~ uplot.type, y_lim = c(-1,6), xlab = "Fall Urban", ylab = "Temperature(C)")
  boxplot(urb_temp_wint ~ uplot.typew, y_lim = c(-1,6), xlab = "Winter Urban", ylab = "Temperature(C)")
  
  
  
  
  #Plot over time
  plot(Fall_rural_temp$Temperature, type = "l", col = "red", ylab = "Temperature (C)", xlab = "Time",
       main = "Fall + Winter Plot Temperatures")
  points(Winter_rural_temp$Temperature, type = "l", col = "blue")
  
  
  

#data analysis of temp data
#Linear model of disease spread change and Temp
mod_temp <- lm(exp_temp1 ~ exp_temp$Date._Time_out, data = exp_temp)
mod
#analysis of variance, summarize f stats, p value, and estimates
anova(mod_temp)
summary(mod_temp)


######
##try to get trace data over the whole winter period

#read in data
WRtemp <- read.csv("Farm_16_MAR_2020 (1).csv")

#remove cols that are air temp or humidity data
WRtemp <- WRtemp[,-c(3:6)]
head(WRtemp)
#change column names to better ones
colnames(WRtemp) <- c("Count", "Date", "USoil_Temp", "RSoil_Temp")

#split Date col into date and time cols with date including unique hour ID
WRtemp1 <- WRtemp %>%
  separate(col = Date,
           into = c("Date", "hour"), sep = " ")
head(WRtemp1)
WRtemp1 <- WRtemp1 %>%
  separate(col = hour, into = c("hour", "m:s", sep = ":"))

#remove the extra column
WRtemp1 <- WRtemp1[,-c(5)]


  pivot_longer(cols = c("hour"), values_to = c("Mean.SoilTemp"))

  

