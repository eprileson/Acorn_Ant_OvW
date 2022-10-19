#######
############ Acorn Ant OvW - MR Batch Processing ############
###
setwd("C:/Users/prile/Box/Research/AcornAntOverwintering2021/MetabolicRate/Eric_OvW")

library(lubridate)
library(car)

#read in data
dat_ovw.3c <- read.csv("AA_labels_edits.csv")
head(dat_ovw.3c)

#adjust dates to R from csv for both date cols
dat_ovw.3c$Trial_date4<- mdy(dat_ovw.3c$Trial_date4)
dat_ovw.3c$Trial_date10<- mdy(dat_ovw.3c$Trial_date10)

#merge gilibrator readings into dataset, adjust dates to R from csv
flow.rates <- read.csv("Gilibrator readings_T.curv.csv")
head(flow.rates)
flow.rates <- flow.rates[!is.na(flow.rates$Incubator_Temp), ]
flow.rates$Date <- mdy(flow.rates$Date)

#make 2 objects to use as flow rates at the two test temps
flow.rates4 <- flow.rates[flow.rates$Incubator_Temp == 4, ]
flow.rates10 <- flow.rates[flow.rates$Incubator_Temp == 10, ]

#create two storage areas for temp data from gilibrator
dat_ovw.3c$Gillibrator4 <- c()
dat_ovw.3c$Gillibrator10 <- c()

#use for loop to merge Gillibrator into MR data; this creates two columns from the 
#gilibrator data set of the actual flow rates to use for the fracflow calcs
for(i in 1:dim(dat_ovw.3c)[1]) { #for every row in the data set,
  
  if(!is.na(dat_ovw.3c$Trial_date4[i])) { #if not NA from trial date in row, make new col w/ actual flow rate from Gil at that temp
    dat_ovw.3c$Gilibrator4[i]<-flow.rates4[which.min(abs(dat_ovw.3c$Trial_date4[i]-flow.rates4$Date)), "Animal_gilibrator_avg."]
  } else(dat_ovw.3c$Gilibrator4[i]<-NA)
  
  if(!is.na(dat_ovw.3c$Trial_date10[i])) { #repeat above for the 10 deg temp
    dat_ovw.3c$Gilibrator10[i]<-flow.rates10[which.min(abs(dat_ovw.3c$Trial_date10[i]-flow.rates10$Date)), "Animal_gilibrator_avg."]
  } else(dat_ovw.3c$Gilibrator10[i]<-NA)
}

#check that the two columns have the merged gillibrator data
head(dat_ovw.3c)

#install / load packages for reading exp files from Sable ExpeData
install.packages("remotes")
remotes::install_github("hawkmoth/sablebase")
install.packages("zoo")
library(zoo)
library(SableBase)

# list all files in the folder / from the working dir
files = dir(pattern="exp") #find all files in directory that have "exp"

files = files[files %in% c(dat_ovw.3c$MR4_name, dat_ovw.3c$MR10_file)==TRUE]

storage=c()


for (i in 1:length(files)) { #for every row in files / each file
  
  # read in files 
  sscf = read.sscf(files[i]) #use sablebase function to read files from exp data
  
  # assign to dataframe object
  sscf.file = data.frame(sscf, check.names=TRUE)
  #temp storage
  temp=c()
  
  j.indices <- unique(sscf.file$I_O_N) #make object w/ unique rows for each unique chamber from within each file
  
  for(j in 1:length(j.indices)) { #nested loop: for each of the unique file rows w/ unique chambers
    
    sscf.sub = sscf.file[sscf.file$I_O_N==j.indices[j],] #make new object for subset that only has unique chambers
    check.length = dim(sscf.sub)[1]
    
    if(check.length >=500) {
      
      # sscf.sub$fracflow = sscf.sub$fracCO2*500 # NEED TO REPLACE W/ ACTUAL FLOW RATE
      sscf.sub$fracCO2 = sscf.sub$CO2.ANIMAL/1000000

      flowrate.4 = match(dat_ovw.3c$MR4_name, files[i])
      flowrate.10 = match(dat_ovw.3c$MR10_name, files[i])
      
      if(sum(flowrate.4, na.rm=T) >=1 ) {
        flowrate = dat_ovw.3c$Gilibrator4[which(is.na(flowrate.4)==F)[1]]
      }
      
      if(sum(flowrate.10, na.rm=T) >=1 ) {
        flowrate = dat_ovw.3c$Gilibrator10[which(is.na(flowrate.10)==F)[1]]
      }
      
      sscf.sub$fracflow = sscf.sub$fracCO2*flowrate        
      
      # Then if you average the last 5 mins of each butterflies fracflow you should get mean CO2ML/min for each of your samples
      # This should be the last 300 readings since there are 10 minutes total in the trial?
      
      sscf.summary = mean(sscf.sub$fracflow[300:600])

      # Use rolling function to see where slope changes are the flattest over the entire trace 
      sscf.sub$index<-seq(1:dim(sscf.sub)[1])
      Coef <- function(Z) coef(lm(fracflow ~ index, as.data.frame(Z)))    
      Slopes <- rollapplyr(zoo(sscf.sub), 300, Coef, by.column = FALSE)
      start.index <- which.min(abs(Slopes$index))[1]
      sscf.summaryRA <- mean(sscf.sub$fracflow[start.index:(start.index+300)])
      
    } else {
      sscf.summary<-NA
      sscf.summaryRA<-NA
      start.index<-NA}
    
    temp = rbind(temp, cbind(sscf.summary, sscf.summaryRA, j.indices[j], files[i], check.length, start.index))
    
  }
  
  storage = rbind(storage, temp)
}

#new data frame w/ mean MR data, exp file
dat.sum<-data.frame(MeanMR=as.numeric(as.character(storage[,1])), MeanMR_RA=as.numeric(as.character(storage[,2])), ChamberIndex=as.numeric(as.character(storage[,3])), ExpeDataFile=storage[,4], DataLength=storage[,5], StartIndex=storage[,6])

head(dat.sum)

#make new csv file w/ mean MR data
setwd("C:/Users/prile/Box/Research/AcornAntOverwintering2021/MetabolicRate/Eric_OvW")
write.csv(dat.sum, "AA_MR_Summary_Stats.csv")


setwd("C:/Users/sarah/Box/DiamondLabResearch/Research/VCardui_Ontogeny_2021/SED_Analyses")

# write.csv(dat.sum, "VCB_MR_Summary_Stats.csv")
#make new column combines data file anme and chamber index #
dat.sum$merge.index<-paste(dat.sum$ExpeDataFile, dat.sum$ChamberIndex, sep="_")

dat_ovw.3c<-read.csv("AA_labels_edits.csv", header=T)
head(dat_ovw.3c)
dat_ovw.3c4 <-dat_ovw.3c[,c(1:9)]
dat_ovw.3c$MRTestTemp <-rep(20, times=dim(dat_ovw.3c4)[1])
colnames(dat_ovw.3c4) <-recode(colnames(dat_ovw.3c4), "'MR4_file'='MR_file';'MR4_name'='MR_name'")
dat_ovw.3c10 <-dat_ovw.3c[,c(1:9)]
dat_ovw.3c$MRTestTemp<-rep(30, times=dim(dat_ovw.3c10)[1])
colnames(dat_ovw.3c10)<-recode(colnames(dat_ovw.3c10), "'MR10_file'='MR_file';'MR10_name'='MR_name'")

#merge the two sets so that the 10 and 4 temp tests are merged in one file column
all.trait<-merge(dat_ovw.3c10, dat_ovw.3c4, all=T)
all.trait$merge.index<-paste(all.trait$MR_name, all.trait$MRchamber_num, sep="_")

#merge the for loop data set (dat_sum) and the newly merged MR file
at.full<-merge(all.trait, dat.sum, by="merge.index", all=T)
at.short<-merge(all.trait, dat.sum, by="merge.index")

#write new csv file and save to box
write.csv(at.full, "AA_all_cases.csv")
write.csv(at.short, "AA_MR_cases.csv")


##################################################################
##################################################################
#Calculating Q10 vals for model + Graphics
install.packages("respirometry")
library(respirometry)

AA_MR <- read.csv("AA_MR_cases.csv", header = TRUE)

#Evolved MR: for each site origin, what is the response of MR across the urbanization gradient?
Q10_val <- Q10(R1 = , R2 = , T1 = 4, T2 = 10) #returns Q10; = 


Q10.vals<- Q10(T_vec=tmp$MRTestTemp, R_vec=tmp$MeanMR)
#model 1 w/ MR as a function of test temp
mod<-lm(MeanMR ~ MRTestTemp, data=tmp)

#model 2 with MR as a function of site origin
mod2 <- lm(MeanMR ~ Site, data = dat1)





