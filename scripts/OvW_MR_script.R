#######
############ Acorn Ant OvW - MR Batch Processing ############
###

#install / load packages for reading exp files from Sable ExpeData
install.packages("remotes")
remotes::install_github("hawkmoth/sablebase")
install.packages("zoo")
library(zoo)
library(SableBase)

library(lubridate)
library(car)

setwd("C:/Users/prile/Box/Research/AcornAntOverwintering2021/MetabolicRate/Eric_OvW")

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
dat_ovw.3c[1:25,]


# list all files in the folder / from the working dir
files = dir(pattern="exp") #find all files in directory that have "exp"

files = files[files %in% c(dat_ovw.3c$MR4_name, dat_ovw.3c$MR10_name)==TRUE] #search files only in these columns

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
      
      # sscf.sub$fracflow = sscf.sub$fracCO2*500 #
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
      # This should be the last 300 readings since there are 10 minutes total in the trial
      
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

dat.sum[1:27,]

#make new csv file w/ mean MR data
setwd("C:/Users/prile/Box/Research/AcornAntOverwintering2021/MetabolicRate/Eric_OvW")
write.csv(dat.sum, "AA_MR_Summary_Stats.new.csv")

## did not do, additional experimental data wrangling
setwd("C:/Users/sarah/Box/DiamondLabResearch/Research/VCardui_Ontogeny_2021/SED_Analyses")

# write.csv(dat.sum, "VCB_MR_Summary_Stats.csv")
#make new column combines data file and and chamber index #
dat.sum$merge.index<-paste(dat.sum$ExpeDataFile, dat.sum$ChamberIndex, sep="_")

dat_ovw.3c<-read.csv("AA_labels_edits.csv", header=T)
head(dat_ovw.3c)
dat_ovw.3c4 <-dat_ovw.3c[,c(1:9)]
dat_ovw.3c$MRTestTemp <-rep(4, times=dim(dat_ovw.3c4)[1])
colnames(dat_ovw.3c4) <-recode(colnames(dat_ovw.3c4), "'MR4_file'='MR_file';'MR4_name'='MR_name'")
dat_ovw.3c10 <-dat_ovw.3c[,c(1:9)]
dat_ovw.3c$MRTestTemp<-rep(10, times=dim(dat_ovw.3c10)[1])
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
## MR data analysis -

#install and load packages
library(ggplot2)
library(nlme)
library(glmmTMB)
library(car)
library(DHARMa) 
library(visreg)
library(ggeffects)
library(dplyr)
library(emmeans)

## Part 1: load data and set working dir

setwd("C:/Users/prile/Box/Research/AcornAntOverwintering2021/Acorn_Ant_OvW/Acorn_Ant_OvW/scripts")
AA_MR <- read.csv("AA_MR_Summary_Stats.csv", header = TRUE)
head(AA_MR)

###Part 2: data wrangling continued
#remove colonies not tested in CCRT
AA_MR <- AA_MR[-c(7,15,38,44,64,71,120,126),]

#change source pop and coldate to factor, check levels for R and U
AA_MR$Source.pop <- as.factor(AA_MR$Source.pop)
levels(AA_MR$Source.pop)
AA_MR$Collection_date <- as.factor(AA_MR$Collection_date)

#change col names to reflect cod date as factor
colnames(AA_MR) <- c("Sec", "MeanMR", "MeanMR_RA", "Cor_SumMR", "Chamber", "ExpeDataFile",
                      "DataLength","StartIndex", "Source.pop", "Colony_mass", "Colony_ID",
                      "facet","Q10","Col_Season")  #changed test temp to facet for visualization help later
head(AA_MR)


#change temp to factor
AA_MR$facet <- as.factor(AA_MR$facet)
class(AA_MR$facet)
levels(AA_MR1$facet)

#remove NAs from the dataset; note, this removes the chamber 8 control tubes
AA_MR1 <- AA_MR[!(AA_MR$Chamber ==8),]

#check
head(AA_MR1)

#log transform MR and colony mass data (don't do for Q10 data)
AA_MR1$LogMeanMR <- log(AA_MR1$MeanMR)
AA_MR1$LogColony_mass <- log(AA_MR1$Colony_mass)

#Part 3: exploratory graphics
#
hist(AA_MR$Q10) #pretty normal
hist(AA_MR1$Log.10MeanMR) #pretty good, log transformed approaches normal

#basic graph to show trend of mass on mean MR w/ source pop
#removed NA's using the notation of only including non-NAs in the dataset (see first ggplot arg)
ggplot(AA_MR1[!is.na(AA_MR1$Test.Temp),], aes(x = Log.10Colony_mass, y = Log.10MeanMR, color = Source.pop))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, level = 0.95)+
  facet_wrap(~Test.Temp)+
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  theme_classic()+
  ylab(bquote("Log 10 Mean CO"[2]))+ 
  xlab("Colony Mass Log (grams)")+
  labs(color = "Source Population")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#basic plot to show Q10s
ggplot(AA_MR1[!is.na(AA_MR1$Test.Temp),], aes(x = Log.10Colony_mass, y = Q10, color = Source.pop))+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, level = 0.68)+
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  theme_classic()+
  labs(y = "Q10 Reaction Rate", x = "Colony Mass (Log grams)", color = "Population")+
  ggtitle('Q10 Reaction Rate Across Source Population')+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#mean x values for segments
# y1 (rural / 4 deg = -3.679414)
# y2 (rural / 10 deg = -3.585731)
# y3 (urban / 4 deg = -3.661369)
# y4 (urban / 10 deg = -3.575597)

#colors for data viz
my_colorsMR <- c("cadetblue", "darkorange")
names(my_colorsMR) <- levels(AA_MR1$Source.pop)
my_colorsMR1 <- c("cadetblue", "darkorange")
names(my_colorsMR1) <- levels(mr_df$Source.pop)

#basic plot to show different test temps and mean MR
ggplot(AA_MR1[!is.na(AA_MR1$Test.Temp),], aes(x = Test.Temp, y = Log.10MeanMR, color = Source.pop)) +
         geom_point(aes(colour = factor(Source.pop)),size = 7, stat = 'summary', fun.y = 'mean') +
         geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                       width=0.05, fun.args = list(mult = 1))+
         scale_colour_manual(values = c("cadetblue", "darkorange"))+
         theme_bw()+
         labs(y = "Mean CO2 Log (ppm)", x = "Test Temp", color = "Source Population")+
         ggtitle('Metabolic Rate Across Source and Temperatures')+
         theme_classic()+
         theme(
             title = element_text(size = 14),
             axis.title = element_text(size = 14),
             axis.text = element_text(size = 10)
           )

         geom_segment(aes(x = 4, y = -3.679414, xend = 10, yend = -3.585731))+
           geom_segment(aes(x = 4, y = -3.661369, xend = 10, yend = -3.575597))+
           

##Part 4: Model Construction
#Focal Model: MR ~ mass + source*temp
mod_MR<-lme(LogMeanMR ~ LogColony_mass + Source.pop*Test.Temp, random = ~1|Colony_ID, data=AA_MR1)

#FINAL models
modMR_control <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
mod_MR1 <- glmmTMB(Q10 ~ LogColony_mass + Source.pop + 
                     (1 | Col_Season),
                   data = AA_MR1, family = gaussian(link = "identity"))
mod_MR2 <- glmmTMB(LogMeanMR ~ LogColony_mass + Source.pop*facet +
                     (1 | Col_Season), data = AA_MR1)

###5 Model Diagnostics: 
#use simulated diagnostic modeling since we have a glmm; works with both glmer and glmmTMB objects
simOutput4 <- simulateResiduals(fittedModel = mod_MR1, plot = F)
plot(simOutput4)  #looks good from qqplot and residuals

simOutput5 <- simulateResiduals(fittedModel = mod_MR2, plot = F) #qqplot looks good, but deviations detected

diagnose(mod_MR2) #error w/ small eigen values and large coefficients detected

#plot basic modeled relationship to see if it matches expl graphics
visreg(mod_MR2)


#Part 6: statistical / hypothesis testing
#summary stats
summary(mod_MR2)
Anova(mod_MR2, type="III") #0.645, p = 0.422 #test temp and colony mass are sign. predictors, source pop is not (p = 0.8118)

#Q10 hypothesis test: Does Q10 differ between urban and rural w/ urban having higher Q10 rate?
summary(mod_MR1) #Q10 model
Anova(mod_MR1, type = "III") #chisq = 0.923, P = 0.337

# Pairwise tests
#Q10 emmeans test, w/ backtransformation
Qmests <- emmeans(mod_MR1, pairwise~LogColony_mass + Source.pop, adjust = "tukey")
summary(Qmests) #contrast = 0.174; rural = 1.22, urban = 1.40

##Part 7: Plotting Model Output
#
#Mean MR plot using ggpredict; modeled 
plot_MRmod1 <- ggpredict(mod_MR2, terms =c("LogColony_mass", 
                                           "Source.pop", "facet"), 
                         ci.lvl = 0.95, colors = c("cadetblue", "darkorange"))

#FINAL PREDICTED mean MR plot w/ raw data avg MR values
#raw
#variable names
variable_names <- list("4" = "4°C","10" = "10°")
#first try; not used, since label is off
plot(plot_MRmod1, show.title=F, alpha = 0.10)+
  geom_point(data = AA_MR1, aes(x = LogColony_mass,y = LogMeanMR, color = Source.pop), alpha = 0.5, inherit.aes = FALSE)+
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  theme_classic()+
  ylab(bquote("Metabolic Rate (Ln"*" CO"["2"]*" ppm)"))+
  xlab(bquote("Colony Mass (Ln"*" grams)"))+ 
  labs(title = " ", color = "Source Population")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )+facet_grid(~ facet, labeller = as_labeller(fac_labels))

#try #2 w AL's help FINAL PLOT!!!
#make quick label df
fac_labels <- c('4' = "4° C", '10' = "10° C")
                        
  ggplot(data = plot_MRmod1,aes(x, y = predicted, group = group))+
    geom_smooth(method = "lm", se = FALSE, aes(colour = group), size = 0.75)+
    geom_ribbon(aes(ymin= conf.low, ymax= conf.high, y= NULL, fill = group), alpha = 0.15)+
    guides(fill = FALSE)+
    geom_point(data = AA_MR1, aes(x = LogColony_mass,y = LogMeanMR, color = Source.pop), alpha = 0.5, inherit.aes = FALSE)+
    scale_colour_manual(values = c("cadetblue", "darkorange"))+
    scale_fill_manual(values = c("cadetblue", "darkorange"))+ #added to correct the fill color match
    theme_classic()+
    ylab(bquote("Metabolic Rate (ln"*" CO"["2"]*" ppm)"))+
    xlab(bquote("Colony Mass (ln"*" grams)"))+ 
    labs(title = " ", color = "Source Population")+
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 10)
    ) + facet_wrap(~ facet, labeller = as_labeller(fac_labels))
    
    
#Q10 model prediction
plot_Qmod1 <- ggpredict(mod_MR1, terms = c("LogColony_mass", "Source.pop"), ci.lvl = 0.95) #predicted values

##FINAL Q10 Plot##
#plot the predicted model w/ smoothed lines at 95% CI, (Q10 ~ LN(Mass) + Source.pop)
ggplot(plot_Qmod1, aes(x, y = predicted, group = group))+
  geom_ribbon(aes(ymin= conf.low, ymax= conf.high, y= NULL, fill = group), alpha = 0.15)+
  geom_smooth(method = "lm",se = FALSE, aes(colour = group), size = 0.75)+
  guides(fill = FALSE)+
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  scale_fill_manual(values = c("cadetblue", "darkorange"))+ #added to correct the fill color match
  geom_point(data = AA_MR1, aes(x = LogColony_mass, y = Q10, color = Source.pop), alpha = 0.5, inherit.aes = FALSE)+
  ylab(bquote("Acute Thermal Sensitivity of Metabolic Rate"))+
  xlab(bquote("Colony Mass (ln"*" grams)"))+
  labs(color = "Source Population")+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

## Old Q10 plot; bands not same color as FINAL plot

plot(plot_Qmod1, show.title=F, facet = FALSE, alpha = 0.15, colors = c("cadet blue", "dark orange"))+
  geom_point(data = AA_MR1, aes(x = LogColony_mass, y = Q10, color = Source.pop), alpha = 0.5, inherit.aes = FALSE)+
  ylab(bquote("Acute Thermal Sensitivity of Metabolic Rate"))+
  xlab(bquote("Colony Mass (ln"*" grams)"))+
  labs(color = "Source Population")+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#Use for italics:  ylab(bquote(italic("Q")["10"]*" Reaction Rate"))+
  


###### additional ###  
###Not used for OvW testing ###
  
#Calculating Q10 vals for model + Graphics
  install.packages("respirometry")
  library(respirometry)
  
  
#Evolved MR: for each site origin, what is the response of MR across the urbanization gradient?
  Q10_val <- Q10(R1 = , R2 = , T1 = 4, T2 = 10) #returns Q10; = 
  Q10_calc <- (R1/R2)^(10/R2-R1)
  
#for loop to try and calc Q10 and add new column to csv file

  j.indices1 <- unique(AA_MR$Chamberindex) #make object w/ unique rows for each unique chamber from within each file
  
  for(j in 1:length(j.indices1)) {
    
    if(AA_MR$Test.temp == 4){
      j = R1
    }
    else {
      R2
    }
    R2 <- AA_MR$Test.temp ==10
    }
    AA_MR$Q10_val <- Q10(R1 = , R2 = , T1 = 4, T2 = 10)
    
    
    
  }
  
  Q10.vals<- Q10(T_vec=tmp$MRTestTemp, R_vec=tmp$MeanMR)
  
  
  
#Establish baseline / control level
  
  # list all files in the folder
  files = dir(pattern="exp")
  
  files = files[files %in% c(dat_ovw.3c$MR4_name, dat_ovw.3c$MR10_name)==TRUE]
  
  
  pdf("Control_MR_plots_timed.pdf", width=15, height=15)
  par(mfrow=c(9,8), mar=c(1,1,1,1)) # count of number of files and divide into a grid;  n = 70 files
  
  for (i in 1:length(files)) { 
    
    # read in files 
    sscf = read.sscf(files[i])
    
    # assign to dataframe object
    sscf.file = data.frame(sscf, check.names=TRUE)
    
    
    baseline.run = sscf.file[sscf.file$I_O_N==8,]
    
    colnames(baseline.run)[3]<-"CO2.ANIMAL"
    baseline.run$fracCO2 = baseline.run$CO2.ANIMAL/1000000 
    
    flowrate.4 = match(dat_ovw.3c$MR4_name, files[i])
    flowrate.10 = match(dat_ovw.3c$MR10_name, files[i])
    
    if(sum(flowrate.20, na.rm=T) >=1 ) {
      flowrate = dat_ovw.3c$Gilibrator40[which(is.na(flowrate.4)==F)[1]]
    }
    
    if(sum(flowrate.30, na.rm=T) >=1 ) {
      flowrate = dat_ovw.3c$Gilibrator10[which(is.na(flowrate.10)==F)[1]]
    }
    
    baseline.run$fracflow = baseline.run$fracCO2*flowrate
    
    baseline.run$index = as.numeric(which(sscf$I_O_N==8)) 
    
    baseline.slope.mod = lm(fracflow ~ index, data=baseline.run)
    
    plot(fracflow ~ index, data=baseline.run, type='l')
    
    abline(baseline.slope.mod, col="blue", lwd=2)
  }
  dev.off()
  

