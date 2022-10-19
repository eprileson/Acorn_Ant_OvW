##### OvW_CCRT Data Analysis ###
##### Acorn Ant Overwintering: Cold Tolerance across the gradient ####
####

#Packages and Libraries
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggeffects)
library(car)
library(visreg)
library(nlme)
library(lme4)
library(glmmTMB)
library(magrittr)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(puniform)
library(knitr) 
library(readxl)

setwd("C:/Users/prile/Box/Research/AcornAntOverwintering2021/CCRT")

### Part 2: Data Wrangling: 

#read in data, take a look
ccrt <- read.csv("ccrt.csv", header = TRUE)
head(ccrt)

#check classes, make sure treatment is a factor (2 levels: rural, urban)
class(ccrt$Treatment) #chr; change to factor, then recheck class and levels
ccrt$Treatment <- as.factor(ccrt$Treatment)
levels(ccrt$Treatment) #may be extra level for blank rows; use droplevels() to exclude

#change date category to factor
ccrt$Collection_date <- as.factor(ccrt$Collection_date)
class(ccrt$Collection_date)
levels(ccrt$Collection_date) #early, mid, late

#change name to col season
colnames(ccrt) <- c("Colony_ID", "Ind_ant", "Site", "Date", "CCRT_min", "Treatment", "Col_Season", "Worker_count")
head(ccrt)

#change date to ymd format for R
ccrt$Date <- ymd(ccrt$Date)
ccrt$Collection_date <- ymd(ccrt$Collection_date) #only if we don't use factor for date

#change ccrt to seconds, change name, round digits to nearest whole #
ccrt$CCRT_sec <- ccrt$CCRT_min*60
ccrt$CCRT_sec <- round(ccrt$CCRT_sec, digits = 0)
#change column names to include seconds as new variable
colnames(ccrt) <- c("Colony_ID", "Ind_ant", "Site", "Date", "CCRT_min", "Treatment", "Col_Season", "Worker_count", "CCRT_sec")
class(ccrt$CCRT_sec)
ccrt$CCRT_sec <- as.integer(ccrt$CCRT_sec) #change ccrt_sec to integer
head(ccrt)

#inspect distribution of ccrt / continuous response variable distribution
#since only 1 continuous response variable, use histogram
par(mfrow= c(1,2)) #show dual panels
hist(ccrt$CCRT_sec) #large rightward skew;
#try log transformed
hist(log(ccrt$CCRT_sec)) #little better, but now skew is flipped, so still not norm

#remove colonies with sample size <= 5
ccrt <- ccrt[-c(42,188:189,294,552:555),] #colony 6951 U, colony 5619, colony 5260, colony 7783

### Part 3: Exploratory Graphics ##
#make initial boxplot to explore relationship of treatment on chill coma recovery time
par(mfrow = c(1,1)) #revert back to widescreen
boxplot(CCRT_sec ~ Treatment, data = ccrt)  #wide variance for both, but not much difference
#between treatment groups initially since the boxes / IQRs overlap, slightly lower median in rural

#make graph showing differences in each site using color as third variable
ggplot(ccrt, aes(x = Treatment, y = CCRT_sec)) +
  geom_point()+
  labs(x = 'Treatment', y = 'Chill Coma Recovery Time') +
  ggtitle('Cold Tolerance Across Urban and Rural Sites')+
  scale_colour_brewer(palette = 'Set1') +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

##calculate means of ind colonies to compare to overall avgs
#then plot to show colony ID as a function of CTmax w/ treatment as color
#save this raw data for the final figure
col_means2 <- aggregate(CCRT_sec ~ Colony_ID + Treatment+Col_Season+Site+Worker_count,
                       data = ccrt, mean)
#round to nearest whole
col_means2$CCRT_sec <- round(col_means2$CCRT_sec, digits = 0) #round digits to nearest 0.5
head(col_means2)
#change site to factor; levels = R1, R2, U1, U2, U3
col_means2$Site <- as.factor(col_means2$Site)
class(col_means2$Site)

#summary of raw data for CCRT colony means
summary(col_means2[c(3,7:8,10:11,15:16,26,29:32,35:41),6]) #meanU = 650.4, medianU = 732, min = 61, max = 1344
summary(col_means2[c(1:2,4:6,9,12:14,17:25,27:28,33:34),6]) #mean R = 646.9, median R = 644.5,min = 72, max = 1585

#create a plot with colony averages so that you can use for data viz later to show
#the variance amongst all sites
first1 <- ggplot(col_means2, aes(x = Treatment, y = CCRT_sec, shape = Site))+
  geom_jitter(width = 0.25, height = 0)+
  scale_color_manual(values = mySite_colors)

  
#use full data set (ccrt) for grand means w/ colony means as jittered pts
first1 + geom_point(data = ccrt, aes(x = Treatment, y = CCRT_sec), size = 7,colour = my_colors1, 
                    stat = 'summary', fun.y = "mean", inherit.aes = FALSE)+
  geom_errorbar(colour = my_colors1, data = ccrt, aes(x = Treatment, y = CCRT_sec), stat = 'summary', fun.data = 'mean_se', 
                width=0.05,  fun.args = list(mult = 1), inherit.aes = FALSE)+
  labs(y = "CCRT (sec))", x = "Treatment")+
  ggtitle('Cold Tolerance Across Urban and Rural Sites')+
  theme_classic()+
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#create colors for Sites
mySite_colors <- c("#2196F3","#4CAF50","#F44336", "#FA28FF", "#FFEB3B")
names(mySite_colors) <- levels(col_means2$Site)

### Part 4: Initial model construction ###

#construct generalized mixed effect model using glmer() with gaussian distribution since the data 
# fit a poisson dist with seconds as count data, log link along with
#random effects of colony ID due to non-independence of those factors and nested
#collection date's interaction with source to control for possible interactions

#make full GLMM model w/ glmer
mod_ccrt1 <- glmmTMB(CCRT_sec ~ Treatment + (1 | Colony_ID), 
                   data = ccrt, family = genpois(link = "log"))

#model refit 1 - FINAL MODEL - use more flexible glmmTMB
mod_control1 <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")) #may not be needed
mod_ccrt2 <- glmmTMB(CCRT_sec ~ Treatment + (1 | Colony_ID) + (1 | Col_Season),
                     data = ccrt, family = genpois(link = "log"))
##original model: glmmTMB(CCRT_sec ~ Treatment + (1 | Colony_ID) + (1 | Collection_date*Treatment))

mod_ccrt3 <- glmmTMB(CCRT_sec ~ Treatment + (1 | Colony_ID) + (1 | Col_Season),
                     data = ccrt, family = poisson(link = "log"))
#FINAL MODEL following model comparison and diagnostics
mod_ccrt4 <- glmmTMB(CCRT_sec ~ Treatment + (1 | Colony_ID) + (1 | Col_Season),
                     data = ccrt, family = nbinom1(link = "log"))


#plot model to see if it matches expl graphics
library(visreg)
visreg(mod_ccrt3)


### Part 5: Model diagnostics ###
diagnose(mod_ccrt4)
library(DHARMa) #use simulated diagnostic modeling since we have a glmm; works with both glmer and glmmTMB objects

#second simulation of refit model diagnostics
simOutput2 <- simulateResiduals(fittedModel = mod_ccrt4, plot = F)
plot(simOutput2)  #less dispersion, so better family, but deviation still significant
plotResiduals(simOutput2, form = ccrt$Treatment)
testOutliers(simOutput2, margin = c("both"), type = c("bootstrap"), plot = T) #outliers are signf. here

#check spread of residuals of each model using base plot
plot(residuals(mod_ccrt4))
abline(h = 0)

#check normality of residuals of each model using base qqplot
q2 <- qqnorm(residuals(mod_ccrt4))

#each model looks pretty decent in terms of normality; go with models #2 and 4 since it contains
#the additional fixed predictor

#use AIC to compare and select model
AIC(mod_ccrt2, mod_ccrt4) #has lowest AIC value
AIC(mod_ccrt1) # AIC value is higher

### Part 6 - Statistical Testing ###

summary(mod_ccrt4) #estimate: 0.0009777(greater than Rural), SE = 0.1676475

#hypothesis testing using analysis of deviance
library(car)
Anova(mod_ccrt4, type = "III") #Treatment 
#chisq -> Treatment = 0.0, p = 0.9953

#emmeans pairwise analysis - use to backtransform data since the data was log transformed in model
library(emmeans)
bt_model <- emmeans(mod_ccrt4, pairwise ~ Treatment, transform = "response", adjust = "fdr")
summary(bt_model) #rural treatment estimate = 472 (SE = 89.8), 
#urban treatment is 482 (SE = 102.8), more than the rural 
#contrasts: Rural - Urban = -.563, SE = 96.6


### Part 7 - Plotting model output ###

#make object from emmeans that can work w/ ggplot for figure
library(magrittr)
bt_df <- bt_model$emmeans %>%
  confint() %>%
  as.data.frame()

#make upper and lower SEs (manually calculated from summary of emmeans)
upper.SEc <- c(562.1, 579.8)
lower.SEc <- c(383.9, 384.2)

#graph w/ emmeans model (to get error bars that match accurately from the backtransformed model)
ggplot(bt_df, aes(x = Treatment, y = rate))+
  geom_point(colour = my_colors1, size = 7,
             stat = 'summary', fun.y = 'mean')+
  geom_errorbar(aes(ymin = lower.SEc, ymax = upper.SEc), #uses upper and lower CI
                width=0.05, colour = my_colors1, fun.args = list(mult = 1))+
  geom_jitter(data = col_means2, aes(x = Treatment, y = CCRT_sec),
              width = 0.2, height = 0, size = 2, inherit.aes = FALSE)+
  labs(y = "CCRT (sec)", x = "Source Population")+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#plot modeled relationships (alternate)
plot_mod1 <- ggpredict(mod_ccrt2, terms = c("Treatment", "Collection_date")) #predicted values
plot_mod1base <- ggpredict(mod_ccrt2, terms = c("Treatment"))

#base plot w/o jittered colony averages
plot(plot_mod1base)+
  geom_point(colour = my_colors1, size = 4, 
             stroke = 1, stat = 'summary', fun.y = 'mean')+
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0.05, colour = my_colors1, fun.args = list(mult = 1))+
  labs(y = "CCRT (sec)", x = "Treatment")+
  ggtitle('Cold Tolerance Across Urban and Rural Sites')+
  theme_classic()+
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

#plot graph - combine modeled total averages w/ jittered colony means and grand mean dot errors w/ ses
ggplot(plot_mod1, aes(x, y = predicted))+
  geom_point(colour = my_colors1, size = 7,
            stat = 'summary', fun.y = 'mean')+
  geom_errorbar(width=0.05, colour = my_colors1, fun.args = list(mult = 1))+
  geom_jitter(data = col_means2, aes(x = Treatment, y = CCRT_sec, color = Collection_date),
              width = 0.2, height = 0, inherit.aes = FALSE)+
  labs(y = "CCRT (sec)", x = "Treatment")+
  ggtitle('Cold Tolerance Across Urban and Rural Sites')+
  annotate("text", x = 2, y = 1600, label = "NS")+
  theme_classic()+
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
#colors to use for graphics
my_colors1 <- c("cadetblue", "darkorange")
names(my_colors1) <- levels(col_means2$Treatment) #matches colors to Rural and Urban

##
### END ###
##



####
##### look at interaction of collection date on colony size  #####
#quick look
class(ccrt$Worker_count)

#convert date to Julian date of collection - changes to the number day out of 365
#use yday() function from lubridate
library(lubridate)
ccrt$Collection_date <- yday(ccrt$Collection_date)

#quick look
plot(Worker_count ~ Collection_date, data = ccrt)

#basic linear model 
library(nlme)
int_mod <- lme(Worker_count ~ Collection_date*Treatment, random = ~1 | Colony_ID, data = ccrt)
library(glmmTMB)
int_mod1 <- glmmTMB(Worker_count ~ Collection_date*Treatment + (1 | Colony_ID),
                data = ccrt, family = poisson(link = "log"))


#stats test and summary
summary(int_mod)  # adjusted R squared = 0.004. Est = 0.07 +/- 0.041 SE
anova(int_mod) # F = 9.22, p = 0.0043

library(car)
Anova(int_mod, type = "III")

#emmeans to find diff between months
library(emmeans)
int_em <- emmeans(int_mod, pairwise ~ month, transform = "response", adjust = "fdr")
int_em #significant diffs between each month except for nov and dec

#visualize
visreg(int_mod) #slight positive correlation, but very wide variance, not great fit
#model predictions with both predictor terms
int_graph <- ggpredict(int_mod, terms = c("Collection_date", "Treatment"), ci.lvl = 0.68, type = "fe")
int_graph$group
#visualize / graph output of model predictions
plot(int_graph)+
  geom_smooth(aes(x, y = predicted, color = group))+
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  scale_fill_manual(values = c("cadetblue", "darkorange"))+
  labs(y = "Worker Number", x = "Date", color = "Source Population")+
  ggtitle('Starting Worker Number Across Collection Period')+
  theme_classic()+
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

mycolorsW <- c("cadetblue", "darkorange")
names(mycolorsW) <- levels(ccrt$Treatment)

