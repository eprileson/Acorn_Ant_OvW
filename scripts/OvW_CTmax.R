##### Acorn Ant Overwintering: Heat Tolerance across the gradient ####
#### CT max data

#Packages and Libraries
install.packages("reporter")
library(reporter)
library(dplyr)
library(ggplot2)
library(lubridate)
library(car)
library(nlme)
library(lme4)
library(glmmTMB)
library(magrittr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(puniform)
library(knitr) 
library(readxl)

setwd("C:/Users/prile/Box/Research/AcornAntOverwintering2021/Acorn_Ant_OvW/Acorn_Ant_OvW/scripts")

#read in data, take a look
ctmax <- read.csv("OvW_ctmax.csv", header = TRUE)
head(ctmax)

#### Part 2: Data Wrangling / Organization

#check classes, make sure treatment is a factor (2 levels: rural, urban)
class(ctmax$Treatment) #chr; change to factor, then recheck class and levels
ctmax$Treatment <- as.factor(ctmax$Treatment)
levels(ctmax$Treatment) 

#check response var class; 
class(ctmax$CTmax_C)
class(ctmax$Worker_count)
ctmax$Worker_count <- as.numeric(ctmax$Worker_count) #change to numeric

#change date to factor for three cohorts
ctmax$Collection_date <- as.factor(ctmax$Collection_date)
levels(ctmax$Collection_date) #change date from char to date using lubridate

#change name to col season
colnames(ctmax) <- c("Colony_ID", "Ind_ant", "Site", "Date", "CTmax_C", "Treatment", "Col_Season", "Worker_count", "Winter_days")
head(ctmax)

#remove colony 2318, 7783, and 3538 from analysis (sample size too small, <= 5)
ctmax <- ctmax[-c(11:12,216:219,345:347),] #2318, #3538, #7783

### Part 3: Exploratory Graphics ####
#inspect distribution of continuous response variable distribution
#since only 1 continuous variable, use histogram
par(mfrow= c(1,2)) #show dual panels
hist(ctmax$CTmax_C) #approaches 

#make initial boxplot to explore relationship of treatment on chill coma recovery time
boxplot(CTmax_C ~ Treatment, data = ctmax)  #wide variance for both, but not much difference
#between treatment groups initially since the boxes / IQRs overlap, slightly lower median in rural

#calculate means of ind colonies to compare
#then aggregate the means by their treatment (ie. rural or urban) and include 
#the other variables - use this for data viz to show variance of colony avgs
col_means <- aggregate(CTmax_C ~ Colony_ID + Treatment + Site + Col_Season,
                        data = ctmax, mean)

#round the colmeans to nearest 0.5
library(plyr)
col_means$CTmax_C <- round_any(col_means$CTmax_C, 0.5) #round digits to nearest 0.5

#avg colonies by their CTmax
mean(col_means[23:41,3]) #mean U = 45.517, mean R = 45.053
mean(col_means[1:22, 3])

#change site to factor
col_means2$Site <- as.factor(col_means2$Site)
class(col_means2$Site)

#make graph w/ jittered averaged colonies and grand means w/ SE from full data set
first <- ggplot(col_means, aes(x = Treatment, y = CTmax_C, color = Site))+
  geom_jitter(width = 0.25, height = 0)+
  scale_color_manual(values = mySite_colors1)#makes first graph w/ jittered average colony pts

first + geom_point(data = ctmax, aes(x = Treatment, y = CTmax_C), colour = "red", size = 4, 
                   fill = "white", stroke = 1, stat = 'summary', fun.y = 'mean',inherit.aes = FALSE)+
  geom_errorbar(data = ctmax, aes(x = Treatment, y = CTmax_C), stat = 'summary', fun.data = 'mean_se', 
                width=0.05, fun.args = list(mult = 1), inherit.aes = FALSE)


#create colors for Sites
mySite_colors1 <- c("#2196F3","#4CAF50","#F44336", "#FA28FF", "#FFEB3B")
names(mySite_colors1) <- levels(col_means2$Site)


### Part 4: Initial model construction ###
#need linear mixed effects model w/ colony ID as random effect since individual ants from same
#colony are not independent, ie. avoid autocorrelation. Include Treatment as main mixed effect
#but also include Collection date as fixed predictor as well

library(visreg)
library(lme4)
library(glmmTMB)
library(nlme)
mod <- lme(CTmax_C ~ Treatment + Collection_date, random =~ 1 | Colony_ID, data = ctmax)

mod1 <- lme(CTmax_C ~ Treatment + Collection_date, random =~ 1 | Colony_ID, data = col_means)

#ORIGINAL FINAL MODEL
mod2 <- glmmTMB(CTmax_C ~ Treatment + (1 | Colony_ID) + (1 | Collection_date:Treatment),
                data = ctmax) 

#NEW FINAL MODEL with new random effect of col_season instead of col_date
mod4 <- glmmTMB(CTmax_C ~ Treatment + (1 | Colony_ID) + (1 | Col_Season),
                data = ctmax)

#plot model to see if matches expl graphics
library(visreg)
visreg(mod4)


### Part 5: Model diagnostics ###

#check spread of residuals of each model using base plot
diagnose(mod4) #liklihood ratio tests may still be ok
plot(residuals(mod2))
abline(h = 0)

library(DHARMa) #use simulated diagnostic modeling since we have a glmm; works with both glmer and glmmTMB objects

#second simulation of refit model diagnostics
simOutput3 <- simulateResiduals(fittedModel = mod4, plot = F)
plot(simOutput3)  #less dispersion, so better family, but deviation still significant
plotResiduals(simOutput3, form = ctmax$Treatment)
testOutliers(simOutput3, margin = c("both"), type = c("bootstrap"), plot = T) #outliers are signf. here

#check normality of residuals of each model using base qqplot
qqnorm(residuals(mod4))

#each model looks pretty decent in terms of normality; 

### Part 6 - Statistical Testing ###
#call w /summary 

summary(mod4) #urban est = 1.2796 deg C higher than Rural; SE = 0.5124

library(car)
#hypothesis testing using analysis of deviance
Anova(mod4, type = "III")  #chisq - Treatment = 6.24, p = 0.013, .ns.

library(emmeans) 
#need to do pairwise comparison
ct_mod <- emmeans(mod4, pairwise ~ Treatment, transform = "response", adjust = "fdr")
ct_mod #estimate: TreatmentRural = 45.3, SE= 0.270; TreatmentUrb = +0.06 higher, SE = 0.288
#CollectionDate estimate = +0.021, SE = 0.0063


### Part 7 - Plotting model output ###
#plot modeled relationships

#change the emmeans object to data frame; use for final plot graphics
library(magrittr)
ct_df <- ct_mod$emmeans %>%
  confint() %>%
  as.data.frame()

#make upper and lower SEs; manually calculated from emmeans from mod4
upper.SE <- c(45.468, 46.772)
lower.SE <- c(44.332, 45.628)

#make plot with the emmeans to get correct SEs that have been backtransformed
## FINAL PLOT = ggemmeans_final.jpeg in Box
ggplot(ct_df, aes(x = Treatment, y= emmean))+
  geom_jitter(data = col_means, aes(x = Treatment, y = CTmax_C), 
              width = 0.2, height = 0, size = 2, inherit.aes = FALSE)+
  geom_point(colour = my_colors2, size = 7)+
  geom_errorbar(aes(ymin = lower.SE, ymax = upper.SE), #uses upper and lower CI
                width=0.05, colour = my_colors2, fun.args = list(mult = 1))+
  ylab(bquote(CT[max]*" °C"))+
  xlab("Source Population")+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#alternate graph - using basic predict
library(ggeffects)
plot.mod <- ggpredict(mod2, terms = c("Treatment", "Collection_date"), ci.lvl = 0.68)

#predict function w/ model
plot.mod_base <- ggpredict(mod2, terms = c("Treatment"))
#base plot w/ grand means and error bars; just to show overall diff w/o the colony avgs
plot(plot.mod_base)+
  geom_point(colour = my_colors, size = 7)+
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0.05, colour = my_colors, fun.args = list(mult = 1))+
  geom_jitter(data = col_means, aes(x = Treatment, y = CTmax_C, color = Collection_date), 
              width = 0.2, height = 0, inherit.aes = FALSE)+
  labs(y = "CTmax (C)", x = "Treatment")+
  annotate("text", x = 1.5, y = 45.5, label = "NS")+
  ggtitle('Heat Tolerance Across Urban and Rural Sites')+
  theme_classic()+
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

#plot graph; use the combined ggplot layers w model first from ggpredict, then
#add in the jittered points from the data exploration
library(viridis)

ggplot(plot.mod, aes(x, y = predicted))+
  geom_jitter(data = col_means, aes(x = Treatment, y = CTmax_C, color = Collection_date), 
              width = 0.2, height = 0, inherit.aes = FALSE)+
  geom_point(colour = my_colors, size = 7, 
           stroke = 1, stat = 'summary', width = 0.05, fun.y = 'mean')+
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0.05, colour = my_colors)+
  annotate("text", x = 2, y = 48, label = "NS")+
  labs(y = "CTmax (C)", x = "Treatment")+
  ggtitle('Heat Tolerance Across Urban and Rural Sites')+
  theme_classic()+
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

#add ind colors for ctmax
my_colors2 <- c("cadetblue", "darkorange")
names(my_colors2) <- levels(ctmax$Treatment)

###
#make pub quality graphic - use tiff() function w/ 480 x 480 pixel size and resolution of 300 ppi
tiff("ggemmeans_final.jpeg", width = 480, height = 480, units = 'px', res = 300)

bitmap("ggemmeans_final.jpeg", height = 480, width = 480, type="tifflzw", res=300)
dev.off()
