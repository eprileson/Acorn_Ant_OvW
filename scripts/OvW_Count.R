################
##################### Acorn Ant Overwintering #####
### Count Data ###
##demographic change / survival count post winter ##

#Packages and Libraries
library(dplyr)
library(tidyr)
library(plyr)
library(lubridate)
library(ggplot2)
library(ggeffects)
library(emmeans)
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
library(cowplot)

#read in data set
setwd("C:/Users/prile/Box/Research/AcornAntOverwintering2021/Acorn_Ant_OvW/Acorn_Ant_OvW/scripts")
count <- read.csv("OvW_count.csv")
head(count)

##Part 2: Data wrangling
#make sure source.pop is factor, check levels
count$Source.pop <- as.factor(count$Source.pop)
class(count$Source.pop)
levels(count$Source.pop)
#R doesn't choose End first

#adjust collection date to col_season factor
count$Collection_date <- as.factor(count$Collection_date)
class(count$Collection_date) #changed to date that R recognizes

#change name to col season
colnames(count) <- c("Colony_ID", "Source.pop", "Worker.start", "Worker.end", "Brood_count", "Col_Season")
head(count)

#create new column of worker demographic change, calculate change in pop. Start - End
count$Worker.loss <- (count$Worker.start - count$Worker.end) #make worker loss positive for count
head(count)

#create new column of worker number change as a proportion
count$worker.prop <- ((count$Worker.end)/(count$Worker.start))

###Part 3: Exploratory Graphics: initial look at data
plot(Worker.loss ~ Source.pop, data = count) #change in worker numb. ~ source pop
plot(worker.prop ~ Source.pop, data = count) #no relationship it seems
plot(Collection_date ~ Source.pop, data = count) #rural collected at start and end; urban in middle period

#shows urban pops have had bigger change in survival than rural
hist(count$Worker.loss)

####Alternate method: wrangle data to long form to show pop change over time on X ~ source
count1 <- read.csv("OvW_count.csv")
head(count1)
#change column name so that the starting pop (beg) is on left and time goes left to right
colnames(count1) <- c("ColonyID", "Source.pop", "Beginning", "End", "Brood_count","Collection_date")

#make source pop a factor & check levels
count1$Source.pop <- as.factor(count1$Source.pop)
levels(count1$Source.pop)

#adjust collection date to col_season factor
count1$Collection_date <- as.factor(count1$Collection_date)
class(count1$Collection_date) #changed to date that R recognizes

#change name to col season
colnames(count1) <- c("Colony_ID", "Source.pop", "Beginning", "End", "Brood_count", "Col_Season")
head(count1)

#pivot the data longer so that we can show the change over time with the data set
count1 <- count1 %>%
  pivot_longer(c("Beginning", "End"), names_to = "Time", values_to = "Worker.pop")

count1$Time <- as.factor(count1$Time)
class(count1$Time)


#exploratory graphics: (worker pop ~ time*source.pop)
ggplot(count1, aes(x = Time, y = Worker.pop, colour = Source.pop))+
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0.05, fun.args = list(mult = 1))+
  geom_point(aes(colour = factor(Source.pop)), size = 7, stat = 'summary', fun.y = "mean")+
  geom_segment(aes(x = 1, y = 45.5, xend = 2, yend = 13), color = "cadetblue")+
  geom_segment(aes(x = 1, y = 62, xend = 2, yend = 15), color = "darkorange")+
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  theme_classic()+
  labs(y = "Worker Population", x = "Census Points", color = "Population")+
  ggtitle('Colony Population Change')+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#calculate median starting pop urban v rural
#base
medR <- median(count[c(1:9, 31:42, 45:47),3]) #rural median = 48.5
sd(count[c(1:9, 31:42, 45:47),3]) #rural sd = 26.325
mean(count[c(1:9, 31:42, 45:47),3]) #rural mean = 45.25
#median urban start
median(count[c(10:30, 43:44),3]) #urban median = 71
sd(count[c(10:30, 43:44),3]) #urban sd = 39.83
mean(count[c(10:30, 43:44),3]) #urban mean = 61.83

#median and mean rural end
median(count[c(1:9, 31:42, 45:47),4]) #rural median = 11.5
sd(count[c(1:9, 31:42, 45:47),4]) #rural sd end = 7.095
mean(count[c(1:9, 31:42, 45:47),4]) #rural mean = 12.21
#median urban end
median(count[c(10:30, 43:44),4]) #urban median = 10
sd(count[c(10:30, 43:44),4]) #urban sd end = 11.27
mean(count[c(10:30, 43:44),4]) #urban mean = 14.04


##Part 4: Model development
#first model for worker loss over 4 mo
mod.control2 <- glmmTMBControl(optimizer=optim,
                               optArgs=list(method="BFGS"))
surv.mod <- glmmTMB(Worker.loss ~ Source.pop + (1 | ColonyID) + (1 | Collection_date:Source.pop), 
                    data = count, control = mod.control2, 
                    family = poisson(link = "log"))

#initial model: Proportional Worker Loss
prop.surv <- glmmTMB(worker.prop ~ Source.pop + (1 | Colony_ID) + (1 | Col_Season),
                    data=count)
#betareg try?
prop.surv2 <- betareg(worker.prop ~ Source.pop, random =  ~1 |Col_Season,
                      weights= Worker.start, link = 'log', data=count) #no
#model with lme built in?
prop.surv3 <- glmmTMB(worker.prop ~ Source.pop, random = ~ 1 | Col_Season,
                      weights= Worker.start, family= quasibinomial(link = 'log'), data=count) #no

#FINAL MODEL 1: Proportion of worker ants remaining after winter
prop.surv1 <- glmmTMB(worker.prop ~ Source.pop + (1 | Col_Season),
                      weights= Worker.start, family = betabinomial(link = "logit"), data=count)

#FINAL MODEL 2: Change model (w/ both end and beg)
#added the colony ID random effect back
change.mod <- glmmTMB(Worker.pop ~ Time*Source.pop + (1 | Colony_ID) + (1 | Col_Season),
                      data = count1, family = poisson(link = "log"))


#plot basic modeled relationship to see if it matches expl graphics
visreg(prop.surv1)

## Part 5: Model diagnostics
#QQplots, residuals
##Diagnostics: warning msgs for models - show model convergence problem w/ small eigen value problems
diagnose(change.mod) #removed the interaction effect due to extremely large SD

diagnose(prop.surv1)  #really good fit both w/ DHARMa and diagnose()

simOutput3 <- simulateResiduals(fittedModel = prop.surv1, plot = F)
plot(simOutput3)

testOutliers(simOutput3, margin = c("both"), type = c("bootstrap"), plot = T) #outliers are not signf. here


##Part 6: Statistical Testing
#gather summary data on model, determine model fit
summary(prop.surv1) #urban = -.093 lower proportion workers left SE = 0.19
summary(change.mod) #time = -1.31, urbanpop = 0.226, time*urbanpop = -0.172

# log lik statistic tests for model significance
Anova(prop.surv1, type = "III") #chisq = 0.236; p = 0.627
Anova(change.mod, type = "III") #source pop = chisq = 67.34, p = 0.399
#time = 396.03, p < 2e-16***, sourcepop*time = chisq = 3.64, p = 0.056

#backtransform log-linked data:
#using emmeans for figure and backtransformed values
#emmeans
change.em <- emmeans(change.mod, pairwise ~ Time*Source.pop, transform = "response", adjust = "fdr")
change.em
##Part 7: Plotting Modeled Output
#
#Make follow up plot of modeled results
prop.graph <- ggpredict(prop.surv1, terms = c("Source.pop"), 
                        type = "fe", allow.new.levels = TRUE,ci.lvl = 0.68, width=0.05, title = FALSE)
  

#make emmeans objects into dfs for visualization
#for survival overall

#for proportions - use summary data since not on Log scale
prop_em <- emmeans(prop.surv1, pairwise ~ Source.pop, transform = "logit", adjust = "fdr")
upper.SEp <- c(0.29, 0.27)
lower.SEp <- c(0.22, 0.21)

#for change data
change_df <- change.em$emmeans %>%
  confint()%>%
  as.data.frame()
upper.SEg <- c(43.12, 11.7, 54.38, 12.41)
lower.SEg <- c(28.92, 7.74, 35.94, 8.11)


#plot colors
my_colors3 <- c("cadetblue", "darkorange")
names(my_colors3) <- levels(count$Source.pop)

my_colors4 <- c("cadetblue", "cadet blue", "dark orange", "dark orange")
names(my_colors4) <- levels(c(count1$Source.pop, count1$Time))

#plot graphs
#change
change <- ggplot(change_df, aes(x = Time, y= rate))+  
  geom_point(data = count1, aes(x = Time, y = Worker.pop, color = Source.pop), 
             alpha = 0.5, position = position_jitterdodge(jitter.width = 0.5, dodge.width = NULL), inherit.aes = FALSE)+
  geom_point(aes(fill = factor(Source.pop)), color = my_colors4, size = 7, position = position_dodge(width = 0.3))+
  geom_errorbar(aes(fill = factor(Source.pop), ymin = lower.SEg, ymax = upper.SEg), color = my_colors4, fun.data = 'mean_se', 
                width=0.05, fun.args = list(mult = 1), position = position_dodge(width = 0.3))+
  scale_color_manual(values = c("cadet blue","dark orange"))+
  guides(fill = "none", color = "none")+  #change back if need legend
  labs(y="Worker Number", x = "Census Points", tag = "A", color = "Source Population")+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    plot.tag = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#workers lost as a proportion of starting size
prop <- ggplot(prop.graph, aes(x, y = predicted))+
  geom_point(data = count, aes(x = Source.pop, y = worker.prop, color = Source.pop), 
             position = position_jitterdodge(jitter.width = 0.5, dodge.width = NULL), inherit.aes = FALSE)+
  geom_point(colour = my_colors3, size = 7)+
  geom_errorbar(aes(ymin = lower.SEp, ymax = upper.SEp), fun.data = 'mean_se', 
                width=0.05, fun.args = list(mult = 1), colour = my_colors3)+
  scale_color_manual(values = c("black","black"))+
  labs(y="Proportion of Workers Remaining", x = "Source Population", tag = "B")+
  guides(color = "none")+  #change back if need legend
  theme_classic()+
  theme(
    title = element_text(size = 14),
    plot.tag = element_text(size = 16, hjust = 0.5),
    axis.text = element_text(size = 10)
  ) 

#make graphs into combined facets
change
total
prop
plot_grid(change, prop, ncol = 2, nrow = 1, 
          rel_widths = c(.5, .5, .5))   #use 1028 x 555 px aspect ratio




#####################
################### END

