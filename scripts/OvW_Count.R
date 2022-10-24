################
#####################
##demographic change / survival count post winter ##

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
count$worker.prop <- ((count$Worker.loss)/(count$Worker.start))

###Part 3: Exploratory Graphics: initial look at data
plot(Worker.loss ~ Source.pop, data = count) #change in worker numb. ~ source pop
plot(Worker.loss ~ Collection_date, data = count) #no relationship it seems
plot(Collection_date ~ Source.pop, data = count) #rural collected at start and end; urban in middle period

#shows urban pops have had bigger change in survival than rural
hist(count$Worker.loss)

#Alternate method: wrangle data to long form to show pop change over time on X ~ source
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
library(magrittr)
library(tidyr)
count1 <- count1 %>%
  pivot_longer(c("Beginning", "End"), names_to = "Time", values_to = "Worker.pop")

count1$Time <- as.factor(count1$Time)
class(count1$Time)


#exploratory graphics: (worker pop ~ time*source.pop)
library(ggplot2)
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

#calculate median starting pop urban v rural - TidyR
count <- count %>%
  group_by(Source.pop)%>%
  mutate(med.R = median(Worker.end))%>%
  mutate(sd = sd(Worker.end))

#OR base
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
library(glmmTMB)
library(nlme)

#first model for worker loss over 4 mo
#FINAL MODEL 1 - 
mod.control2 <- glmmTMBControl(optimizer=optim,
                               optArgs=list(method="BFGS"))
surv.mod <- glmmTMB(Worker.loss ~ Source.pop + (1 | ColonyID) + (1 | Collection_date:Source.pop), 
                    data = count, control = mod.control2, 
                    family = poisson(link = "log"))

#FINAL MODEL 2: Proportional Worker Loss
prop.surv <- glmmTMB(worker.prop ~ Source.pop + (1 | Colony_ID) + (1 | Col_Season),
                    data=count)
#model w/ quasi binomial - doesn't work with glmmTMB
library(lme4)
library(MASS)
prop.surv1 <- glmmTMB(worker.prop ~ Source.pop + (1 | Colony_ID) + (1 | Col_Season),
                      weights= Worker.start, family=binomial, data=count)
prop.surv2 <- betareg(worker.prop ~ Source.pop, random = ~1 |Colony_ID, random =  ~1 |Col_Season,
                      weights= Worker.start, link = 'log', data=count)

prop.surv3 <- glmmTMB(worker.prop ~ Source.pop + (1 |Colony_ID) + (1 |Col_Season),
                      weights= Worker.start, beta_family(link = 'log'), data=count)

#FINAL MODEL 3: Colony size change. +  beginning
#colony beginning population model
beg.mod <- glmmTMB(Worker.start ~ Source.pop + Collection_date + (1 | ColonyID),
                   data = count, control = mod.control2, 
                   family = poisson(link = "log"))
#FINAL MODEL 4: Colony size end
end.mod <- glmmTMB(Worker.end ~ Source.pop + (1 | ColonyID), 
                  data = count, control = mod.control2, 
                  family = poisson(link = "log")) 

#FINAL MODEL 5: Change model (w/ both end and beg)
change.mod <- glmmTMB(Worker.pop ~ Time*Source.pop + (1 | Colony_ID) + (1 | Col_Season),
                      data = count1, family = poisson(link = "log"))


#plot basic modeled relationship to see if it matches expl graphics
library(visreg)
visreg(prop.surv)

## Part 5: Model diagnostics
#QQplots, residuals, AIC
##Diagnostics: warning msgs for models - show model convergence problem w/ small eigen value problems
diagnose(end.mod) #removed the interaction effect due to extremely large SD

diagnose(prop.surv1) 

library(DHARMa) #use simulated diagnostic modeling since we have a glmm; works with both glmer and glmmTMB objects
simOutput3 <- simulateResiduals(fittedModel = change.mod, plot = F)
plot(simOutput3)

testOutliers(simOutput3, margin = c("both"), type = c("bootstrap"), plot = T) #outliers are not signf. here


##Part 6: Statistical Testing
#gather summary data on model, determine model fit
summary(surv.mod)
summary(beg.mod) #pop start - Urban = 61.826, rural = 45.25
summary(end.mod) #pop end - Urban = 14.04, rural = 12.21
summary(prop.surv1)
summary(change.mod) #time = -1.31, urbanpop = 0.226, time*urbanpop = -0.172

# log lik statistic tests for model significance
library(car)
Anova(surv.mod, type = "III") #log likelihood
Anova(prop.surv, type = "III")
Anova(beg.mod, type = "III") #beginning #s: effect of source pop = 0.951, p = 0.33
#date effect = 0.57, p = 0.45, interaction effect = 0.954, p = 0.33
Anova(end.mod, type = "III") 
Anova(change.mod, type = "III") #source pop = chisq = 0.71, p = 0.399
#time = 396.03, p < 2e-16***, sourcepop*time = chisq = 3.64, p = 0.056

#backtransform log-linked data:
#using emmeans for figure and backtransformed values
#emmeans
library(emmeans)
surv.em <- emmeans(surv.mod, pairwise ~ Source.pop, transform = "response", adjust = "fdr")
surv.em
em_beg <- emmeans(beg.mod, pairwise ~ Source.pop, transform = "response", adjust = "fdr")
em_beg
#urban pop start = 49.3, SE = 6.16, rural start = 35.1, SE = 9.08
em_end <- emmeans(end.mod, pairwise ~ Source.pop, transform = "response", adjust = "fdr")
em_end #urban pop end = 9.6, SE = 2.09, rural pop end = 9.84, SE = 1.93
change.em <- emmeans(change.mod, pairwise ~ Time*Source.pop, transform = "response", adjust = "fdr")
change.em
##Part 7: Plotting Modeled Output
#
#Make follow up plot of modeled results
library(ggeffects)
prop.graph <- ggpredict(prop.surv1, terms = c("Source.pop"), 
                        type = "fe", allow.new.levels = TRUE,ci.lvl = 0.68, width=0.05, colour = my_colors3)
plot(prop.graph) #basic modeled plot


#make emmeans objects into dfs for visualization
#for survival overall
surv_df <- surv.em$emmeans %>%
  confint() %>%
  as.data.frame()
#make upper and lower SEs for Surv.em; manually calculated from emmeans
upper.SES <- c(29.56, 42.51)
lower.SES <- c(20.44, 29.23)

#for pop start
beg_df <- em_beg$emmeans %>%
  confint()%>%
  as.data.frame()
upper.SEb <- c(41.26, 58.38)
lower.SEb <- c(28.94, 40.22)

#for pop end
end_df <- em_end$emmeans %>%
  confint()%>%
  as.data.frame()
upper.SEd <- c(11.77, 11.69)
lower.SEd <- c(7.91, 7.51)

#for proportions - use summary data since not on Log scale

prop_em <- emmeans(prop.surv1, pairwise ~ Source.pop, transform = "log", adjust = "fdr")
upper.SEp <- c(0.738, 0.8321)
lower.SEp <- c(0.674, 0.7679)

prop_df <- prop_em$emmeans %>%
  confint()%>%
  as.data.frame()

#for change data
change_df <- change.em$emmeans %>%
  confint()%>%
  as.data.frame()
upper.SEg <- c(43.12, 11.7, 54.38, 12.41)
lower.SEg <- c(28.92, 7.74, 35.94, 8.11)

#plot graphs
#start graph
beg <- ggplot(beg_df, aes(x = Source.pop, y= rate))+  
  geom_point(colour = my_colors3, size = 7)+
  geom_errorbar(aes(ymin = lower.SEb, ymax = upper.SEb), fun.data = 'mean_se', 
                width=0.05, colour = my_colors3, fun.args = list(mult = 1))+
  labs(y="Worker Population Start", x = "Source Population")+
  ggtitle('Starting Population in Urban and Rural Colonies')+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )
#end graph
end <- ggplot(end_df, aes(x = Source.pop, y= rate))+  
  geom_point(colour = my_colors3, size = 7)+
  geom_errorbar(aes(ymin = lower.SEd, ymax = upper.SEd), fun.data = 'mean_se', 
                width=0.05, colour = my_colors3, fun.args = list(mult = 1))+
  labs(y="Worker Population End", x = "Source Population")+
  ggtitle('Ending Winter Population in Urban and Rural Colonies')+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )
library(ggplot2)
#change
change <- ggplot(change_df, aes(x = Time, y= rate))+  
  geom_point(color = my_colors4, size = 7)+
  geom_errorbar(aes(ymin = lower.SEg, ymax = upper.SEg), color = my_colors4, fun.data = 'mean_se', 
                width=0.05, fun.args = list(mult = 1))+
  labs(y="Worker Population", x = "Census Points", tag = "A")+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    plot.tag = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#total lost
total <- ggplot(surv_df, aes(x = Source.pop, y= rate))+  
  geom_point(colour = my_colors3, size = 7)+
  geom_errorbar(aes(ymin = lower.SES, ymax = upper.SES), fun.data = 'mean_se', 
                width=0.05, colour = my_colors3, fun.args = list(mult = 1))+
  labs(y="Worker Population Loss", x = "Source Population")+
  ggtitle('Winter Population Loss')+
  theme_classic()+
  theme(
    title = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10)
  )

#workers lost as a proportion of starting size
prop <- ggplot(prop_df, aes(x = Source.pop, y = emmean))+  
  geom_point(colour = my_colors3, size = 7)+
  geom_errorbar(aes(ymin = lower.SEp, ymax = upper.SEp), fun.data = 'mean_se', 
                width=0.05, colour = my_colors3, fun.args = list(mult = 1))+
  labs(y="Proportional Loss", x = "Source Population", tag = "B")+
  theme_classic()+
  theme(
    title = element_text(size = 14),
    plot.tag = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

#make graphs into combined facets
library(cowplot)
change
total
prop
plot_grid(change, prop, ncol = 2, nrow = 1, 
          rel_widths = c(.5, .5, .5))   #use 1028 x 555 px aspect ratio


my_colors3 <- c("cadetblue", "darkorange")
names(my_colors3) <- levels(count$Source.pop)

my_colors4 <- c("cadetblue", "cadet blue", "dark orange", "darkorange")
names(my_colors4) <- levels(c(count1$Source.pop, count1$Time))


#####################
################### END

