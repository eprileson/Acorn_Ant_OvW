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

#replace 1s in prop data with 0.999
count$worker.prop[count$worker.prop==1] <- 0.9999999
head(count)
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

#pivot the data longer so that we can show the change over time with the data set
library(magrittr)
library(tidyr)
count1 <- count1 %>%
  pivot_longer(c("Beginning", "End"), names_to = "Time", values_to = "Worker.pop")

count1$Time <- as.factor(count1$Time)
class(count1$Time)
head(count1)

#exploratory graphics: (worker pop ~ time*source.pop)
library(ggplot2)
ggplot(count, aes(x = Source.pop, y = worker.prop, colour = Source.pop))+
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', 
                width=0.05, fun.args = list(mult = 1))+
  geom_point(size = 7, stat = 'summary', fun.y = "mean")+
  scale_colour_manual(values = c("cadetblue", "darkorange"))+
  theme_classic()+
  labs(y = "Proportion of Worker Loss", x = "Source Population")+
  ggtitle('Proportion Worker Change')+
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
library(lme4)
library(MASS)
#first model for worker loss over 4 mo
#FINAL MODEL 1 - 
mod.control2 <- glmmTMBControl(optimizer=optim,
                               optArgs=list(method="BFGS"))
surv.mod <- glmmTMB(Worker.loss ~ Source.pop + (1 | ColonyID) + (1 | Collection_date:Source.pop), 
                    data = count, control = mod.control2, 
                    family = poisson(link = "log"))

#FINAL MODEL 2: Proportional Worker Loss
#adjust prop data = 1 to be 0.999
prop.surv <- glmmTMB(worker.prop ~ Source.pop  + (1 | Col_Season),
                    data=count, family = beta_family(link = "logit"),
                    weights = Worker.start)

#model with binomial dist to try
prop.surv2 <- glmmTMB(worker.prop ~ Source.pop + (1 |Col_Season),
                      weights= Worker.start, data=count, family = binomial(link = "logit"))

#model with quasi-binomial using glm
prop.surv3 <- glm(worker.prop ~ Source.pop*Col_Season,
                  data=count, family = quasibinomial(link = "logit"),
                  weights = Worker.start)

#FINAL MODEL: Change model (w/ both end and beg)
change.mod <- glmmTMB(Worker.pop ~ Time*Source.pop + (1 | Col_Season),
                      data = count1, family = poisson(link = "log"))


## Part 5: Model diagnostics
#QQplots, residuals, AIC
##Diagnostics: warning msgs for models - show model convergence problem w/ small eigen value problems
diagnose(change.mod) #removed the interaction effect due to extremely large SD

diagnose(prop.surv3) 

library(DHARMa) #use simulated diagnostic modeling since we have a glmm; works with both glmer and glmmTMB objects
simOutput3 <- simulateResiduals(fittedModel = prop.surv3)
plot(simOutput3)
plotResiduals(simOutput3)
testOutliers(simOutput3, margin = c("both"), type = c("bootstrap"), plot = T) #outliers are not signf. here

#regular diagnostic tests for glm
plot(prop.surv3)

#plot basic modeled relationship to see if it matches expl graphics
library(visreg)
visreg(prop.surv)


##Part 6: Statistical Testing
#gather summary data on model, determine model fit
summary(surv.mod)
summary(prop.surv3)
summary(change.mod) #time = -1.31, urbanpop = 4.177, time*urbanpop = -0.172

# log lik statistic tests for model significance
library(car)
Anova(surv.mod, type = "III") #log likelihood
Anova(prop.surv3, type = "III")
Anova(change.mod, type = "III") #source pop = chisq = 0.71, p = 0.399
#time = 396.03, p < 2e-16***, sourcepop*time = chisq = 3.64, p = 0.056

#backtransform log-linked data:
#using emmeans for figure and backtransformed values
#emmeans
library(emmeans)
surv.em <- emmeans(surv.mod, pairwise ~ Source.pop, transform = "response", adjust = "fdr")
surv.em
change.em <- emmeans(change.mod, pairwise ~ Time*Source.pop, transform = "response", adjust = "fdr")
change.em
##Part 7: Plotting Modeled Output
#
#Make follow up plot of modeled results
library(ggeffects)
prop.graph <- ggpredict(prop.surv3, terms = c("Source.pop"), 
                        type = "fe",ci.lvl = 0.68, width=0.05, colour = my_colors3)
plot(prop.graph) #basic modeled plot


#make emmeans objects into dfs for visualization
#for survival overall
surv_df <- surv.em$emmeans %>%
  confint() %>%
  as.data.frame()
#make upper and lower SEs for Surv.em; manually calculated from emmeans
upper.SES <- c(29.56, 42.51)
lower.SES <- c(20.44, 29.23)

#for proportions - use summary data since not on Log scale

prop_em <- emmeans(prop.surv, pairwise ~ Source.pop, transform = "response", adjust = "fdr")
upper.SEp <- c(0.7958, 0.8057)
lower.SEp <- c(0.732, 0.744)

prop_df <- prop_em$emmeans %>%
  confint()%>%
  as.data.frame()

#for change data
change_df <- change.em$emmeans %>%
  confint()%>%
  as.data.frame()
upper.SEg <- c(48.56, 13.22, 71.43, 16.41)
lower.SEg <- c(39.84, 10.58, 58.77, 13.19)

#plot graphs
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
prop <- ggplot(prop_df, aes(x = Source.pop, y = response))+  
  geom_point(colour = my_colors3, size = 7)+
  geom_errorbar(aes(ymin = lower.SEp, ymax = upper.SEp), fun.data = 'mean_se', 
                width=0.05, colour = my_colors3, fun.args = list(mult = 1))+
  labs(y="Proportional Worker Loss", x = "Source Population", tag = "B")+
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

