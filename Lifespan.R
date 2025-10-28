#---------------------Chapter 2.3 GLMM Lifespan analysis-----------------
#Lets doing a GLMM on lifespan.

#--------------------------------libraries-------------
if (!is_installed("coxme")) { 
  install.packages(coxme)}

if (!is_installed("arm")) { 
  install.packages(arm)}

if (!is_installed("survival")) { 
  install.packages(survival)}

if (!is_installed("glmmTMB")) { 
  install.packages(glmmTMB)}

if (!is_installed("DHARMa")) { 
  install.packages(DHARMa)}

if (!is_installed("beepr")) { 
  install.packages(beepr)}

if (!is_installed("performance")) { 
  install.packages(performance)}

if (!is_installed("survminer")) { 
  install.packages(survminer)}

if (!is_installed("data.table")) { 
  install.packages(data.table)}

if (!is_installed("car")) { 
  install.packages(car)}

if (!is_installed("MuMIn")) { 
  install.packages("MuMIn")}

if (!is_installed("lme4")) { 
  install.packages("lme4")}

if (!is_installed("patchwork")) { 
  install.packages("patchwork")}

install.packages("ggsignif")
install.packages("performance")
install.packages("see")
install.packages("stargazer")
install.packages("viridis")

library(ggplot2)
library(ggeffects)
library("ggsignif")   
library(MuMIn)
library(car)
library(coxme)
library(arm)
library(survival)
library(glmmTMB)
library(DHARMa)
library(beepr)
library(performance)
library(data.table)
library(survminer)
library(dplyr)
library(stats)
library(lme4)
library(patchwork)
library(viridis)
library(modelsummary)
library(gtsummary)
library(ggiraphExtra)
library(ggiraph)
library(interactions)

#----------------------------------import data--------------

##Factorise
stat.birds <- read.csv("Data/clean_data.csv")
stat.birds$h_countF <- as.factor(stat.birds$h_count)
stat.birds$BreedSeason <- as.factor(stat.birds$BreedSeason)
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$DeathEventL <- as.logical(stat.birds$DeathEvent)
stat.birds$mum <- as.factor(stat.birds$mum)
stat.birds$dad <- as.factor(stat.birds$dad)
stat.birds$birth_yearF <- as.factor(stat.birds$birth_year)
stat.birds$adulthood <- as.factor(stat.birds$adulthood)
stat.birds$sib_pres <- as.factor(stat.birds$sib_pres)


# select males and females, filter for coverage

just_birds <- stat.birds %>% filter(event == 0,
                                    Coverage >= 0.25,
                                    # remove some NAs
                                    !is.na(mean_total_rain),
                                    !is.na(mean_rain_cv),
                                    !is.na(mum),
                                    !is.na(dad)
)





female_birds <- stat.birds %>% filter(Sex == 0,
                                      event == 0,
                                      Coverage >= 0.25 ,
                                      CovMum >= 0.25 ,
                                      CovDad >= 0.25,
                                      CovBrM >= 0.25,
                                      !is.na(mean_total_rain),
                                      !is.na(mean_rain_cv),
                                      !is.na(mum),
                                      !is.na(dad))

male_birds <- stat.birds %>% filter(Sex == 1, 
                                    event == 0,
                                    Coverage >= 0.25 ,
                                    CovMum >= 0.25 ,
                                    CovDad >= 0.25,
                                    CovBrM >= 0.25,
                                    !is.na(mean_total_rain),
                                    !is.na(mean_rain_cv),
                                    !is.na(mum),
                                    !is.na(dad))


both_birds <- stat.birds %>% filter(event == 0,
                                    Coverage >= 0.25 ,
                                    CovMum >= 0.25 ,
                                    CovDad >= 0.25,
                                    CovBrM >= 0.25,
                                    !is.na(mean_total_rain),
                                    !is.na(mean_rain_cv),
                                    !is.na(mum),
                                    !is.na(dad))


dim(both_birds) #n = 499


#begin testing for families

both.lifespan.poisson <- glmmTMB(lifespan ~ rescale(LargeFROH) + rescale(BrMLargeFroh)
                                 +rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                   Sex + h_countF + 
                                   rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                   sib_pres +
                                   natal_group_size +
                                   EPP +
                                   rescale(birth_year)+
                                   (1 | mum) + (1 | dad) ,
                                 data = both_birds,
                                 ziformula=~0,
                                 family=poisson
)


simulateResiduals(both.lifespan.poisson, plot = T) # looks a bit wobbly

check_collinearity(both.lifespan.poisson) #vifs around 1

check_overdispersion(both.lifespan.poisson) #none detected

summary(both.lifespan.poisson) #converges ok


#nbinom1

both.lifespan.nbinom1 <- update(both.lifespan.poisson, family = nbinom1)

simulateResiduals(both.lifespan.nbinom1, plot = T) # looks slightly better

check_collinearity(both.lifespan.nbinom1) #vifs around 1

check_overdispersion(both.lifespan.nbinom1) #none detected

summary(both.lifespan.nbinom1) #converges ok

#nbinom2

both.lifespan.nbinom2 <- update(both.lifespan.poisson, family = nbinom2) #step failure

simulateResiduals(both.lifespan.nbinom2, plot = T) # outliers detected but line looks ok

check_collinearity(both.lifespan.nbinom2) #vifs around 1

check_overdispersion(both.lifespan.nbinom2) #none detected

summary(both.lifespan.nbinom2) #converges ok

#check AIC

AIC(both.lifespan.poisson, both.lifespan.nbinom1, both.lifespan.nbinom2)

#both.lifespan.poisson 17 2573.122
#both.lifespan.nbinom1 18 2487.459
#both.lifespan.nbinom2 18 2452.439

#nbinom 2 is lowest AIC 

#frohs w/ sex

#brm * sex

both.lifespan.nbinom2.a <- update(both.lifespan.nbinom2, ~ . + rescale(BrMLargeFroh)*Sex) 

simulateResiduals(both.lifespan.nbinom2.a, plot = T) # looks ok

check_collinearity(both.lifespan.nbinom2.a) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.a) #none detected

summary(both.lifespan.nbinom2.a) #ns

#dad * sex

both.lifespan.nbinom2.b <- update(both.lifespan.nbinom2, ~ . + rescale(dadLargeFroh)*Sex) 

simulateResiduals(both.lifespan.nbinom2.b, plot = T) # looks a bit weird

check_collinearity(both.lifespan.nbinom2.b) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.b) #none detected

summary(both.lifespan.nbinom2.b) #it's almost significant!

anova(both.lifespan.nbinom2,both.lifespan.nbinom2.b) #anova says almost?

#mum * sex

both.lifespan.nbinom2.c <- update(both.lifespan.nbinom2, ~ . + rescale(mumLargeFroh)*Sex) 

simulateResiduals(both.lifespan.nbinom2.c, plot = T) # looks a bit weird

check_collinearity(both.lifespan.nbinom2.c) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.c) #none detected

summary(both.lifespan.nbinom2.c) #ns

#dadFROH * sex interaction!

#what happens in this interaction???

#quick plot to see

ggplot(just_birds, aes(dadLargeFroh,lifespan,colour = Sex))+
  geom_smooth(method = "lm")

#from this know that males react differently from females. But are they different from zero?
#split model by sex to find out

females.lifespan.nbinom1 <- glmmTMB(lifespan ~ rescale(LargeFROH) + rescale(BrMLargeFroh)
                                 +rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                    h_countF + 
                                   rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                   rescale(birth_year) + 
                                   (1 | mum) + (1 | dad) ,
                                 data = female_birds,
                                 ziformula=~0,
                                 family=nbinom1
)

summary(females.lifespan.nbinom1) #no effect of dadFROH

males.lifespan.nbinom1 <- glmmTMB(lifespan ~ rescale(LargeFROH) + rescale(BrMLargeFroh)
                                    +rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                      h_countF + 
                                      rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                      rescale(birth_year) + 
                                      (1 | mum) + (1 | dad) ,
                                    data = male_birds,
                                    ziformula=~0,
                                    family=nbinom1
)


summary(males.lifespan.nbinom1) #a marginal, but not quite significant POSITIVE effect of genetic dad and mum

#quick plot to look at this

ggplot(both_birds, aes(x= dadLargeFroh,y= lifespan))+
  geom_smooth(method = "lm")


#try help

#brm * help

both.lifespan.nbinom2.d<- update(both.lifespan.nbinom2, ~ . + rescale(BrMLargeFroh)*h_countF) 

simulateResiduals(both.lifespan.nbinom2.d, plot = T) # looks ok

check_collinearity(both.lifespan.nbinom2.d) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.d) #none detected

summary(both.lifespan.nbinom2.d) #not quite ns

#dad * help

both.lifespan.nbinom2.e <- update(both.lifespan.nbinom2, ~ . + rescale(dadLargeFroh)*h_countF) 

simulateResiduals(both.lifespan.nbinom2.e, plot = T) # looks a bit weird

check_collinearity(both.lifespan.nbinom2.e) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.e) #none detected

summary(both.lifespan.nbinom2.e) #ns

anova(both.lifespan.nbinom2, both.lifespan.nbinom2.e) #ns

#mum * help

both.lifespan.nbinom2.f <- update(both.lifespan.nbinom2, ~ . + rescale(mumLargeFroh)*h_countF) 

simulateResiduals(both.lifespan.nbinom2.f, plot = T) # looks a bit weird

check_collinearity(both.lifespan.nbinom2.f) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.f) #none detected

summary(both.lifespan.nbinom2.f) #ns


#interactions with total rain

#brm * rescale(mean_total_rain

both.lifespan.nbinom2.g<- update(both.lifespan.nbinom2, ~ . + rescale(BrMLargeFroh)*rescale(mean_total_rain) )

simulateResiduals(both.lifespan.nbinom2.g, plot = T) # looks ok

check_collinearity(both.lifespan.nbinom2.g) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.g) #none detected

summary(both.lifespan.nbinom2.g) #ns

#dad * rescale(mean_total_rain

both.lifespan.nbinom2.h <- update(both.lifespan.nbinom2, ~ . + rescale(dadLargeFroh)*rescale(mean_total_rain)) 

simulateResiduals(both.lifespan.nbinom2.h, plot = T) # looks a bit weird

check_collinearity(both.lifespan.nbinom2.h) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.h) #none detected

summary(both.lifespan.nbinom2.h) #ns

#mum * rescale(mean_total_rain

both.lifespan.nbinom2.i <- update(both.lifespan.nbinom2, ~ . + rescale(mumLargeFroh)*rescale(mean_total_rain) )

simulateResiduals(both.lifespan.nbinom2.i, plot = T) # looks a bit weird

check_collinearity(both.lifespan.nbinom2.i) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.i) #none detected

summary(both.lifespan.nbinom2.i) #ns


#interactions with rain variance

#brm * rescale(mean_rain_cv)

both.lifespan.nbinom2.j<- update(both.lifespan.nbinom2, ~ . + rescale(BrMLargeFroh)*rescale(mean_rain_cv) )

simulateResiduals(both.lifespan.nbinom2.j, plot = T) # looks ok

check_collinearity(both.lifespan.nbinom2.j) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.j) #none detected

summary(both.lifespan.nbinom2.j) #ns

#dad * rescale(mean_rain_cv

both.lifespan.nbinom2.k <- update(both.lifespan.nbinom2, ~ . + rescale(dadLargeFroh)*rescale(mean_rain_cv)) 

simulateResiduals(both.lifespan.nbinom2.k, plot = T) # looks a bit weird

check_collinearity(both.lifespan.nbinom2.k) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.k) #none detected

summary(both.lifespan.nbinom2.k) #ns

#mum * rescale mean_rain_cv

both.lifespan.nbinom2.l <- update(both.lifespan.nbinom2, ~ . + rescale(mumLargeFroh)*rescale(mean_rain_cv) )

simulateResiduals(both.lifespan.nbinom2.l, plot = T) # looks a bit weird

check_collinearity(both.lifespan.nbinom2.l) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.l) #none detected

summary(both.lifespan.nbinom2.l) #not

anova(both.lifespan.nbinom2, both.lifespan.nbinom2.l) #

#so there is no significant interactions. dad x sex is marginal.

summary(both.lifespan.nbinom2.b) #dadFROH * sex


#but base model is still best

summary(both.lifespan.nbinom2)

#output######

#lets plot the base model first

summary(both.lifespan.nbinom2)

tbl_regression(both.lifespan.nbinom2, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH", 
                            "h_countF" = "Helper in natal territory",
                            "rescale(mean_total_rain)"	 = "Mean annual rainfall during lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance annual rainfall during lifespan",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH",
                            "rescale(birth_year)" = "Hatch Year",
                            "sib_pres" = "Sibling presence in nest",
                            "natal_group_size" = "Natal group size"),
               estimate_fun = label_style_sigfig(digits = 3),
               pvalue_fun = label_style_pvalue(digits = 3)
)%>% 
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) %>%
  
  modify_fmt_fun(
    std.error  = function(x) style_sigfig(x, digits = 3),
    statistic  = function(x) style_sigfig(x, digits = 3)
  )

#plot base model


#genetic male Froh
plot.both.lifespan.gendad <-ggplot(both_birds, aes(dadLargeFroh, lifespan)) +
  geom_smooth(method = "lm")+
  geom_point()+
  theme_classic2()+
  xlim(0,0.5)+
  labs( x = bquote("Genetic Father's F"[ROH]),
        y= " ")+
  theme(text = element_text(size = 20))+
  ggtitle("Fig. 5.1.
a")



#social dad Froh
plot.both.lifespan.brm <-ggplot(both_birds, aes(BrMLargeFroh, lifespan)) +
  geom_smooth(method = "lm")+
  geom_point()+
  theme_classic2()+
  labs( x = bquote("Social Father's F"[ROH]),
        y= "Lifespan")+
  theme(text = element_text(size = 20))+
  ggtitle("b")



#mumFroh
  plot.both.lifespan.mum<-ggplot(both_birds, aes(mumLargeFroh, lifespan)) +
  geom_smooth(method = "lm")+
  geom_point()+
    xlim(0,0.5)+
  theme_classic2()+
  labs( x = bquote("Mother's F"[ROH]),
        y = "")+
  theme(text = element_text(size = 20))+
  ggtitle("c")


(plot.both.lifespan.gendad / plot.both.lifespan.brm / plot.both.lifespan.mum)


#output model with both sig interactions

summary(both.lifespan.nbinom2.b)

tbl_regression(both.lifespan.nbinom2.b, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH",
                            "rescale(mean_total_rain)"	 = "Mean annual rainfall during lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance annual rainfall during lifespan",
                            "rescale(birth_year)" = "Birth Year",
                            "rescale(dadLargeFroh)*Sex" = "Genetic father's FROH * Sex"),
               estimate_fun = label_style_sigfig(digits = 3),
               pvalue_fun = label_style_pvalue(digits = 3)
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  %>%
  
  modify_fmt_fun(
    std.error  = function(x) style_sigfig(x, digits = 3),
    statistic  = function(x) style_sigfig(x, digits = 3)
  )

#plot dad * sex


ggplot(both_birds,aes(dadLargeFroh, lifespan, group = Sex, fill = Sex))+
  geom_smooth(method = "lm")+
  geom_point(aes(colour = Sex))+
  theme_classic2()+
  labs( x = bquote("Genetic Father's F"[ROH]),
        y= "Lifespan",)+
  theme(text = element_text(size = 20)) +
  ggtitle("Fig. 5.2. ")



#output ables for sex split models


#for females
tbl_regression(females.lifespan.nbinom1 , intercept = T,
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH",
                            "rescale(mean_total_rain)"	 = "Mean annual rainfall during lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance annual rainfall during lifespan",
                            "rescale(birth_year)" = "Birth Year",
                            "rescale(dadLargeFroh)*Sex" = "Genetic father's FROH * Sex"
                            
               )
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  


#for males
tbl_regression(males.lifespan.nbinom1 , intercept = T,
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH",
                            "rescale(mean_total_rain)"	 = "Mean annual rainfall during lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance annual rainfall during lifespan",
                            "rescale(birth_year)" = "Birth Year",
                            "rescale(dadLargeFroh)*Sex" = "Genetic father's FROH * Sex"
                            
               )
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  


















#trying to use interaction plot
#try un rescaling?

int.both.lifespan.nbinom2.l <-  glmmTMB(lifespan ~ LargeFROH + BrMLargeFroh
                                                                +dadLargeFroh + mumLargeFroh +
                                                                  Sex + help + 
                                                                  mean_total_rain + mean_rain_cv +
                                          mumLargeFroh * mean_rain_cv +
                                                                  birth_year + 
                                                                  (1 | mum) + (1 | dad) ,
                                                                data = both_birds,
                                                                ziformula=~1,
                                                                family=poisson
)



plot_birds <- na.omit(both_birds)

#

interact_plot(int.both.lifespan.nbinom2.l,
              pred = mumLargeFroh,
              modx = mean_rain_cv,
              modx.values = "terciles",
              int.type = "confidence",
              plot.points = T,
              data = plot_birds)+
              theme_classic2()+
              labs( x = "Mother FROH",
              y= "Lifespan",
              fill = "Mean annual variance in rainfall over lifespan")+
              theme(text = element_text(size = 20)) +
              ggtitle("Fig. 6.3 - Predicted fit of GLMM lifespan ~ Mother FROH
              with Mean Annual Variance in rainfall interaction")



#------------plot males only subset-------------


























##-------------just inbreeding---------

just_birds <- stat.birds %>% filter(event == 0, Coverage >= 1 )
dim(just_birds)

just.lifespan.1 <- glmmTMB(lifespan ~ FROH  + Sex + help + MumAge + meanInsects +
                              (1 | birth_year) + (1 | mum) + (1 | dad) ,
                            data = just_birds,
                            ziformula=~1,
                            family=poisson
)


summary(just.lifespan.1)
simulateResiduals(just.lifespan.1, plot = TRUE)


hist(early_birds$lifespan)

ggplot(data = just_birds, aes(x = FROH, y = lifespan))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

##----------just inbreeding, filtered for first year sampling----------------




#lifespan model

early_birds <- stat.birds %>% filter(event == 0, Coverage >= 1 , SamplingAge == 0 )

dim(early_birds)

early.lifespan.1 <- glmmTMB(lifespan ~ rescale(FROH) + rescale(meanInsects) + Sex + 
                              (1 | birth_year) + (1 | mum) + (1 | dad) ,
                          data = early_birds,
                          ziformula=~1,
                          family=poisson
                          )

summary(early.lifespan.1)

hist(early_birds$lifespan)

ggplot(data = early_birds, aes(x = FROH, y = lifespan))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()


vif(early.lifespan.1)
simulateResiduals(early.lifespan.1, plot = TRUE)



#first year survival success

early.recruitment.1 <- glmer(adulthood ~ FROH + meanInsects + MumAge + FROH*meanInsects + (1 | birth_year) + (1 | mum) + (1 | dad),
                             data = early_birds,
                             family=binomial,
                             control=glmerControl())

summary(early.recruitment.1)

vif(early.recruitment.1)
simulateResiduals(early.recruitment.1, plot = TRUE)


ggplot(data = early_birds, aes(x = FROH, y = adulthood))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

plot(ggpredict(early.recruitment.1, terms="FROH [all]"), show_data = TRUE)



#lifetime reproductive success

#find out what makes FROH sig or not

early.lrs.1 <- glmmTMB(n_off ~ FROH  + lifespan ,
                     data = early_birds,
                     ziformula = ~0,
                     family = poisson())

summary(early.lrs.1)
vif(early.lrs.1)
simulateResiduals(early.lrs.1, plot = TRUE)
#FROH  -4.15345    1.15648  -3.591 0.000329 *** when only including FRHO and lifespan



ggplot(data = early_birds, aes(x = n_off, y = FROH))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  labs(x = "Lifetime Reproductive Success")

#add birth year

early.lrs.2 <- glmmTMB(n_off ~ FROH  + lifespan + (1 | birth_year),
                       data = early_birds,
                       ziformula = ~0,
                       family = poisson())

summary(early.lrs.2)#  FROH        -2.17750    1.23338  -1.765 0.077484 . not quite sig when birth year included


#add mum id 
early.lrs.3 <- glmmTMB(n_off ~ FROH  + lifespan + (1 | mum),
                       data = early_birds,
                       ziformula = ~0,
                       family = poisson())

summary(early.lrs.3)# FROH        -0.55333    2.00049  -0.277    0.782 not sig at all when mum included. n = 314


#-----------------select birds first sampled under 6 months-----


#filter birds 
earlier_birds <- stat.birds %>% filter(event == 0, Coverage >= 1 , first_catch > 0 )
dim(earlier_birds)

#lifespan model

earlier.lifespan.1 <- glmmTMB(lifespan ~ FROH + meanInsects + Sex + helpF + MumAge +
                              (1 | birth_year) + (1 | mum) + (1 | dad) ,
                            data = earlier_birds,
                            ziformula=~1,
                            family=poisson
)

summary(earlier.lifespan.1) #FROH        -1.178263   0.621354  -1.896   0.0579 .  

hist(early_birds$lifespan)

ggplot(data = early_birds, aes(x = FROH, y = lifespan))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()


vif(early.lifespan.1)
simulateResiduals(early.lifespan.1, plot = TRUE)



#first year survival success

earlier.recruitment.1 <- glmer(adulthood ~ FROH + meanInsects + MumAge +  
                               (1 | birth_year) + (1 | mum) + (1 | dad),
                             data = earlier_birds,
                             family=binomial,
                             control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(earlier.recruitment.1)

vif(earlier.recruitment.1)
simulateResiduals(earlier.recruitment.1, plot = TRUE)


ggplot(data = early_birds, aes(x = FROH, y = adulthood))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

plot(ggpredict(early.recruitment.1, terms="FROH [all]"), show_data = TRUE)



#lifetime reproductive success

#find out what makes FROH sig or not

earlier.lrs.1 <- glmmTMB(n_off ~ FROH  + lifespan + meanInsects + MumAge +  
                           (1 | birth_year) + (1 | mum) + (1 | dad),
                       data = earlier_birds,
                       ziformula = ~1,
                       family = poisson())

summary(earlier.lrs.1) #FROH        -2.396378   0.813840  -2.945  0.00323 **
vif(earlier.lrs.1)
simulateResiduals(earlier.lrs.1, plot = TRUE)
hist(earlier_birds$n_off)

ggplot(data = early_birds, aes(x = n_off, y = FROH))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  labs(x = "Lifetime Reproductive Success")





#-----------------select birds first sampled under 1 months-----


#filter birds 
earliest_birds <- stat.birds %>% filter(event == 0, Coverage >= 1 , first_catch == 1 )
dim(earliest_birds)



#lifespan model

earliest.lifespan.1 <- glmmTMB(lifespan ~ FROH + meanInsects + Sex + helpF + MumAge +
                                (1 | birth_year) + (1 | mum) + (1 | dad) ,
                              data = earliest_birds,
                              ziformula=~1,
                              family=poisson
)

summary(earliest.lifespan.1) #FROH        -1.081535   0.650367  -1.663   0.0963 . 

hist(early_birds$lifespan)

ggplot(data = early_birds, aes(x = FROH, y = lifespan))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()


vif(early.lifespan.1)
simulateResiduals(early.lifespan.1, plot = TRUE)



#first year survival success

earlier.recruitment.1 <- glmer(adulthood ~ FROH + meanInsects + MumAge +  
                                 (1 | birth_year) + (1 | mum) + (1 | dad),
                               data = earlier_birds,
                               family=binomial,
                               control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(earlier.recruitment.1)

vif(earlier.recruitment.1)
simulateResiduals(earlier.recruitment.1, plot = TRUE)


ggplot(data = early_birds, aes(x = FROH, y = adulthood))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

plot(ggpredict(early.recruitment.1, terms="FROH [all]"), show_data = TRUE)



#lifetime reproductive success

#find out what makes FROH sig or not

earlier.lrs.1 <- glmmTMB(n_off ~ FROH  + lifespan + meanInsects + MumAge +  
                           (1 | birth_year) + (1 | mum) + (1 | dad),
                         data = earlier_birds,
                         ziformula = ~1,
                         family = poisson())

summary(earlier.lrs.1) #FROH        -2.396378   0.813840  -2.945  0.00323 **
vif(earlier.lrs.1)
simulateResiduals(earlier.lrs.1, plot = TRUE)
hist(earlier_birds$n_off)

ggplot(data = early_birds, aes(x = n_off, y = FROH))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  labs(x = "Lifetime Reproductive Success")






##--------------------------joint sex models----------------

#filter out translocated/alive birds

both_birds <- stat.birds %>% filter(event == 0, Coverage >= 1 , CovMum >= 1 , CovDad >= 1, CovBrM >= 1)

dim(both_birds) #534 samples

#check distribution
hist(both_birds$lifespan) #it's pretty zero inflated

#so we need to use a zim

both.lifespan.1 <- glmer(lifespan ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                          data = both_birds,
                          family= gaussian(link = "identity"),
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

summary(both.lifespan.1)  #nothing significant when dataset is not split by EPP

vif(both.lifespan.1)
simulateResiduals(both.lifespan.1, plot = TRUE)




#split birds into EPP and not EPP 

cuck_both_birds <- both_birds %>% filter(EPP == 1)
dim(cuck_both_birds) #222 birds from EPP

non_cuck_both_birds <- both_birds %>% filter(EPP == 0)
dim(non_cuck_both_birds) #312 birds not from EPP


#not EPP tests.
#We need to test if sex has interaction with mumF and dadF

#basic model, no interactions
non_cuck_both_lifespan_base <- glmer(lifespan ~ mumF + dadF + FROH + help + MumAge + Sex + (1 | birth_year) + (1 | mum) + (1 | dad),
                                      data = non_cuck_both_birds,
                                      family=binomial,
                                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

summary(non_cuck_both_lifespan_base) #mumage and sex are significant


#-------------------------tests for female birds------------

female.lifespan.1 <- glmer.nb(lifespan ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                            data = female_birds                    )

summary(female.lifespan.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = female.lifespan.1, plot = T)

#check for collinearirty with vif
vif(female.lifespan.1)

##plot basic model with no interactions
#plot recruitment to lifespan for different Fs

female.lifespan.indF.plot <- ggplot(female_birds, aes(lifespan, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.lifespan.mumF.plot <- ggplot(female_birds, aes(lifespan,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.lifespan.dadF.plot <- ggplot(female_birds, aes(lifespan,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Lifespan", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.lifespan.indF.plot <- ggplot(male_birds, aes(lifespan, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.lifespan.mumF.plot <- ggplot(male_birds, aes(lifespan,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.lifespan.dadF.plot <- ggplot(male_birds, aes(lifespan,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Lifespan", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


lifespan.global <- (female.lifespan.indF.plot + male.lifespan.indF.plot) / (female.lifespan.mumF.plot + male.lifespan.mumF.plot) / (female.lifespan.dadF.plot + male.lifespan.dadF.plot)

ggsave("Figures/lifespan_global.png", lifespan.global)

#----------------female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
female.lifespan.2 <- update(female.lifespan.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.lifespan.2) #interaction term not significant.
vif(female.lifespan.2)
anova(female.lifespan.1,female.lifespan.2) #models not significantly different


#check dadF help interaction. 
female.lifespan.3 <- update(female.lifespan.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.lifespan.3) #interaction term very almost significant p = 0.0503 .
vif(female.lifespan.3)
anova(female.lifespan.1,female.lifespan.3)  #models are not significantly different



#check FROH help interaction. 
female.lifespan.4 <- update(female.lifespan.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.lifespan.4) #interaction term not sig.
vif(female.lifespan.3)
anova(female.lifespan.1,female.lifespan.3)  #models are significantly different

#----------------------------------female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
female.lifespan.4 <- update(female.lifespan.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.lifespan.4) #interaction term not significant.
vif(female.lifespan.2)
anova(female.lifespan.1,female.lifespan.2) #models not significantly different


#check dadF MumAge interaction. 
female.lifespan.5 <- update(female.lifespan.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.lifespan.5) #interaction term not sig .
vif(female.lifespan.3)
anova(female.lifespan.1,female.lifespan.5)  #models are not significantly different


#check FROH MumAge interaction. 
female.lifespan.6 <- update(female.lifespan.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.lifespan.6) #interaction term not sig.
vif(female.lifespan.3)
anova(female.lifespan.1,female.lifespan.6)  #models are not significantly different

#-------------------------tests for male birds------------
#basic model
male.lifespan.1 <- glmer(lifespan ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                          data = male_birds,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa"))

summary(male.lifespan.1) #here FROH and mumage arte sig

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = male.lifespan.1, plot = T)

#check for collinearirty with vif
vif(male.lifespan.1)

#-----------------male helper interactions-----------

#check mumF help interaction. 
male.lifespan.2 <- update(male.lifespan.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.lifespan.2) #interaction term not significant.
vif(male.lifespan.2)
anova(male.lifespan.1,male.lifespan.2) #models not significantly different


#check dadF help interaction. 
male.lifespan.3 <- update(male.lifespan.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.lifespan.3) #interaction term not significant
vif(male.lifespan.3)
anova(male.lifespan.1,male.lifespan.3)  #models are not significantly different


#check FROH help interaction. 
male.lifespan.4 <- update(male.lifespan.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.lifespan.4) #interaction term not sig. Nearly though
vif(male.lifespan.3)
anova(male.lifespan.1,male.lifespan.3)  #models are not significantly different


#-------------------male mumage interaction-------

#check mumF MumAge interaction. 
male.lifespan.4 <- update(male.lifespan.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.lifespan.4) #interaction term not significant.
vif(male.lifespan.2)
anova(male.lifespan.1,male.lifespan.2) #models not significantly different


#check dadF MumAge interaction. 
male.lifespan.5 <- update(male.lifespan.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.lifespan.5) #interaction term not sig .
vif(male.lifespan.3)
anova(male.lifespan.1,male.lifespan.5)  #models are not significantly different


#check FROH MumAge interaction. 
male.lifespan.6 <- update(male.lifespan.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.lifespan.6) #interaction term not sig.
vif(male.lifespan.3)
anova(male.lifespan.1,male.lifespan.6)  #models are not significantly different


##--------------------------------EPP tests-------------------------------------


#select non cuck offspring
non_cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for non_cucks_female  birds------------
non_cucks_female.lifespan.1 <- glmer(lifespan ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                      data = non_cucks_female_birds,
                                      family=binomial,
                                      control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(non_cucks_female.lifespan.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_female.lifespan.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_female.lifespan.1)

##plot basic model with no interactions
#plot recruitment to lifespan for different Fs

non_cucks_female.lifespan.indF.plot <- ggplot(non_cucks_female_birds, aes(lifespan, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.lifespan.mumF.plot <- ggplot(non_cucks_female_birds, aes(lifespan,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.lifespan.dadF.plot <- ggplot(non_cucks_female_birds, aes(lifespan,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Lifespan", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.lifespan.indF.plot <- ggplot(male_birds, aes(lifespan, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.lifespan.mumF.plot <- ggplot(male_birds, aes(lifespan,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.lifespan.dadF.plot <- ggplot(male_birds, aes(lifespan,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Lifespan", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


lifespan.global <- (non_cucks_female.lifespan.indF.plot + male.lifespan.indF.plot) / (non_cucks_female.lifespan.mumF.plot + male.lifespan.mumF.plot) / (non_cucks_female.lifespan.dadF.plot + male.lifespan.dadF.plot)

ggsave("Figures/lifespan_global.png", lifespan.global)

#----------------non_cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
non_cucks_female.lifespan.2 <- update(non_cucks_female.lifespan.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                         optCtrl = list(maxfun = 100000000)))

summary(non_cucks_female.lifespan.2) #interaction term not significant.
vif(non_cucks_female.lifespan.2)
anova(non_cucks_female.lifespan.1,non_cucks_female.lifespan.2) #models not significantly different


#check dadF help interaction. 
non_cucks_female.lifespan.3 <- update(non_cucks_female.lifespan.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.lifespan.3) #interaction term very almost significant p = 0.0503 .
vif(non_cucks_female.lifespan.3)
anova(non_cucks_female.lifespan.1,non_cucks_female.lifespan.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_female.lifespan.4 <- update(non_cucks_female.lifespan.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.lifespan.4) #interaction term not sig.
vif(non_cucks_female.lifespan.3)
anova(non_cucks_female.lifespan.1,non_cucks_female.lifespan.3)  #models are significantly different

#----------------------------------non_cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
non_cucks_female.lifespan.4 <- update(non_cucks_female.lifespan.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.lifespan.4) #interaction term not significant.
vif(non_cucks_female.lifespan.2)
anova(non_cucks_female.lifespan.1,non_cucks_female.lifespan.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_female.lifespan.5 <- update(non_cucks_female.lifespan.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.lifespan.5) #interaction term not sig .
vif(non_cucks_female.lifespan.3)
anova(non_cucks_female.lifespan.1,non_cucks_female.lifespan.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_female.lifespan.6 <- update(non_cucks_female.lifespan.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.lifespan.6) #interaction term not sig.
vif(non_cucks_female.lifespan.3)
anova(non_cucks_female.lifespan.1,non_cucks_female.lifespan.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select non bastard
non_cucks_male_birds <-filter(male_birds, EPP == 0) 

#basic model
non_cucks_male.lifespan.1 <- glmer(lifespan ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                                    data = non_cucks_male_birds,
                                    family=binomial,
                                    control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.lifespan.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_male.lifespan.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_male.lifespan.1)


#-----------------non_cucks_male helper interactions-----------





#check mumF help interaction. 
non_cucks_male.lifespan.2 <- update(non_cucks_male.lifespan.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.lifespan.2) #interaction term not significant.
vif(non_cucks_male.lifespan.2)
anova(non_cucks_male.lifespan.1,non_cucks_male.lifespan.2) #models not significantly different


#check dadF help interaction. 
non_cucks_male.lifespan.3 <- update(non_cucks_male.lifespan.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.lifespan.3) #interaction term not significant
vif(non_cucks_male.lifespan.3)
anova(non_cucks_male.lifespan.1,non_cucks_male.lifespan.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_male.lifespan.4 <- update(non_cucks_male.lifespan.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.lifespan.4) #interaction term not sig. Nearly though
vif(non_cucks_male.lifespan.3)
anova(non_cucks_male.lifespan.1,non_cucks_male.lifespan.3)  #models are not significantly different



#-------------------non_cucks_male mumage interaction-------


#check mumF MumAge interaction. 
non_cucks_male.lifespan.4 <- update(non_cucks_male.lifespan.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.lifespan.4) #interaction term not significant.
vif(non_cucks_male.lifespan.2)
anova(non_cucks_male.lifespan.1,non_cucks_male.lifespan.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_male.lifespan.5 <- update(non_cucks_male.lifespan.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.lifespan.5) #interaction term not sig .
vif(non_cucks_male.lifespan.3)
anova(non_cucks_male.lifespan.1,non_cucks_male.lifespan.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_male.lifespan.6 <- update(non_cucks_male.lifespan.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.lifespan.6) #interaction term not sig.
vif(non_cucks_male.lifespan.3)
anova(non_cucks_male.lifespan.1,non_cucks_male.lifespan.6)  #models are not significantly different



##-------------------cuck tests--------------------------

#select cuck offspring
cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for cucks_female  birds------------
cucks_female.lifespan.1 <- glmer(lifespan ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                  data = cucks_female_birds,
                                  family=binomial,
                                  control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(cucks_female.lifespan.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_female.lifespan.1, plot = T)

#check for collinearirty with vif
vif(cucks_female.lifespan.1)

##plot basic model with no interactions
#plot recruitment to lifespan for different Fs

cucks_female.lifespan.indF.plot <- ggplot(cucks_female_birds, aes(lifespan, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.lifespan.mumF.plot <- ggplot(cucks_female_birds, aes(lifespan,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.lifespan.dadF.plot <- ggplot(cucks_female_birds, aes(lifespan,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to lifespan", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.lifespan.indF.plot <- ggplot(male_birds, aes(lifespan, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.lifespan.mumF.plot <- ggplot(male_birds, aes(lifespan,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.lifespan.dadF.plot <- ggplot(male_birds, aes(lifespan,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to lifespan", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


lifespan.global <- (cucks_female.lifespan.indF.plot + male.lifespan.indF.plot) / (cucks_female.lifespan.mumF.plot + male.lifespan.mumF.plot) / (cucks_female.lifespan.dadF.plot + male.lifespan.dadF.plot)

ggsave("Figures/lifespan_global.png", lifespan.global)

#----------------cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
cucks_female.lifespan.2 <- update(cucks_female.lifespan.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                 optCtrl = list(maxfun = 100000000)))

summary(cucks_female.lifespan.2) #interaction term not significant.
vif(cucks_female.lifespan.2)
anova(cucks_female.lifespan.1,cucks_female.lifespan.2) #models not significantly different


#check dadF help interaction. 
cucks_female.lifespan.3 <- update(cucks_female.lifespan.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.lifespan.3) #interaction term very almost significant p = 0.0503 .
vif(cucks_female.lifespan.3)
anova(cucks_female.lifespan.1,cucks_female.lifespan.3)  #models are not significantly different


#check FROH help interaction. 
cucks_female.lifespan.4 <- update(cucks_female.lifespan.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.lifespan.4) #interaction term not sig.
vif(cucks_female.lifespan.3)
anova(cucks_female.lifespan.1,cucks_female.lifespan.3)  #models are significantly different

#----------------------------------cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
cucks_female.lifespan.4 <- update(cucks_female.lifespan.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.lifespan.4) #interaction term not significant.
vif(cucks_female.lifespan.2)
anova(cucks_female.lifespan.1,cucks_female.lifespan.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_female.lifespan.5 <- update(cucks_female.lifespan.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.lifespan.5) #interaction term not sig .
vif(cucks_female.lifespan.3)
anova(cucks_female.lifespan.1,cucks_female.lifespan.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_female.lifespan.6 <- update(cucks_female.lifespan.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.lifespan.6) #interaction term not sig.
vif(cucks_female.lifespan.3)
anova(cucks_female.lifespan.1,cucks_female.lifespan.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select bastard
cucks_male_birds <-filter(na.omit(male_birds, EPP == 0) )

#basic model
cucks_male.lifespan.1 <- glmer(lifespan ~ mumF + dadF +FROH + BrMF+  help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                                data = cucks_male_birds,
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.lifespan.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_male.lifespan.1, plot = T)

#check for collinearirty with vif
vif(cucks_male.lifespan.1)


#-----------------cucks_male helper interactions-----------





#check mumF help interaction. 
cucks_male.lifespan.2 <- update(cucks_male.lifespan.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.lifespan.2) #interaction term not significant.
vif(cucks_male.lifespan.2)
anova(cucks_male.lifespan.1,cucks_male.lifespan.2) #models not significantly different


#check dadF help interaction. 
cucks_male.lifespan.3 <- update(cucks_male.lifespan.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.lifespan.3) #interaction term not significant
vif(cucks_male.lifespan.3)
anova(cucks_male.lifespan.1,cucks_male.lifespan.3)  #models are not significantly different


#check FROH help interaction. 
cucks_male.lifespan.4 <- update(cucks_male.lifespan.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.lifespan.4) #interaction term not sig. Nearly though
vif(cucks_male.lifespan.3)
anova(cucks_male.lifespan.1,cucks_male.lifespan.3)  #models are not significantly different


#-------------------cucks_male mumage interaction-------

#check mumF MumAge interaction. 
cucks_male.lifespan.4 <- update(cucks_male.lifespan.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.lifespan.4) #interaction term not significant.
vif(cucks_male.lifespan.2)
anova(cucks_male.lifespan.1,cucks_male.lifespan.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_male.lifespan.5 <- update(cucks_male.lifespan.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.lifespan.5) #interaction term not sig .
vif(cucks_male.lifespan.3)
anova(cucks_male.lifespan.1,cucks_male.lifespan.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_male.lifespan.6 <- update(cucks_male.lifespan.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.lifespan.6) #interaction term not sig.
vif(cucks_male.lifespan.3)
anova(cucks_male.lifespan.1,cucks_male.lifespan.6)  #models are not significantly different









