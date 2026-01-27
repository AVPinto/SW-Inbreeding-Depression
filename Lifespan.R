#---------------------Chapter 2.3 GLMM Lifespan analysis-----------------
#Lets doing a GLMM on lifespan.

#--------------------------------libraries-------------


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
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$mum <- as.factor(stat.birds$mum)
stat.birds$dad <- as.factor(stat.birds$dad)
stat.birds$birth_yearF <- as.factor(stat.birds$birth_year)
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


dim(both_birds) 


#begin testing for families

both.lifespan.nbinom2 <- glmmTMB(lifespan ~ rescale(LargeFROH) + rescale(BrMLargeFroh)
                                 +rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                   Sex + h_countF + 
                                   rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                   sib_pres +
                                   natal_group_size +
                                   EPP +
                                   rescale(birth_year)+
                                   (1 | mum) ,
                                 data = both_birds,
                                 ziformula=~0,
                                 family=nbinom2
)



#nbinom2

simulateResiduals(both.lifespan.nbinom2, plot = T) # best?

check_collinearity(both.lifespan.nbinom2) #vifs around 1

check_overdispersion(both.lifespan.nbinom2) #none detected

check_convergence(both.lifespan.nbinom2)#yes!

check_singularity(both.lifespan.nbinom2) #false



summary(both.lifespan.nbinom2) #converges ok


#test interactions

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

summary(both.lifespan.nbinom2.b) #it's marginal

anova(both.lifespan.nbinom2,both.lifespan.nbinom2.b) #anova says almost?

#mum * sex

both.lifespan.nbinom2.c <- update(both.lifespan.nbinom2, ~ . + rescale(mumLargeFroh)*Sex) 

simulateResiduals(both.lifespan.nbinom2.c, plot = T) # looks a bit weird

check_collinearity(both.lifespan.nbinom2.c) #vifs around 1

check_overdispersion(both.lifespan.nbinom2.c) #none detected

summary(both.lifespan.nbinom2.c) #ns

#dadFROH * sex interaction maybe?

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
check_singularity(both.lifespan.nbinom2)


r2(both.lifespan.nbinom2)

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











#------------epp only models-------------


epp_birds <- filter(both_birds, EPP == 1)

str(epp_birds) #n = 204


epp_lifespan_glmm <- glmmTMB(lifespan ~ rescale(LargeFROH) + rescale(BrMLargeFroh)
                             +rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                               Sex + h_countF + 
                               rescale(mean_total_rain) + rescale(mean_rain_cv) +
                               sib_pres +
                               natal_group_size +
                               rescale(birth_year)+
                               (1 | mum) + (1 | dad) ,
                             data = epp_birds,
                             ziformula=~0,
                             family=poisson
)

simulateResiduals(epp_lifespan_glmm, plot = T)
check_convergence(epp_lifespan_glmm)
check_singularity(epp_lifespan_glmm)
check_collinearity(epp_lifespan_glmm)

#all looks ok

summary(epp_lifespan_glmm)

#still no parental effect significant

#output


tbl_regression(epp_lifespan_glmm, intercept = T,
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


