#---------------------Chapter 2.4 Lifetime reproductive success-----------------
#Do a GLMM on parental Fs and Lifetime reproductive success(LRS).
#LRS can be measured with the total number of offspring that survived > 1 year

#--------------------------------libraries-------------

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
#----------------------------------import data--------------

##Factorise
stat.birds <- read.csv("Data/clean_data.csv")
stat.birds$h_countF <- as.factor(stat.birds$h_count)
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$mum <- as.factor(stat.birds$mum)
stat.birds$dad <- as.factor(stat.birds$dad)
stat.birds$birth_yearF <- as.factor(stat.birds$birth_year)
stat.birds$sib_pres <- as.factor(stat.birds$sib_pres)



# select males and females, birds that are already dead, never translocated


female_birds <- stat.birds %>% filter(Sex == 0,
                                      event == 0,
                                      Coverage >= 0.25 ,
                                      CovMum >= 0.25 ,
                                      CovDad >= 0.25,
                                      CovBrM >= 0.25)

male_birds <- stat.birds %>% filter(Sex == 1, 
                                    event == 0,
                                    Coverage >= 0.25 ,
                                    CovMum >= 0.25 ,
                                    CovDad >= 0.25,
                                    CovBrM >= 0.25)


both_birds <- stat.birds %>% filter(event == 0,
                                    Coverage >= 0.25 ,
                                    CovMum >= 0.25 ,
                                    CovDad >= 0.25,
                                    CovBrM >= 0.25)



#-------------------modelling---------------

#test if maternal or paternal or social roh affects the lrs
#do this using 'n_off' variable, which is number of offspring that survived >1 year

hist(female_birds$n_off) #data has lots of zeros


both.lrs.poisson <- glmmTMB(n_off ~ rescale(LargeFROH) + 
                              rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                            Sex +
                            h_countF + 
                            rescale(mean_total_rain) + rescale(mean_rain_cv) +
                            EPP +  
                            sib_pres+
                              natal_group_size +
                            rescale(birth_year) + 
                            lifespan  ,
                            data = both_birds,
                            ziformula=~rescale(LargeFROH) + 
                              rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                              Sex +
                              h_countF + 
                              rescale(mean_total_rain) + rescale(mean_rain_cv) +
                              EPP +  
                              sib_pres+
                              natal_group_size+
                              rescale(birth_year) + 
                              lifespan,
                            family=poisson)


summary(both.lrs.poisson)

simulateResiduals(both.lrs.poisson, plot = T)  
#fails dispersion test but qq plots looks ok

#check for collinearity and overdispersion
check_collinearity(both.lrs.poisson) #vif around one
check_zeroinflation(both.lrs.poisson) #seems ok


#--------------test interactions-----------

#test for sex

#dadLargeFroh * Sex
both.lrs.poisson.a <- update(both.lrs.poisson, ~ . + rescale(dadLargeFroh)*Sex)

summary(both.lrs.poisson.a) #ns

simulateResiduals(both.lrs.poisson.a, plot = T) #seems ok

#brMfroh * Sex

both.lrs.poisson.b <- update(both.lrs.poisson, ~ . + rescale(BrMLargeFroh)*Sex)

summary(both.lrs.poisson.b) #ns but marginal 

simulateResiduals(both.lrs.poisson.b, plot = T) #seems ok

#mumfroh * Sex

both.lrs.poisson.c <- update(both.lrs.poisson, ~ . + rescale(mumLargeFroh)*Sex)

summary(both.lrs.poisson.c) #significant!

anova(both.lrs.poisson, both.lrs.poisson.c) #adds to model fit great

simulateResiduals(both.lrs.poisson.c, plot = T) #seems ok


#what is happening

ggplot(both_birds,aes(mumLargeFroh, n_off,group = Sex,colour = Sex))+
  geom_smooth(method = "lm")

#seems that females have +ve, males -ve, but both seem close to zero? will test seperately


#test against help

#rescale(dadLargeFroh) * help
both.lrs.poisson.d<- update(both.lrs.poisson.c, ~ . + rescale(dadLargeFroh)*h_countF ) #model convergence prob

summary(both.lrs.poisson.d) #ns

simulateResiduals(both.lrs.poisson.d, plot = T) #seems ok

#brMfroh * help

both.lrs.poisson.e <- update(both.lrs.poisson.c, ~ . + rescale(BrMLargeFroh)*h_countF)

summary(both.lrs.poisson.e) #ns

simulateResiduals(both.lrs.poisson.e, plot = T) #seems ok

#mumfroh * help

both.lrs.poisson.f <- update(both.lrs.poisson.c, ~ . + rescale(mumLargeFroh)*h_countF)

summary(both.lrs.poisson.f) #ns

anova(both.lrs.poisson, both.lrs.poisson.f) #ns

simulateResiduals(both.lrs.poisson.f, plot = T) #seems ok


#int with total rain

#rescale(dadLargeFroh) * total rain
both.lrs.poisson.g <- update(both.lrs.poisson.c, ~ . + rescale(dadLargeFroh)*rescale(mean_total_rain))

summary(both.lrs.poisson.g) #ns

simulateResiduals(both.lrs.poisson.g, plot = T) #seems ok

#brMfroh * total rain

both.lrs.poisson.h <- update(both.lrs.poisson.c, ~ . + rescale(BrMLargeFroh)*rescale(mean_total_rain))

summary(both.lrs.poisson.h) #ns  

simulateResiduals(both.lrs.poisson.h, plot = T) #seems ok

#mumfroh * mean_total_rain

both.lrs.poisson.i <- update(both.lrs.poisson.c, ~ . + rescale(mumLargeFroh)*rescale(mean_total_rain))

summary(both.lrs.poisson.i) #small effect size, marginal

simulateResiduals(both.lrs.poisson.i, plot = T) #seems ok



#test rain variance


#rescale(dadLargeFroh) * rain variance
both.lrs.poisson.j <- update(both.lrs.poisson.c, ~ . + rescale(dadLargeFroh)*rescale(mean_rain_cv))

summary(both.lrs.poisson.j) #ns
anova(both.lrs.poisson, both.lrs.poisson.j) #ns

simulateResiduals(both.lrs.poisson.j, plot = T) #seems ok


#brMfroh * mean_rain_cv

both.lrs.poisson.k <- update(both.lrs.poisson.c, ~ . + rescale(BrMLargeFroh)*rescale(mean_rain_cv))

summary(both.lrs.poisson.k) #ns 
anova(both.lrs.poisson, both.lrs.poisson.k) #ns

simulateResiduals(both.lrs.poisson.j, plot = T) #seems ok

#mumfroh * mean_rain_cv

both.lrs.poisson.l <- update(both.lrs.poisson.c, ~ . + rescale(mumLargeFroh)*rescale(mean_rain_cv))

summary(both.lrs.poisson.l) #ns

simulateResiduals(both.lrs.poisson.l, plot = T) #seems ok


#in summary, mumFROH x SEX interaction IS ALL

#--------output---


#output mumfroh x sex interaction
  
summary(both.lrs.poisson.c)

r2(both.lrs.poisson.c)

tbl_regression(both.lrs.poisson.c, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH", 
                            "h_count" = "Helper in natal territory",
                            "rescale(mean_total_rain)"	 = "Mean annual rainfall during lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance annual rainfall during lifespan",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH",
                            "rescale(birth_year)" = "Hatch Year"),
               estimate_fun = label_style_sigfig(digits = 3),
               pvalue_fun = label_style_pvalue(digits = 3)
)%>% 
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) %>%
  
  modify_fmt_fun(
    std.error  = function(x) style_sigfig(x, digits = 3),
    statistic  = function(x) style_sigfig(x, digits = 3)
  )

#plot basic parental froh with ggplot


#social dad Froh
plot.both.lrs.brm <-ggplot(both_birds, aes(BrMLargeFroh, n_off)) +
  geom_smooth(method = "lm")+
  geom_point()+
  theme_classic2()+
  labs( x = bquote("Social Father's F"[ROH]),
        y= "")+
  theme(text = element_text(size = 20))+
  ggtitle(" Fig. 6. 
a")


#gen male Froh
plot.both.lrs.gendad <-ggplot(both_birds, aes(dadLargeFroh, n_off)) +
  geom_smooth(method = "lm")+
  geom_point()+
  theme_classic2()+
  xlim(0,0.5)+
  labs( x = bquote("Genetic Father's F"[ROH]),
        y= "LRS")+
  theme(text = element_text(size = 20))+
  ggtitle("b")

#mumFroh
plot.both.lrs.mum<-ggplot(both_birds, aes(mumLargeFroh, n_off)) +
  geom_smooth(method = "lm")+
  geom_point()+
  xlim(0,0.5)+
  theme_classic2()+
  labs( x = bquote("Mother's F"[ROH]), 
        y= " ")+
  theme(text = element_text(size = 20))+
  ggtitle("c")


(plot.both.lrs.brm / plot.both.lrs.gendad) / (plot.both.lrs.mum )



#plot mother froh Sex

ggplot(both_birds,aes(mumLargeFroh, n_off , group = Sex, fill = Sex))+
  geom_smooth(method = "lm")+
  geom_jitter(aes(color  = Sex))+
  theme_classic2()+
 labs( x = bquote("Mother's F"[ROH]),
        y= "Productivity" ) +
  theme(text = element_text(size = 20)) +
  ggtitle("")



#------------epp only models-------------


epp_birds <- filter(both_birds, EPP == 1)

str(epp_birds) #n = 204


epp_lrs_glmm <- glmmTMB(n_off ~ rescale(LargeFROH) + rescale(BrMLargeFroh)
                             +rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                               Sex + 
                               rescale(mean_total_rain) + rescale(mean_rain_cv) +
                               sib_pres +
                               natal_group_size +
                               rescale(birth_year)+
                               (1 | mum) + (1 | dad) ,
                             data = epp_birds,
                             ziformula=~0,
                             family=poisson
)


simulateResiduals(epp_lrs_glmm, plot = T)

summary(epp_lrs_glmm)

#still no parental FROH effects significant



#------------test lrs split by sex-------

dim(female_birds)

female.lrs.poisson.1 <- glmmTMB(n_off ~ rescale(LargeFROH) + 
                                          rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                          rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                          sib_pres+
                                          rescale(birth_year) + 
                                          lifespan  ,
                                        data = female_birds,
                                        ziformula=~. ,
                                        family=poisson)


summary(female.lrs.poisson.1)

tbl_regression(female.lrs.poisson.1, intercept = T,
               show_single_row = "sib_pres",
               label = list("rescale(LargeFROH)" = "FROH", 
                            "h_count" = "Helper in natal territory",
                            "rescale(mean_total_rain)"	 = "Mean annual rainfall during lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance annual rainfall during lifespan",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH",
                            "rescale(birth_year)" = "Hatch Year"),
               estimate_fun = label_style_sigfig(digits = 3),
               pvalue_fun = label_style_pvalue(digits = 3)
)%>% 
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) %>%
  
  modify_fmt_fun(
    std.error  = function(x) style_sigfig(x, digits = 3),
    statistic  = function(x) style_sigfig(x, digits = 3)
  )

dim(male_birds)

male.lrs.poisson.1 <- glmmTMB(n_off ~ rescale(LargeFROH) + 
                                  rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                  rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                  sib_pres+
                                  rescale(birth_year) + 
                                  lifespan  ,
                                data = male_birds,
                                ziformula=~rescale(LargeFROH) + 
                                  rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                  rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                  sib_pres+
                                  rescale(birth_year) + 
                                  lifespan,
                                family=poisson)


summary(male.lrs.poisson.1)

tbl_regression(male.lrs.poisson.1, intercept = T,
               show_single_row = "sib_pres",
               label = list("rescale(LargeFROH)" = "FROH", 
                            "h_count" = "Helper in natal territory",
                            "rescale(mean_total_rain)"	 = "Mean annual rainfall during lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance annual rainfall during lifespan",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH",
                            "rescale(birth_year)" = "Hatch Year"),
               estimate_fun = label_style_sigfig(digits = 3),
               pvalue_fun = label_style_pvalue(digits = 3)
)%>% 
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) %>%
  
  modify_fmt_fun(
    std.error  = function(x) style_sigfig(x, digits = 3),
    statistic  = function(x) style_sigfig(x, digits = 3)
  )


