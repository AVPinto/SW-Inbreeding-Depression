#---------------------Chapter 2.1 Recruitment to Adulthood Stats-----------------
#Where does mortality happen in the warbler lifespan? Start at the beginning.
#Do a GLMM on recruitment to adulthood

#--------------------------------libraries-------------


library(sjPlot)
library(simpleboot)
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
library(ggtext)

#----------------------------------import data--------------

##Factorise
stat.birds <- read.csv("Data/clean_data.csv")
stat.birds$h_countF <- as.factor(stat.birds$h_count)
stat.birds$help <- 0
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$mum <- as.factor(stat.birds$mum)
stat.birds$dad <- as.factor(stat.birds$dad)
stat.birds$birth_yearF <- as.factor(stat.birds$birth_year)
stat.birds$sib_pres <- as.factor(stat.birds$sib_pres)


# select males and females, birds that are already dead and never translocated

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

#first year survival---------------

#test if maternal or paternal roh affects the first year survival
#do this using 'adulthood' variable, which is 0 for bird died before 1 year old
#this is a success or failure, so should use binomial distribution

##--------------------------joint sex models----------------

#filter out translocated/alive birds

both_birds <- stat.birds %>% filter(event == 0,
                                    Coverage >= 0.25 ,
                                    CovMum >= 0.25 ,
                                    CovDad >= 0.25,
                                    CovBrM >= 0.25,
                                    )

dim(both_birds) #503 samples

#rescale variables for model convergence

both.fys.b <- glmer(unix_fys ~ rescale(LargeFROH) + 
                      rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                      h_countF +  rescale(birth_total_rain) + rescale(BirthRainCV) + 
                      Sex + rescale(natal_group_size) + sib_pres + EPP +
                      (1 | birth_year) + (1 | mum) + (1 | dad),
                    data = both_birds,
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

#is singular
summary(both.fys.b)
#random effects are incredibly small, remove them

both.fys.b <- glm(unix_fys ~ rescale(LargeFROH)+ 
                      rescale(BrMLargeFroh) + 
                      rescale(dadLargeFroh) + 
                      rescale(mumLargeFroh) +  
                      rescale(birth_total_rain) + rescale(BirthRainCV) + 
                      Sex + 
                      h_countF +
                      rescale(natal_group_size) + 
                      sib_pres +
                      EPP ,
                    data = both_birds,
                    family=binomial)

check_collinearity(both.fys.b)

check_overdispersion(both.fys.b)

check_singularity(both.fys.b)
#no longer singular

vif(both.fys.b) #vifs are low. BrM is at 1.3. Can include EPP

simulateResiduals(both.fys.b, plot = T) #looks good.



##------check for interactions---

#HELP

#BrM x help

both.fys.c <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*h_countF )

simulateResiduals(both.fys.b, plot = T) #looks good.

summary(both.fys.c) #ns

#dad x help

both.fys.d <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*h_countF )

simulateResiduals(both.fys.d, plot = T) #looks good.

summary(both.fys.d) #ns

#mum x help

both.fys.e <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*h_countF )

simulateResiduals(both.fys.e, plot = T) #looks good.

summary(both.fys.e) #not significant

#BIRTH ANNUAL RAIN

#BrM x birth annual rain

both.fys.f <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*rescale(birth_total_rain) )


simulateResiduals(both.fys.f, plot = T) #looks good.

summary(both.fys.f) #ns

#dad x birth annual rain

both.fys.g <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*rescale(birth_total_rain) )


simulateResiduals(both.fys.g, plot = T) #looks good.

summary(both.fys.g) #ns

#mum x birth annual rain

both.fys.h <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*rescale(birth_total_rain) )


simulateResiduals(both.fys.h, plot = T) #looks good.

summary(both.fys.h) #ns

#BIRTH RAIN VARIANCE

#BrM x birth rain variance

both.fys.i <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*rescale(BirthRainCV) )

simulateResiduals(both.fys.i, plot = T) #looks good.

summary(both.fys.i) #ns

#dad x birth rain variance

both.fys.j <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*rescale(BirthRainCV) )

simulateResiduals(both.fys.f, plot = T) #looks good.

summary(both.fys.j) #ns

#mum x birth rain variance

both.fys.k <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*rescale(BirthRainCV) )

simulateResiduals(both.fys.k, plot = T) #looks good.


summary(both.fys.k) #ns

#SEX interactions

#BrM x sex

both.fys.l <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*Sex )

simulateResiduals(both.fys.l, plot = T) #looks good.

summary(both.fys.l) #ns

#dad x sex

both.fys.m <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*Sex )

simulateResiduals(both.fys.m, plot = T) #looks good.

summary(both.fys.m) #not significant

anova(both.fys.b, both.fys.m)

#mum x sex

both.fys.n <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*Sex )

simulateResiduals(both.fys.n, plot = T) #looks good.

summary(both.fys.n) #ns

#NATAL GROUP SIZE


#BrM x natal group size

both.fys.o <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*rescale(natal_group_size) )

simulateResiduals(both.fys.o, plot = T) #slightly bad

summary(both.fys.o) #ns

#dad x natal_group_size

both.fys.p <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*rescale(natal_group_size) )

simulateResiduals(both.fys.p, plot = T) #looks good.

summary(both.fys.p) #ns!

anova(both.fys.p, both.fys.b) #ns!

#mum x rescale(natal_group_size)

both.fys.q <- update(both.fys.b,~ . + rescale(mumLargeFroh)*rescale(natal_group_size) )

simulateResiduals(both.fys.q, plot = T) #looks good.

summary(both.fys.q) #ns

#SIB PRESENCE

#BrM x sib_pres

both.fys.u <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*sib_pres )

simulateResiduals(both.fys.u, plot = T) #slightly bad

summary(both.fys.u) #ns

#dad x sib_pres

both.fys.v <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*sib_pres)

simulateResiduals(both.fys.v, plot = T) #looks good.

summary(both.fys.v) #ns

#mum x sib_pres

both.fys.w <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*sib_pres )

simulateResiduals(both.fys.w, plot = T) #looks good.

summary(both.fys.w) #significant!!!!!

anova(both.fys.b, both.fys.w) #imporves model fit!

#so in summary,  a mum froh * sib_presence(both.fys.t) interaction

##--output------------- 
#just the base model

summary(both.fys.b)

r2(both.fys.b)

tbl_regression(both.fys.b, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall during hatch year",
                            "rescale(BirthRainCV)" = "Variance in annual rainfall during hatch year",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH"),
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

#plot it with ggplot

fys.mum.ggplot <- ggplot(just_birds, aes(mumLargeFroh, unix_fys))+
  geom_smooth(method = "lm",
              col = "black",
              size = 1)+
  geom_point(shape = 1,
             size = 2)+
  scale_y_continuous(breaks = c(0,1))+
  labs(title = "Fig. 4 Effect of parental FROH on FYS",
       y = "",
       x = "Mother's FROH > 3.3MB")+
  theme_classic()+
  theme(text = element_text(size = 20)) 


fys.gendad.ggplot <- ggplot(just_birds, aes(dadLargeFroh, unix_fys))+
  geom_smooth(method = "lm",
              col = "black",
              size = 2)+
  geom_point(shape = 1,
             size = 1)+
  scale_y_continuous(breaks = c(0,1))+
  labs(title = "Fig. 4.1.
a",
       y = "",
       x = bquote("Genetic Father's F"[ROH]))+
  theme_classic()+
  theme(text = element_text(size = 20)) 

fys.socdad.ggplot <- ggplot(just_birds, aes(BrMLargeFroh, unix_fys))+
  geom_smooth(method = "lm",
              col = "black",
              size = 2)+
  geom_point(shape = 1,
             size = 1)+
  scale_y_continuous(breaks = c(0,1))+
  labs(title = "b",
       y = "Probability of FYS",
       x = bquote("Social Father's F"[ROH]))+
  theme_classic()+
  theme(text = element_text(size = 20)) 


fys.mum.ggplot <- ggplot(just_birds, aes(mumLargeFroh, unix_fys))+
  geom_smooth(method = "lm",
              col = "black",
              size = 2)+
  geom_point(shape = 1,
             size = 1)+
  scale_y_continuous(breaks = c(0,1))+
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5))+
  labs(title = "c",
       y = "",
       x = bquote("Mother's F"[ROH]))+
  theme_classic()+
  theme(text = element_text(size = 20)) 




(fys.gendad.ggplot / fys.socdad.ggplot / fys.mum.ggplot)




#significant interactions model - model w
summary(both.fys.w)

r2(both.fys.w)

tbl_regression(both.fys.w, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            "h_countF" = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall during hatch year",
                            "rescale(BirthRainCV)" = "Variance in annual rainfall during hatch year",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH",
                            "rescale(natal_group_size)"= "Natal group size",
                            "sib_pres"  = "Sibling presence",
                            "rescale(MumAge)" = "Mother's age at hatching"),
               estimate_fun = label_style_sigfig(digits = 3),
               pvalue_fun = label_style_pvalue(digits = 3)
)%>% 
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) %>%
  
  modify_fmt_fun(
    std.error  = function(x) style_sigfig(x, digits = 3),
    statistic  = function(x) style_sigfig(x, digits = 3)
  )


#plot interactions

#mum * sib


fys.mumXsib.ggplot <- ggplot(both_birds, aes(x =mumLargeFroh,
                                               y = unix_fys,
                                               group = sib_pres,
                                               colour = sib_pres))+
  geom_smooth(alpha = 0.2,
              method = "lm",
              size = 2,
              se = T)+
  geom_point(
             size = 2)+
  scale_y_continuous(breaks = c(0,1))+
  labs(title = "",
       y = "First Year Survival",
       x = bquote("Mother's F"[ROH]))+
  theme_classic()+
  theme(text = element_text(size = 20)) +
  theme(legend.position = "none")

fys.mumXsib.ggplot 

ggsave("Figures/Figure_3.jpg",
       plot = fys.mumXsib.ggplot,
       width = 10,
       height = 10,
       dpi = 600)

