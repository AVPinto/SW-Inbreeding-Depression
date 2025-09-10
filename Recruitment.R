#---------------------Chapter 2.1 Recruitment to Adulthood Stats-----------------
#Where does mortality happen in the warbler lifespan? Start at the beginning.
#Do a GLMM on recruitment to adulthood

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
  install.packages(MuMIn)}

if (!is_installed("lme4")) { 
  install.packages(lmr4)}

if (!is_installed("patchwork")) { 
  install.packages("patchwork")}

install.packages("ggsignif")
install.packages("ggeffects")
install.packages("simpleboot")
install.packages("ggtext")

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
stat.birds$BreedSeason <- as.factor(stat.birds$BreedSeason)
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$DeathEventL <- as.logical(stat.birds$DeathEvent)
stat.birds$mum <- as.factor(stat.birds$mum)
stat.birds$dad <- as.factor(stat.birds$dad)
stat.birds$birth_yearF <- as.factor(stat.birds$birth_year)
stat.birds$adulthood <- as.factor(stat.birds$adulthood)
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

both.fys.b <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(BrMLargeFroh) + rescale(dadLargeFroh)
                    + rescale(mumLargeFroh) +
                      h_countF +  rescale(birth_total_rain) + rescale(BirthRainCV) + 
                      Sex + rescale(natal_group_size) + sib_pres + EPP +
                      (1 | birth_year) + (1 | mum) + (1 | dad),
                    data = both_birds,
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

#is singular
#mum and dad as random effects are incredibly small, remove them

both.fys.b <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(BrMLargeFroh) + rescale(dadLargeFroh)
                    + rescale(mumLargeFroh) +  
                        rescale(birth_total_rain) + rescale(BirthRainCV) + 
                      Sex + rescale(natal_group_size) + sib_pres + EPP +
                      (1 | birth_year),
                    data = both_birds,
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e9)))

#no longer singular

vif(both.fys.b) #vifs are low. BrM is at 1.3. Can include EPP

simulateResiduals(both.fys.b, plot = T) #looks good.

summary(both.fys.b)

##------check for interactions---

#HELP

#BrM x help

both.fys.c <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*h_count )

vif(both.fys.c) #vifs are ok

simulateResiduals(both.fys.b, plot = T) #looks good.

summary(both.fys.c) #ns

#dad x help

both.fys.d <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*help )

vif(both.fys.d) #vifs are ok

simulateResiduals(both.fys.d, plot = T) #looks good.

summary(both.fys.d) #ns

#mum x help

both.fys.e <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*help )

vif(both.fys.e) #vifs are ok


simulateResiduals(both.fys.e, plot = T) #looks good.

summary(both.fys.e) #not significant

#BIRTH ANNUAL RAIN

#BrM x birth annual rain

both.fys.f <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*rescale(birth_total_rain) )

vif(both.fys.f) #vifs are ok

simulateResiduals(both.fys.f, plot = T) #looks good.

summary(both.fys.f) #ns

#dad x birth annual rain

both.fys.g <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*rescale(birth_total_rain) )

vif(both.fys.g) #vifs are ok

simulateResiduals(both.fys.g, plot = T) #looks good.

summary(both.fys.g) #ns

#mum x birth annual rain

both.fys.h <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*rescale(birth_total_rain) )

vif(both.fys.h) #vifs are ok

simulateResiduals(both.fys.h, plot = T) #looks good.

summary(both.fys.h) #ns

#BIRTH RAIN VARIANCE

#BrM x birth rain variance

both.fys.i <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*rescale(BirthRainCV) )

vif(both.fys.i) #vifs are ok

simulateResiduals(both.fys.i, plot = T) #looks good.

summary(both.fys.i) #ns

#dad x birth rain variance

both.fys.j <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*rescale(BirthRainCV) )

vif(both.fys.j) #vifs are ok

simulateResiduals(both.fys.f, plot = T) #looks good.

summary(both.fys.j) #ns

#mum x birth rain variance

both.fys.k <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*rescale(BirthRainCV) )

vif(both.fys.k) #vifs are ok

simulateResiduals(both.fys.k, plot = T) #looks good.

summary(both.fys.k) #ns

#SEX interactions

#BrM x sex

both.fys.l <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*Sex )

vif(both.fys.l) #vifs are ok

simulateResiduals(both.fys.l, plot = T) #looks good.

summary(both.fys.l) #ns

#dad x sex

both.fys.m <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*Sex )

vif(both.fys.m) #vifs are ok

simulateResiduals(both.fys.m, plot = T) #looks good.

summary(both.fys.m) #not significant

anova(both.fys.b, both.fys.m)

#mum x sex

both.fys.n <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*Sex )

vif(both.fys.n) #vifs are ok

simulateResiduals(both.fys.n, plot = T) #looks good.

summary(both.fys.n) #ns

#NATAL GROUP SIZE


#BrM x natal group size

both.fys.o <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*rescale(natal_group_size) )

vif(both.fys.o) #vifs are ok

simulateResiduals(both.fys.o, plot = T) #slightly bad

summary(both.fys.o) #ns

#dad x natal_group_size

both.fys.p <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*rescale(natal_group_size) )

vif(both.fys.p) #vifs are ok

simulateResiduals(both.fys.p, plot = T) #looks good.

summary(both.fys.p) #ns!

anova(both.fys.p, both.fys.b) #ns!

#mum x rescale(natal_group_size)

both.fys.q <- update(both.fys.b,~ . + rescale(mumLargeFroh)*rescale(natal_group_size) )

vif(both.fys.q) #vifs are ok

simulateResiduals(both.fys.q, plot = T) #looks good.

summary(both.fys.q) #ns

#SIB PRESENCE

#BrM x sib_pres

both.fys.u <- update(both.fys.b, ~ . + rescale(BrMLargeFroh)*sib_pres )

vif(both.fys.u) #vifs are ok

simulateResiduals(both.fys.u, plot = T) #slightly bad

summary(both.fys.u) #ns

#dad x sib_pres

both.fys.v <- update(both.fys.b, ~ . + rescale(dadLargeFroh)*sib_pres)

vif(both.fys.v) #vifs are ok

simulateResiduals(both.fys.v, plot = T) #looks good.

summary(both.fys.v) #ns

#mum x sib_pres

both.fys.w <- update(both.fys.b, ~ . + rescale(mumLargeFroh)*sib_pres )

vif(both.fys.w) #vifs are ok

simulateResiduals(both.fys.w, plot = T) #looks good.

summary(both.fys.w) #significant!!!!!

#so in summary,  a mum froh * mum status(both.fys.t) interaction

##--output------------- 
#lets try just the base model first

summary(both.fys.b)

vif(both.fys.b)

tbl_regression(both.fys.b, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall during hatch year",
                            "rescale(BirthRainCV)" = "Variance in annual rainfall during hatch year",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH"
                            )
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  


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

tbl_regression(both.fys.w, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall during hatch year",
                            "rescale(BirthRainCV)" = "Variance in annual rainfall during hatch year",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH",
                            "rescale(natal_group_size)"= "Natal group size",
                            "sib_pres"  = "Sibling presence",
                            "rescale(MumAge)" = "Mother's age at hatching"
               )
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  


#plot interactions

#mum * sib


fys.mumXgroup.ggplot <- ggplot(both_birds, aes(x =mumLargeFroh,
                                               y = unix_fys,
                                               group = sib_pres,
                                               colour = sib_pres))+
  geom_smooth(alpha = 0.2,
              method = "lm",
              size = 2,
              se = T)+
  geom_point(shape = 1,
             size = 2)+
  scale_y_continuous(breaks = c(0,1))+
  labs(title = "Fig. 4.2.",
       y = "FYS",
       x = bquote("Mother's F"[ROH]))+
  theme_classic()+
  theme(text = element_text(size = 20)) 

fys.mumXgroup.ggplot 




#plots split by sex

plot.fys.males <- ggplot(filter(both_birds, Sex == 1), aes(dadLargeFroh, unix_fys))+
  geom_boxplot(fill = "blue")+
  geom_jitter(fill = "blue")+
  theme_classic2()+
  labs( x = "", y= " Males ")+
  xlim(0,0.5)+
  theme(text = element_text(size = 20))+
  ggtitle("Fig 5.2 FYS split by sex")+
  geom_signif(comparisons = list(c("0","1")),
              map_signif_level = TRUE)

plot.fys.females <- ggplot(filter(both_birds, Sex == 0), aes(dadLargeFroh, unix_fys))+
  geom_boxplot(fill = "pink")+
  geom_jitter(fill = "pink")+
  theme_classic2()+
  labs( x = "Genetic Father FROH", y= " Females ")+
  xlim(0,0.5)+
  theme(text = element_text(size = 20))+
  geom_signif(comparisons = list(c("0","1")),
              map_signif_level = TRUE)

plot.fys.males / plot.fys.females

#it seems that males are affected by gen dads, but not females



#--------------------do these models split by sex------------

#males, fys ~ FROH 


males.fys.a <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                      help + rescale(birth_total_rain) + rescale(BirthRainCV) + natal_group_size + mum_status + sib_pres + 
                      (1 | birth_year) + (1 | mum) + (1 | dad),
                    data = male_birds,
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

#fit is singular, try removing some random effect structure
males.fys.a <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                       help + rescale(birth_total_rain) + rescale(BirthRainCV) + natal_group_size + mum_status + sib_pres + EPP+ rescale(MumAge) +
                       (1 | birth_year) ,
                     data = male_birds,
                     family=binomial,
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

check_singularity(males.fys.a) #no more singularity


simulateResiduals(males.fys.a, plot = T)

check_overdispersion(males.fys.a)

check_collinearity(males.fys.a)
#looks ok

summary(males.fys.a)

#gen dad  significant.
#natal group size almost significant
#sib presence is almost significant
#mum status is significant

#output males 
summary(males.fys.a)


tbl_regression(males.fys.a, intercept = T,
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall during hatch year",
                            "rescale(BirthRainCV)" = "Variance in annual rainfall during hatch year",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH"
               )
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  





#try interactions

#with birth_total_rain birth_total_rainbirth_total_rainbirth_total_rainbirth_total_rainbirth_total_rain

#gendad with birth total rain
males.fys.b <- update(males.fys.a, ~ . + rescale(dadLargeFroh)*rescale(birth_total_rain))

summary(males.fys.b) #ns

#socdad with birth total rain
males.fys.c <- update(males.fys.a, ~ . + rescale(BrMLargeFroh)*rescale(birth_total_rain))

summary(males.fys.c) #ns

#mum with birth total rain
males.fys.d <- update(males.fys.a, ~ . + rescale(mumLargeFroh)*rescale(birth_total_rain))

summary(males.fys.d) #ns



#with rain variance rain variance rain variance  rain variance rain variance rain variance

#gendad with BirthRainCV
males.fys.e <- update(males.fys.a, ~ . + rescale(dadLargeFroh)*rescale(BirthRainCV))

summary(males.fys.e) #ns

#socdad with BirthRainCV
males.fys.f <- update(males.fys.a, ~ . + rescale(BrMLargeFroh)*rescale(BirthRainCV))

summary(males.fys.f) #ns

#mum with BirthRainCV
males.fys.g <- update(males.fys.a, ~ . + rescale(mumLargeFroh)*rescale(BirthRainCV))

summary(males.fys.g) #ns


# help helphelphelphelphelp help  help help help help help help help HELP

#gendad with help
males.fys.h <- update(males.fys.a, ~ . + rescale(dadLargeFroh)*rescale(help))

summary(males.fys.h) #ns

#socdad with help
males.fys.i <- update(males.fys.a, ~ . + rescale(BrMLargeFroh)*rescale(help))

summary(males.fys.i) #ns

#mum with help
males.fys.j <- update(males.fys.a, ~ . + rescale(mumLargeFroh)*rescale(help))

summary(males.fys.j) #ns




# natal_group_size  natal_group_size natal_group_size natal_group_size natal_group_size natal_group_size natal_group_size 

#gendad with natal_group_size 
males.fys.h <- update(males.fys.a, ~ . + rescale(dadLargeFroh)*rescale(natal_group_size ))

summary(males.fys.h) #ns

#socdad with natal_group_size 
males.fys.i <- update(males.fys.a, ~ . + rescale(BrMLargeFroh)*rescale(natal_group_size ))

summary(males.fys.i) #ns

#mum with natal_group_size 
males.fys.j <- update(males.fys.a, ~ . + rescale(mumLargeFroh)*rescale(natal_group_size ))

summary(males.fys.j) #ns



#mum_status vmum_status mum_status mum_status mum_status v mum_status

#gendad with mum_status
males.fys.k <- update(males.fys.a, ~ . + rescale(dadLargeFroh)*mum_status)

summary(males.fys.k) #ns

#socdad with mum_status
males.fys.l <- update(males.fys.a, ~ . + rescale(BrMLargeFroh)*mum_status)

summary(males.fys.l) #ns

#mum with mum_status
males.fys.m <- update(males.fys.a, ~ . + rescale(mumLargeFroh)*mum_status)

summary(males.fys.m) #ns


#sib_pres sib_pres sib_pres sib_pres sib_pres sib_pres sib_pres sib_pres

#gendad with sib_pres
males.fys.n <- update(males.fys.a, ~ . + rescale(dadLargeFroh)*sib_pres)

summary(males.fys.n) #ns

#socdad with sib_pres
males.fys.o <- update(males.fys.a, ~ . + rescale(BrMLargeFroh)*sib_pres)

summary(males.fys.o) #ns

#mum with sib_pres
males.fys.p <- update(males.fys.a, ~ . + rescale(mumLargeFroh)*sib_pres)

summary(males.fys.p) #marginally significant!!!!


#EPP EPP EPP EPP EPP EPP EPP EPP!

#gendad with EPP
males.fys.q <- update(males.fys.a, ~ . + rescale(dadLargeFroh)*EPP)

summary(males.fys.q) #ns

#socdad with EPP
males.fys.r <- update(males.fys.a, ~ . + rescale(BrMLargeFroh)*EPP)

summary(males.fys.r) #ns

#mum with EPP
males.fys.s <- update(males.fys.a, ~ . + rescale(mumLargeFroh)*EPP)

summary(males.fys.s) #ns


#MumAge MumAge MumAge MumAge MumAge MumAge MumAge MumAge

#gendad with MumAge
males.fys.t <- update(males.fys.a, ~ . + rescale(dadLargeFroh)*MumAge)

summary(males.fys.t) #ns

#socdad with MumAge
males.fys.u <- update(males.fys.a, ~ . + rescale(BrMLargeFroh)*MumAge)

summary(males.fys.u) #ns

#mum with MumAge
males.fys.v <- update(males.fys.a, ~ . + rescale(mumLargeFroh)*MumAge)

summary(males.fys.v) #ns


#in summary, none of these interactions are significant in males only.

#----------------------EPP interaction messaround-----------

#regarding the EPP interactions....

summary(males.fys.q)

vif(males.fys.q) #it's a bit high with BrMFrohLarge and the interaction

vif(males.fys.r)




#so we could try this model with with only one kind of dad when testing for interactions

#try this model only with genetic dad

males.fys.gendad <- glmer(unix_fys ~ rescale(LargeFROH)  + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                       help + rescale(birth_total_rain) + rescale(BirthRainCV) + natal_group_size + mum_status + sib_pres + EPP+ rescale(MumAge) +
                       (1 | birth_year),
                     data = male_birds,
                     family=binomial,
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

summary(males.fys.gendad) #gen dad is significant, EPP not.

#try with interaction

males.fys.gendad.int <- update(males.fys.gendad, ~ .  + rescale(dadLargeFroh)*EPP )

summary(males.fys.gendad.int) #interaction NOT significant


#try with only social dad

males.fys.socdad <- glmer(unix_fys ~ rescale(LargeFROH)  + rescale(BrMLargeFroh) + rescale(mumLargeFroh) +
                       help + rescale(birth_total_rain) + rescale(BirthRainCV) + natal_group_size + mum_status + sib_pres + EPP+ rescale(MumAge) +
                       (1 | birth_year),
                     data = male_birds,
                     family=binomial,
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

summary(males.fys.socdad) #BrM not significant

#try with interaction

males.fys.socdad.int <- update(males.fys.socdad, ~ .  + rescale(BrMLargeFroh)*EPP )

summary(males.fys.socdad.int) #interaction IS significant!

tbl_regression(males.fys.j, intercept = T,
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall during hatch year",
                            "rescale(BirthRainCV)" = "Variance in annual rainfall during hatch year",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH"
               )
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  


#PLOT  BRM ONLY , EP + WP included

plot.fys.males.BRM.EPP <- ggplot(filter(male_birds, EPP == 1), aes(BrMLargeFroh, unix_fys))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic2()+
  labs( x = "", y= " EP offspring ")+
  xlim(0,0.5)+
  theme(text = element_text(size = 20))+
  ggtitle("Fig 5.4 Male offspring FYS")+
  geom_signif(comparisons = list(c("0","1")),
              map_signif_level = TRUE)

plot.fys.males.BRM.EPP


plot.fys.males.BRM.WPP <- ggplot(filter(male_birds, EPP == 0), aes(BrMLargeFroh, unix_fys))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic2()+
  labs( x = "Social Father FROH", y= " WP ofspring ")+
  xlim(0,0.5)+
  theme(text = element_text(size = 20))+
  geom_signif(comparisons = list(c("0","1")),
              map_signif_level = TRUE)

plot.fys.males.BRM.WPP

plot.fys.males.BRM.EPP  / plot.fys.males.BRM.WPP









#--------EPP split-----------

#split model into non EPP and EPP

wp_males <- filter(male_birds, EPP == 0)
dim(wp_males) #n = 201

epp_males <- filter(male_birds, EPP == 1)
dim(epp_males) #n =159


wp_males.fys.1 <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                        help + rescale(birth_total_rain) + rescale(BirthRainCV) +
                                        natal_group_size + mum_status + sib_pres + rescale(MumAge) +
                                        (1 | birth_year) ,
                                      data = wp_males,
                                      family=binomial,
                                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
summary(wp_males.fys.1)
simulateResiduals(wp_males.fys.1, plot = T)


#in this wp_males dataset. dads are both social and genetic. DadFROH has an effect here


#try epp, removing the mum and dad random effects

epp_males.fys.1 <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(dadLargeFroh) + rescale(BrMLargeFroh) + rescale(mumLargeFroh) +
                          help + rescale(birth_total_rain) + rescale(BirthRainCV) +sib_pres + natal_group_size + mum_status +
                          (1 | birth_year) ,
                        data = epp_males,
                        family=binomial,
                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))


summary(epp_males.fys.1)

check_collinearity(epp_males.fys.1)
isSingular(epp_males.fys.1)
simulateResiduals(epp_males.fys.1, plot = T)


#no significance of genetic dad in this dataset of extra pair offspring?

#test if there is an interaction between them males?

epp_males.fys.2 <- update(epp_males.fys.1, ~ . + rescale(dadLargeFroh)*rescale(BrMLargeFroh))

summary(epp_males.fys.2) #ns

#interestingly, sib_pres also becomes unsignificant. it was the same direction and magnitude as gen, dad but know both are unsig.

#is this a power thing??? lets find out......


#try bootstrapping at the same sample size of the ep dataset.

#toy bootstrap model
rand_birds.1 <- male_birds[sample(nrow(male_birds), 150),]

rand.fys.1 <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(dadLargeFroh) + rescale(BrMLargeFroh) + rescale(mumLargeFroh) +
                      help + rescale(birth_total_rain) + rescale(BirthRainCV) +sib_pres + natal_group_size + mum_status +
                      (1 | birth_year) ,
                    data = rand_birds.1,
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))


summary(rand.fys.1)


dim(rand_birds.1)


#some ai bullshit

library(lme4)  # for glmer
library(glmmboot)  # if you still want to use parts of it

# Set bootstrap parameters
n_resamples <- 100  # number of bootstrap iterations
sample_size <- 150   # number of rows to sample in each iteration

# Manual bootstrap
boot_coefs <- lapply(1:n_resamples, function(i) {
  # Sample rows with replacement
  sampled_rows <- sample(nrow(male_birds), size = sample_size, replace = FALSE)
  boot_data <- male_birds[sampled_rows, ]
  
  # Refit model (with error handling)
  tryCatch({
    boot_model <- rand.fys.1 <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(dadLargeFroh) + rescale(BrMLargeFroh) + rescale(mumLargeFroh) +
                                        help + rescale(birth_total_rain) + rescale(BirthRainCV) +sib_pres + natal_group_size + mum_status +
                                        (1 | birth_year) ,
                                      data = boot_data,
                                      family=binomial,
                                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))
    fixef(boot_model)  # return fixed effects
  }, error = function(e) NULL)  # skip if model fails
})


# Remove NULLs (failed fits)
boot_coefs <- do.call(rbind, boot_coefs[!sapply(boot_coefs, is.null)]) 



boot_coefs <- data.frame(boot_coefs)


t.test(boot_coefs$rescale.dadLargeFroh.)

#so after 5000 tests

# mean of x 
#-0.7058043
#95 percent confidence interval:
#  -0.7185442 -0.6930643

range(boot_coefs$rescale.dadLargeFroh.)

hist(boot_coefs$rescale.dadLargeFroh.)













#what about females?

  
wp_females <- filter(female_birds, EPP == 0)
dim(wp_females) #n = 188

epp_females <- filter(female_birds, EPP == 1)
dim(epp_females) #n =151

#try the slimmest model
wp_females.fys.1 <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                           + (1 | birth_year) ,
                        data = wp_females,
                        family=binomial,
                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

summary(wp_females.fys.1) #genetic dad not significant. individual is though
simulateResiduals(wp_males.fys.1, plot = T)
  
#epp females  
  
epp_females.fys.1 <- glmer(unix_fys ~ rescale(LargeFROH) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                            + (1 | birth_year) ,
                          data = epp_females,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

summary(epp_females.fys.1) #genetic dad not significant
simulateResiduals(epp_males.fys.1, plot = T)

  
#output and plot EPP split

#table 5.3 - males wp
  
summary(wp_males.fys.2)  

tbl_regression(wp_males.fys.2, intercept = T,
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall during hatch year",
                            "rescale(BirthRainCV)" = "Variance in annual rainfall during hatch year",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH"
               )
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  


#table 5.4 - males epp

summary(epp_males.fys.1)  


tbl_regression(epp_males.fys.1, intercept = T,
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall during hatch year",
                            "rescale(BirthRainCV)" = "Variance in annual rainfall during hatch year",
                            "rescale(BrMLargeFroh)" = "Social Father FROH",
                            "rescale(dadLargeFroh)" = "Genetic Father FROH",
                            "rescale(mumLargeFroh)" = "Mother FROH"
               )
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  


#plot these

plot.fys.wp.males <- ggplot(wp_males, aes(dadLargeFroh, unix_fys))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic2()+
  labs( x = "Father FROH", y= "WP Offspring")+
  xlim(0,0.5)+
  theme(text = element_text(size = 20))+
  geom_signif(comparisons = list(c("0","1")),
              map_signif_level = TRUE)+
  ggtitle("Fig 5.3 Effects of paternal FROH on FYS in EPP and WP offspring")


plot.fys.ep.males.genetic <- ggplot(epp_males, aes(dadLargeFroh, unix_fys))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic2()+
  labs( x = "Genetic Father FROH", y= " EPP Offspring")+
  xlim(0,0.5)+
  theme(text = element_text(size = 20))+
  geom_signif(comparisons = list(c("0","1")),
              map_signif_level = TRUE)


plot.fys.ep.males.social <- ggplot(epp_males, aes(BrMLargeFroh, unix_fys))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic2()+
  labs( x = "Social Father FROH", y= "")+
  xlim(0,0.5)+
  theme(text = element_text(size = 20))+
  geom_signif(comparisons = list(c("0","1")),
              map_signif_level = TRUE)


plot.fys.wp.males / (plot.fys.ep.males.genetic | plot.fys.ep.males.social)



  








































  

#sgeom_violindot()#shit model


shit.unix_fys.1 <- glmer(unix_fys ~ LargeFROH + meanInsects +  (1 | birth_year) + (1 | mum) + (1 | dad),
                          data = both_birds,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

summary(shit.unix_fys.1)

#no interactions model

both.unix_fys.1 <- glmer(unix_fys ~ mumF + dadF + FROH + help + MumAge + Sex + BirthInsect + (1 | birth_year) + (1 | mum) + (1 | dad),
                          data = both_birds,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

simulateResiduals(both.adulthood.1, plot = T)
vif(both.adulthood.1)

summary(both.adulthood.1)  #nothing significant when dataset is not split by EPP


#interactions models

#mumfroh x birthinsect
both.adulthood.BirthInsectXMumF <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + Sex + BirthInsect + BirthInsect*mumF + (1 | birth_year) + (1 | mum) + (1 | dad),
                          data = both_birds,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

simulateResiduals(both.adulthood.BirthInsectXMumF, plot = T)
vif(both.adulthood.BirthInsectXMumF)

summary(both.adulthood.BirthInsectXMumF)

#froh x birthinsect

both.adulthood.BirthInsectXFROH <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + Sex + BirthInsect + BirthInsect*FROH + (1 | birth_year) + (1 | mum) + (1 | dad),
                          data = both_birds,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

simulateResiduals(both.adulthood.BirthInsectXFROH, plot = T)
vif(both.adulthood.BirthInsectXFROH)

summary(both.adulthood.BirthInsectXFROH)


#dadfroh x birthinsect

both.adulthood.BirthInsectXDadF <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + Sex + BirthInsect + dadF*BirthInsect + (1 | birth_year) + (1 | mum) + (1 | dad),
                                         data = both_birds,
                                         family=binomial,
                                         control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

simulateResiduals(both.adulthood.BirthInsectXDadF, plot = T)
vif(both.adulthood.BirthInsectXDadF)

summary(both.adulthood.BirthInsectXDadF)


#plot
ggplot(both_birds, aes(x = FROH, y = adulthood, group = Sex, color = Sex))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "glm")+
  labs( x = "Individual FROH", y = "First Year Survival", title = "All types of Parentage")+
  scale_color_discrete(labels=c('Female', 'Male'))+
  scale_x_continuous(limits = c(0.05,0.5), breaks = c(0.1,0.2,0.3,0.4))




#split birds into EPP and not EPP 

cuck_both_birds <- both_birds %>% filter(EPP == 1)
dim(cuck_both_birds) #222 birds from EPP

non_cuck_both_birds <- both_birds %>% filter(EPP == 0)
dim(non_cuck_both_birds) #312 birds not from EPP


#not EPP tests.
#We need to test if sex has interaction with mumF and dadF

#basic model, no interactions
non_cuck_both_adulthood_base <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + Sex + BirthInsect + (1 | birth_year) + (1 | mum) + (1 | dad),
                         data = non_cuck_both_birds,
                         family=binomial,
                         control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

summary(non_cuck_both_adulthood_base) #mumage and sex are significant

#plot just this 
ggplot(non_cuck_both_birds, aes(y = adulthood, x = mumF))+
  geom_point()+
  geom_()

#plot sex difference




ggplot(non_cuck_both_birds %>% 
  mutate(adulthoodN = as.numeric(as.character(adulthood)) ) %>%
  group_by(Sex) %>%
  summarise(mean = mean(adulthoodN),  sd(adulthoodN), se = (sd(adulthoodN)/sqrt(n()))),
  aes(x= Sex, y = mean) )+
  geom_bar(stat = "identity")


#test interactions

#test FROH x Sex interaction

non_cuck_both_adulthood_frohXsex <- update(non_cuck_both_adulthood_base, ~. +FROH*Sex ,control=glmerControl(optimizer="bobyqa"))

summary(non_cuck_both_adulthood_frohXsex) #interaction not significant

#test mumF x Sex interaction

non_cuck_both_adulthood_mumFXsex <- update(non_cuck_both_adulthood_base, ~. +mumF*Sex ,control=glmerControl(optimizer="bobyqa"))

summary(non_cuck_both_adulthood_mumFXsex) #interaction not significant

#test dadF x Sex interaction

non_cuck_both_adulthood_dadFXsex <- update(non_cuck_both_adulthood_base, ~. +dadF*Sex ,control=glmerControl(optimizer="bobyqa"))

summary(non_cuck_both_adulthood_dadFXsex) #interaction is significant!

#plot this whole deal 

dadfxSex_adulthood_plot <- ggplot(non_cuck_both_birds, aes(x = dadF, y = adulthood, group = Sex, color = Sex))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "glm")+
  labs( x = "Sire FROH", y = "First Year Survival")+
  scale_color_discrete(labels=c('Female', 'Male'))+
  scale_x_continuous(limits = c(0.05,0.5), breaks = c(0.1,0.2,0.3,0.4))

dadfxSex_adulthood_plot

mumfxSex_adulthood_plot <- ggplot(non_cuck_both_birds, aes(x = mumF, y = adulthood, group = Sex, color = Sex))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "glm")+
  labs( x = "Dam FROH", y = "First Year Survival")+
  scale_color_discrete(labels=c('Female', 'Male'))+
  scale_x_continuous(limits = c(0.05,0.5), breaks = c(0.1,0.2,0.3,0.4))

mumfxSex_adulthood_plot

indfxSex_adulthood_plot <- ggplot(non_cuck_both_birds, aes(x = FROH, y = adulthood, group = Sex, color = Sex))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "glm")+
  labs( x = "Individual FROH", y = "First Year Survival", title = "No EPP")+
  scale_color_discrete(labels=c('Female', 'Male'))+
  scale_x_continuous(limits = c(0.05,0.5), breaks = c(0.1,0.2,0.3,0.4))

indfxSex_adulthood_plot

indfxSex_adulthood_plot / mumfxSex_adulthood_plot / dadfxSex_adulthood_plot 



#EPP tests.
#here we will also test with BrM as its is different to sire


#basic model, no interactions
cuck_both_adulthood_base <- glmer(adulthood ~ mumF + dadF + FROH + BrMF + help + MumAge + Sex + BirthInsect + (1 | birth_year) + (1 | mum) + (1 | dad),
                                      data = cuck_both_birds,
                                      family=binomial,
                                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8)))

summary(cuck_both_adulthood_base) #mumage and sex are significant

#test interactions

#test FROH x Sex interaction

cuck_both_adulthood_frohXsex <- update(cuck_both_adulthood_base, ~. +FROH*Sex ,control=glmerControl(optimizer="bobyqa"))

summary(cuck_both_adulthood_frohXsex) #interaction not significant

#test mumF x Sex interaction

cuck_both_adulthood_mumFXsex <- update(cuck_both_adulthood_base, ~. +mumF*Sex ,control=glmerControl(optimizer="bobyqa"))

summary(cuck_both_adulthood_mumFXsex) #interaction not significant

#test dadF x Sex interaction

cuck_both_adulthood_dadFXsex <- update(cuck_both_adulthood_base, ~. +dadF*Sex ,control=glmerControl(optimizer="bobyqa"))

summary(cuck_both_adulthood_dadFXsex) #interaction is not significant!

#test BrMF x Sex interaction

cuck_both_adulthood_BrMFXsex <- update(cuck_both_adulthood_base, ~. +BrMF*Sex ,control=glmerControl(optimizer="bobyqa"))

summary(cuck_both_adulthood_dadFXsex) #interaction is not significant!


#plot this whole deal 

cuck_dadfxSex_adulthood_plot <- ggplot(cuck_both_birds, aes(x = dadF, y = adulthood, group = Sex, color = Sex))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "glm")+
  labs( x = "Genetic sire FROH", y = "First Year Survival")+
  scale_color_discrete(labels=c('Female', 'Male'))+
  scale_x_continuous(limits = c(0.05,0.5), breaks = c(0.1,0.2,0.3,0.4))+
  theme(legend.position="none")

cuck_dadfxSex_adulthood_plot

cuck_mumfxSex_adulthood_plot <- ggplot(cuck_both_birds, aes(x = mumF, y = adulthood, group = Sex, color = Sex))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "glm")+
  labs( x = "Dam FROH", y = "")+
  scale_color_discrete(labels=c('Female', 'Male'))+
  scale_x_continuous(limits = c(0.05,0.5), breaks = c(0.1,0.2,0.3,0.4))

cuck_mumfxSex_adulthood_plot

cuck_indfxSex_adulthood_plot <- ggplot(cuck_both_birds, aes(x = FROH, y = adulthood, group = Sex, color = Sex))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "glm")+
  labs( x = "Individual FROH", y = "First Year Survival")+
  scale_color_discrete(labels=c('Female', 'Male'))+
  scale_x_continuous(limits = c(0.05,0.5), breaks = c(0.1,0.2,0.3,0.4))+
  theme(legend.position="none")

cuck_indfxSex_adulthood_plot

cuck_BrMfxSex_adulthood_plot <- ggplot(cuck_both_birds, aes(x = BrMF, y = adulthood, group = Sex, color = Sex))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "glm")+
  labs( x = "Social Sire FROH", y = "")+
  scale_color_discrete(labels=c('Female', 'Male'))+
  scale_x_continuous(limits = c(0.05,0.5), breaks = c(0.1,0.2,0.3,0.4))

cuck_BrMfxSex_adulthood_plot

(cuck_indfxSex_adulthood_plot | cuck_mumfxSex_adulthood_plot )/ (cuck_dadfxSex_adulthood_plot | cuck_BrMfxSex_adulthood_plot)





##---------------------------split models by sex------
#
#
#
#
#
#
#-------------------------tests for female birds------------
female.adulthood.1 <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
      data = female_birds,
      family=binomial,
      control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(female.adulthood.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = female.adulthood.1, plot = T)

#check for collinearirty with vif
vif(female.adulthood.1)

##plot basic model with no interactions
#plot recruitment to adulthood for different Fs

female.adulthood.indF.plot <- ggplot(female_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.adulthood.mumF.plot <- ggplot(female_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.adulthood.dadF.plot <- ggplot(female_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to Adulthood", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.indF.plot <- ggplot(male_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.mumF.plot <- ggplot(male_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.dadF.plot <- ggplot(male_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to Adulthood", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


adulthood.global <- (female.adulthood.indF.plot + male.adulthood.indF.plot) / (female.adulthood.mumF.plot + male.adulthood.mumF.plot) / (female.adulthood.dadF.plot + male.adulthood.dadF.plot)

ggsave("Figures/adulthood_global.png", adulthood.global)

#----------------female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
female.adulthood.2 <- update(female.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.2) #interaction term not significant.
vif(female.adulthood.2)
anova(female.adulthood.1,female.adulthood.2) #models not significantly different


#check dadF help interaction. 
female.adulthood.3 <- update(female.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.3) #interaction term very almost significant p = 0.0503 .
vif(female.adulthood.3)
anova(female.adulthood.1,female.adulthood.3)  #models are not significantly different

#visually asses the difference.
#It seems that having helpers means that birds are LESS likely to recruit to adulthood when dads are inbred
female.adulthood.dadF.helpYES.plot <- ggplot(female_birds[female_birds$help == 1,], aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  scale_y_continuous(limits = c(0, 0.25))+
  labs( x = "Recruitment to Adulthood", y = "Paternal F", title = "Helpers in Nest")+
  theme_bw()+
  theme(text = element_text(size = 20))

female.adulthood.dadF.helpNO.plot <- ggplot(female_birds[female_birds$help == 0,], aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  scale_y_continuous(limits = c(0, 0.25))+
  labs( x = "Recruitment to Adulthood", y = "Paternal F", title = "No Helpers")+
  theme_bw()+
  theme(text = element_text(size = 20))

adulthood.female.help <- female.adulthood.dadF.helpYES.plot | female.adulthood.dadF.helpNO.plot 

ggsave("Figures/adulthood_female_help.png", adulthood.female.help)

#check FROH help interaction. 
female.adulthood.4 <- update(female.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.4) #interaction term not sig.
vif(female.adulthood.3)
anova(female.adulthood.1,female.adulthood.3)  #models are significantly different

#----------------------------------female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
female.adulthood.4 <- update(female.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.4) #interaction term not significant.
vif(female.adulthood.2)
anova(female.adulthood.1,female.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
female.adulthood.5 <- update(female.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.5) #interaction term not sig .
vif(female.adulthood.3)
anova(female.adulthood.1,female.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
female.adulthood.6 <- update(female.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.6) #interaction term not sig.
vif(female.adulthood.3)
anova(female.adulthood.1,female.adulthood.6)  #models are not significantly different

#-------------------------tests for male birds------------
#basic model
male.adulthood.1 <- glmer(adulthood ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                            data = male_birds,
                            family=binomial,
                            control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.1) #here FROH and mumage arte sig

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = male.adulthood.1, plot = T)

#check for collinearirty with vif
vif(male.adulthood.1)

#-----------------male helper interactions-----------

#check mumF help interaction. 
male.adulthood.2 <- update(male.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.2) #interaction term not significant.
vif(male.adulthood.2)
anova(male.adulthood.1,male.adulthood.2) #models not significantly different


#check dadF help interaction. 
male.adulthood.3 <- update(male.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.3) #interaction term not significant
vif(male.adulthood.3)
anova(male.adulthood.1,male.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
male.adulthood.4 <- update(male.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.4) #interaction term not sig. Nearly though
vif(male.adulthood.3)
anova(male.adulthood.1,male.adulthood.3)  #models are not significantly different


#-------------------male mumage interaction-------

#check mumF MumAge interaction. 
male.adulthood.4 <- update(male.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.4) #interaction term not significant.
vif(male.adulthood.2)
anova(male.adulthood.1,male.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
male.adulthood.5 <- update(male.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.5) #interaction term not sig .
vif(male.adulthood.3)
anova(male.adulthood.1,male.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
male.adulthood.6 <- update(male.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.6) #interaction term not sig.
vif(male.adulthood.3)
anova(male.adulthood.1,male.adulthood.6)  #models are not significantly different


##--------------------------------EPP tests-------------------------------------


#select non cuck offspring
non_cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for non_cucks_female  birds------------
non_cucks_female.adulthood.1 <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                            data = non_cucks_female_birds,
                            family=binomial,
                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(non_cucks_female.adulthood.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_female.adulthood.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_female.adulthood.1)

##plot basic model with no interactions
#plot recruitment to adulthood for different Fs

non_cucks_female.adulthood.indF.plot <- ggplot(non_cucks_female_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.adulthood.mumF.plot <- ggplot(non_cucks_female_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.adulthood.dadF.plot <- ggplot(non_cucks_female_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to Adulthood", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.indF.plot <- ggplot(male_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.mumF.plot <- ggplot(male_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.dadF.plot <- ggplot(male_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to Adulthood", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


adulthood.global <- (non_cucks_female.adulthood.indF.plot + male.adulthood.indF.plot) / (non_cucks_female.adulthood.mumF.plot + male.adulthood.mumF.plot) / (non_cucks_female.adulthood.dadF.plot + male.adulthood.dadF.plot)

ggsave("Figures/adulthood_global.png", adulthood.global)

#----------------non_cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
non_cucks_female.adulthood.2 <- update(non_cucks_female.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                         optCtrl = list(maxfun = 100000000)))

summary(non_cucks_female.adulthood.2) #interaction term not significant.
vif(non_cucks_female.adulthood.2)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.2) #models not significantly different


#check dadF help interaction. 
non_cucks_female.adulthood.3 <- update(non_cucks_female.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.3) #interaction term very almost significant p = 0.0503 .
vif(non_cucks_female.adulthood.3)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_female.adulthood.4 <- update(non_cucks_female.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.4) #interaction term not sig.
vif(non_cucks_female.adulthood.3)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.3)  #models are significantly different

#----------------------------------non_cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
non_cucks_female.adulthood.4 <- update(non_cucks_female.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.4) #interaction term not significant.
vif(non_cucks_female.adulthood.2)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_female.adulthood.5 <- update(non_cucks_female.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.5) #interaction term not sig .
vif(non_cucks_female.adulthood.3)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_female.adulthood.6 <- update(non_cucks_female.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.6) #interaction term not sig.
vif(non_cucks_female.adulthood.3)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select non bastard
non_cucks_male_birds <-filter(male_birds, EPP == 0) 

#basic model
non_cucks_male.adulthood.1 <- glmer(adulthood ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                          data = non_cucks_male_birds,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_male.adulthood.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_male.adulthood.1)



#----------plot sig---------
plot1 <- ggpredict(non_cucks_male.adulthood.1, terms = "dadF [all]")

plot(plot1,show_data = TRUE ) +
 theme_bw()+




#-----------------non_cucks_male helper interactions-----------





#check mumF help interaction. 
non_cucks_male.adulthood.2 <- update(non_cucks_male.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.2) #interaction term not significant.
vif(non_cucks_male.adulthood.2)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.2) #models not significantly different


#check dadF help interaction. 
non_cucks_male.adulthood.3 <- update(non_cucks_male.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.3) #interaction term not significant
vif(non_cucks_male.adulthood.3)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_male.adulthood.4 <- update(non_cucks_male.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.4) #interaction term not sig. Nearly though
vif(non_cucks_male.adulthood.3)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.3)  #models are not significantly different



#-------------------non_cucks_male mumage interaction-------


#check mumF MumAge interaction. 
non_cucks_male.adulthood.4 <- update(non_cucks_male.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.4) #interaction term not significant.
vif(non_cucks_male.adulthood.2)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_male.adulthood.5 <- update(non_cucks_male.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.5) #interaction term not sig .
vif(non_cucks_male.adulthood.3)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_male.adulthood.6 <- update(non_cucks_male.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.6) #interaction term not sig.
vif(non_cucks_male.adulthood.3)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.6)  #models are not significantly different



##-------------------cuck tests--------------------------

#select cuck offspring
cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for cucks_female  birds------------
cucks_female.adulthood.1 <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                      data = cucks_female_birds,
                                      family=binomial,
                                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(cucks_female.adulthood.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_female.adulthood.1, plot = T)

#check for collinearirty with vif
vif(cucks_female.adulthood.1)

##plot basic model with no interactions
#plot recruitment to adulthood for different Fs

cucks_female.adulthood.indF.plot <- ggplot(cucks_female_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.adulthood.mumF.plot <- ggplot(cucks_female_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.adulthood.dadF.plot <- ggplot(cucks_female_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to Adulthood", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.indF.plot <- ggplot(male_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.mumF.plot <- ggplot(male_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.dadF.plot <- ggplot(male_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to Adulthood", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


adulthood.global <- (cucks_female.adulthood.indF.plot + male.adulthood.indF.plot) / (cucks_female.adulthood.mumF.plot + male.adulthood.mumF.plot) / (cucks_female.adulthood.dadF.plot + male.adulthood.dadF.plot)

ggsave("Figures/adulthood_global.png", adulthood.global)

#----------------cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
cucks_female.adulthood.2 <- update(cucks_female.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                         optCtrl = list(maxfun = 100000000)))

summary(cucks_female.adulthood.2) #interaction term not significant.
vif(cucks_female.adulthood.2)
anova(cucks_female.adulthood.1,cucks_female.adulthood.2) #models not significantly different


#check dadF help interaction. 
cucks_female.adulthood.3 <- update(cucks_female.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.3) #interaction term very almost significant p = 0.0503 .
vif(cucks_female.adulthood.3)
anova(cucks_female.adulthood.1,cucks_female.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
cucks_female.adulthood.4 <- update(cucks_female.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.4) #interaction term not sig.
vif(cucks_female.adulthood.3)
anova(cucks_female.adulthood.1,cucks_female.adulthood.3)  #models are significantly different

#----------------------------------cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
cucks_female.adulthood.4 <- update(cucks_female.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.4) #interaction term not significant.
vif(cucks_female.adulthood.2)
anova(cucks_female.adulthood.1,cucks_female.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_female.adulthood.5 <- update(cucks_female.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.5) #interaction term not sig .
vif(cucks_female.adulthood.3)
anova(cucks_female.adulthood.1,cucks_female.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_female.adulthood.6 <- update(cucks_female.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.6) #interaction term not sig.
vif(cucks_female.adulthood.3)
anova(cucks_female.adulthood.1,cucks_female.adulthood.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select bastard
cucks_male_birds <-filter(na.omit(male_birds, EPP == 1) )

#basic model
cucks_male.adulthood.1 <- glmer(adulthood ~ mumF + dadF +FROH + BrMF+  help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                                    data = cucks_male_birds,
                                    family=binomial,
                                    control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_male.adulthood.1, plot = T)

#check for collinearirty with vif
vif(cucks_male.adulthood.1)

plot2 <- ggpredict(cucks_male.adulthood.1, terms = "dadF [all]")

plot(plot2,show_data = TRUE ) +
  theme_bw()

plot3 <- ggpredict(cucks_male.adulthood.1, terms = "BrMF [all]")

plot(plot3,show_data = TRUE ) +
  theme_bw()


#-----------------cucks_male helper interactions-----------





#check mumF help interaction. 
cucks_male.adulthood.2 <- update(cucks_male.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.2) #interaction term not significant.
vif(cucks_male.adulthood.2)
anova(cucks_male.adulthood.1,cucks_male.adulthood.2) #models not significantly different


#check dadF help interaction. 
cucks_male.adulthood.3 <- update(cucks_male.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.3) #interaction term not significant
vif(cucks_male.adulthood.3)
anova(cucks_male.adulthood.1,cucks_male.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
cucks_male.adulthood.4 <- update(cucks_male.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.4) #interaction term not sig. Nearly though
vif(cucks_male.adulthood.3)
anova(cucks_male.adulthood.1,cucks_male.adulthood.3)  #models are not significantly different


#-------------------cucks_male mumage interaction-------

#check mumF MumAge interaction. 
cucks_male.adulthood.4 <- update(cucks_male.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.4) #interaction term not significant.
vif(cucks_male.adulthood.2)
anova(cucks_male.adulthood.1,cucks_male.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_male.adulthood.5 <- update(cucks_male.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.5) #interaction term not sig .
vif(cucks_male.adulthood.3)
anova(cucks_male.adulthood.1,cucks_male.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_male.adulthood.6 <- update(cucks_male.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.6) #interaction term not sig.
vif(cucks_male.adulthood.3)
anova(cucks_male.adulthood.1,cucks_male.adulthood.6)  #models are not significantly different



##-----------------coverage subsample--------------


# select males and females, birds that are already dead and never translocated, and coverage above 2X

cov_female_birds <- stat.birds %>% filter(Sex == 0, event == 0, Coverage > 2, CovMum > 2, CovDad > 2)

cov_male_birds <- stat.birds %>% filter(Sex == 1, event == 0, Coverage > 2, CovMum > 2, CovDad > 2)

#--------------------------juvenile survival---------------

#test if maternal or paternal roh affects the juvenile survival
#do this using 'adulthood' variable, which is 0 for bird died before 1 year old
#this is a success or failure, so should use binomial distribution

#-------------------------tests for cov_female birds------------
cov_female.adulthood.1 <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                            data = cov_female_birds,
                            family=binomial,
                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(cov_female.adulthood.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cov_female.adulthood.1, plot = T)

#check for collinearirty with vif
vif(cov_female.adulthood.1)

##plot basic model with no interactions
#plot recruitment to adulthood for different Fs

female.adulthood.indF.plot <- ggplot(female_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.adulthood.mumF.plot <- ggplot(female_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.adulthood.dadF.plot <- ggplot(female_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to Adulthood", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.indF.plot <- ggplot(male_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.mumF.plot <- ggplot(male_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.dadF.plot <- ggplot(male_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to Adulthood", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


adulthood.global <- (female.adulthood.indF.plot + male.adulthood.indF.plot) / (female.adulthood.mumF.plot + male.adulthood.mumF.plot) / (female.adulthood.dadF.plot + male.adulthood.dadF.plot)

ggsave("Figures/adulthood_global.png", adulthood.global)

#----------------female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
female.adulthood.2 <- update(female.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.2) #interaction term not significant.
vif(female.adulthood.2)
anova(female.adulthood.1,female.adulthood.2) #models not significantly different


#check dadF help interaction. 
female.adulthood.3 <- update(female.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.3) #interaction term very almost significant p = 0.0503 .
vif(female.adulthood.3)
anova(female.adulthood.1,female.adulthood.3)  #models are not significantly different

#visually asses the difference.
#It seems that having helpers means that birds are LESS likely to recruit to adulthood when dads are inbred
female.adulthood.dadF.helpYES.plot <- ggplot(female_birds[female_birds$help == 1,], aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  scale_y_continuous(limits = c(0, 0.25))+
  labs( x = "Recruitment to Adulthood", y = "Paternal F", title = "Helpers in Nest")+
  theme_bw()+
  theme(text = element_text(size = 20))

female.adulthood.dadF.helpNO.plot <- ggplot(female_birds[female_birds$help == 0,], aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  scale_y_continuous(limits = c(0, 0.25))+
  labs( x = "Recruitment to Adulthood", y = "Paternal F", title = "No Helpers")+
  theme_bw()+
  theme(text = element_text(size = 20))

adulthood.female.help <- female.adulthood.dadF.helpYES.plot | female.adulthood.dadF.helpNO.plot 

ggsave("Figures/adulthood_female_help.png", adulthood.female.help)

#check FROH help interaction. 
female.adulthood.4 <- update(female.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.4) #interaction term not sig.
vif(female.adulthood.3)
anova(female.adulthood.1,female.adulthood.3)  #models are significantly different

#----------------------------------female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
female.adulthood.4 <- update(female.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.4) #interaction term not significant.
vif(female.adulthood.2)
anova(female.adulthood.1,female.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
female.adulthood.5 <- update(female.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.5) #interaction term not sig .
vif(female.adulthood.3)
anova(female.adulthood.1,female.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
female.adulthood.6 <- update(female.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.adulthood.6) #interaction term not sig.
vif(female.adulthood.3)
anova(female.adulthood.1,female.adulthood.6)  #models are not significantly different

#-------------------------tests for male birds------------
#basic model
male.adulthood.1 <- glmer(adulthood ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                          data = male_birds,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.1) #here FROH and mumage arte sig

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = male.adulthood.1, plot = T)

#check for collinearirty with vif
vif(male.adulthood.1)

#-----------------male helper interactions-----------

#check mumF help interaction. 
male.adulthood.2 <- update(male.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.2) #interaction term not significant.
vif(male.adulthood.2)
anova(male.adulthood.1,male.adulthood.2) #models not significantly different


#check dadF help interaction. 
male.adulthood.3 <- update(male.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.3) #interaction term not significant
vif(male.adulthood.3)
anova(male.adulthood.1,male.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
male.adulthood.4 <- update(male.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.4) #interaction term not sig. Nearly though
vif(male.adulthood.3)
anova(male.adulthood.1,male.adulthood.3)  #models are not significantly different


#-------------------male mumage interaction-------

#check mumF MumAge interaction. 
male.adulthood.4 <- update(male.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.4) #interaction term not significant.
vif(male.adulthood.2)
anova(male.adulthood.1,male.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
male.adulthood.5 <- update(male.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.5) #interaction term not sig .
vif(male.adulthood.3)
anova(male.adulthood.1,male.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
male.adulthood.6 <- update(male.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.adulthood.6) #interaction term not sig.
vif(male.adulthood.3)
anova(male.adulthood.1,male.adulthood.6)  #models are not significantly different


##--------------------------------EPP tests-------------------------------------


#select non cuck offspring
non_cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for non_cucks_female  birds------------
non_cucks_female.adulthood.1 <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                      data = non_cucks_female_birds,
                                      family=binomial,
                                      control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(non_cucks_female.adulthood.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_female.adulthood.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_female.adulthood.1)

##plot basic model with no interactions
#plot recruitment to adulthood for different Fs

non_cucks_female.adulthood.indF.plot <- ggplot(non_cucks_female_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.adulthood.mumF.plot <- ggplot(non_cucks_female_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.adulthood.dadF.plot <- ggplot(non_cucks_female_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to Adulthood", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.indF.plot <- ggplot(male_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.mumF.plot <- ggplot(male_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.dadF.plot <- ggplot(male_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to Adulthood", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


adulthood.global <- (non_cucks_female.adulthood.indF.plot + male.adulthood.indF.plot) / (non_cucks_female.adulthood.mumF.plot + male.adulthood.mumF.plot) / (non_cucks_female.adulthood.dadF.plot + male.adulthood.dadF.plot)

ggsave("Figures/adulthood_global.png", adulthood.global)

#----------------non_cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
non_cucks_female.adulthood.2 <- update(non_cucks_female.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                         optCtrl = list(maxfun = 100000000)))

summary(non_cucks_female.adulthood.2) #interaction term not significant.
vif(non_cucks_female.adulthood.2)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.2) #models not significantly different


#check dadF help interaction. 
non_cucks_female.adulthood.3 <- update(non_cucks_female.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.3) #interaction term very almost significant p = 0.0503 .
vif(non_cucks_female.adulthood.3)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_female.adulthood.4 <- update(non_cucks_female.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.4) #interaction term not sig.
vif(non_cucks_female.adulthood.3)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.3)  #models are significantly different

#----------------------------------non_cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
non_cucks_female.adulthood.4 <- update(non_cucks_female.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.4) #interaction term not significant.
vif(non_cucks_female.adulthood.2)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_female.adulthood.5 <- update(non_cucks_female.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.5) #interaction term not sig .
vif(non_cucks_female.adulthood.3)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_female.adulthood.6 <- update(non_cucks_female.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.adulthood.6) #interaction term not sig.
vif(non_cucks_female.adulthood.3)
anova(non_cucks_female.adulthood.1,non_cucks_female.adulthood.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select non bastard
non_cucks_male_birds <-filter(male_birds, EPP == 0) 

#basic model
non_cucks_male.adulthood.1 <- glmer(adulthood ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                                    data = non_cucks_male_birds,
                                    family=binomial,
                                    control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_male.adulthood.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_male.adulthood.1)


#-----------------non_cucks_male helper interactions-----------





#check mumF help interaction. 
non_cucks_male.adulthood.2 <- update(non_cucks_male.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.2) #interaction term not significant.
vif(non_cucks_male.adulthood.2)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.2) #models not significantly different


#check dadF help interaction. 
non_cucks_male.adulthood.3 <- update(non_cucks_male.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.3) #interaction term not significant
vif(non_cucks_male.adulthood.3)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_male.adulthood.4 <- update(non_cucks_male.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.4) #interaction term not sig. Nearly though
vif(non_cucks_male.adulthood.3)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.3)  #models are not significantly different



#-------------------non_cucks_male mumage interaction-------


#check mumF MumAge interaction. 
non_cucks_male.adulthood.4 <- update(non_cucks_male.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.4) #interaction term not significant.
vif(non_cucks_male.adulthood.2)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_male.adulthood.5 <- update(non_cucks_male.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.5) #interaction term not sig .
vif(non_cucks_male.adulthood.3)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_male.adulthood.6 <- update(non_cucks_male.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.adulthood.6) #interaction term not sig.
vif(non_cucks_male.adulthood.3)
anova(non_cucks_male.adulthood.1,non_cucks_male.adulthood.6)  #models are not significantly different



##-------------------cuck tests--------------------------

#select cuck offspring
cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for cucks_female  birds------------
cucks_female.adulthood.1 <- glmer(adulthood ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                  data = cucks_female_birds,
                                  family=binomial,
                                  control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(cucks_female.adulthood.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_female.adulthood.1, plot = T)

#check for collinearirty with vif
vif(cucks_female.adulthood.1)

##plot basic model with no interactions
#plot recruitment to adulthood for different Fs

cucks_female.adulthood.indF.plot <- ggplot(cucks_female_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.adulthood.mumF.plot <- ggplot(cucks_female_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.adulthood.dadF.plot <- ggplot(cucks_female_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to Adulthood", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.indF.plot <- ggplot(male_birds, aes(adulthood, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.mumF.plot <- ggplot(male_birds, aes(adulthood,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.adulthood.dadF.plot <- ggplot(male_birds, aes(adulthood,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to Adulthood", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


adulthood.global <- (cucks_female.adulthood.indF.plot + male.adulthood.indF.plot) / (cucks_female.adulthood.mumF.plot + male.adulthood.mumF.plot) / (cucks_female.adulthood.dadF.plot + male.adulthood.dadF.plot)

ggsave("Figures/adulthood_global.png", adulthood.global)

#----------------cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
cucks_female.adulthood.2 <- update(cucks_female.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                 optCtrl = list(maxfun = 100000000)))

summary(cucks_female.adulthood.2) #interaction term not significant.
vif(cucks_female.adulthood.2)
anova(cucks_female.adulthood.1,cucks_female.adulthood.2) #models not significantly different


#check dadF help interaction. 
cucks_female.adulthood.3 <- update(cucks_female.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.3) #interaction term very almost significant p = 0.0503 .
vif(cucks_female.adulthood.3)
anova(cucks_female.adulthood.1,cucks_female.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
cucks_female.adulthood.4 <- update(cucks_female.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.4) #interaction term not sig.
vif(cucks_female.adulthood.3)
anova(cucks_female.adulthood.1,cucks_female.adulthood.3)  #models are significantly different

#----------------------------------cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
cucks_female.adulthood.4 <- update(cucks_female.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.4) #interaction term not significant.
vif(cucks_female.adulthood.2)
anova(cucks_female.adulthood.1,cucks_female.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_female.adulthood.5 <- update(cucks_female.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.5) #interaction term not sig .
vif(cucks_female.adulthood.3)
anova(cucks_female.adulthood.1,cucks_female.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_female.adulthood.6 <- update(cucks_female.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.adulthood.6) #interaction term not sig.
vif(cucks_female.adulthood.3)
anova(cucks_female.adulthood.1,cucks_female.adulthood.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select bastard
cucks_male_birds <-filter(na.omit(male_birds, EPP == 0) )

#basic model
cucks_male.adulthood.1 <- glmer(adulthood ~ mumF + dadF +FROH + BrMF+  help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                                data = cucks_male_birds,
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_male.adulthood.1, plot = T)

#check for collinearirty with vif
vif(cucks_male.adulthood.1)


#-----------------cucks_male helper interactions-----------





#check mumF help interaction. 
cucks_male.adulthood.2 <- update(cucks_male.adulthood.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.2) #interaction term not significant.
vif(cucks_male.adulthood.2)
anova(cucks_male.adulthood.1,cucks_male.adulthood.2) #models not significantly different


#check dadF help interaction. 
cucks_male.adulthood.3 <- update(cucks_male.adulthood.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.3) #interaction term not significant
vif(cucks_male.adulthood.3)
anova(cucks_male.adulthood.1,cucks_male.adulthood.3)  #models are not significantly different


#check FROH help interaction. 
cucks_male.adulthood.4 <- update(cucks_male.adulthood.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.4) #interaction term not sig. Nearly though
vif(cucks_male.adulthood.3)
anova(cucks_male.adulthood.1,cucks_male.adulthood.3)  #models are not significantly different


#-------------------cucks_male mumage interaction-------

#check mumF MumAge interaction. 
cucks_male.adulthood.4 <- update(cucks_male.adulthood.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.4) #interaction term not significant.
vif(cucks_male.adulthood.2)
anova(cucks_male.adulthood.1,cucks_male.adulthood.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_male.adulthood.5 <- update(cucks_male.adulthood.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.5) #interaction term not sig .
vif(cucks_male.adulthood.3)
anova(cucks_male.adulthood.1,cucks_male.adulthood.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_male.adulthood.6 <- update(cucks_male.adulthood.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.adulthood.6) #interaction term not sig.
vif(cucks_male.adulthood.3)
anova(cucks_male.adulthood.1,cucks_male.adulthood.6)  #models are not significantly different





