#---------------------Chapter 2.4 Lifetime reproductive success-----------------
#Do a GLMM on parental Fs and Lifetime reproductive success(LRS).
#LRS can be measured with the total number of offspring that survived > 1 year

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
stat.birds$BreedSeason <- as.factor(stat.birds$BreedSeason)
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$DeathEventL <- as.logical(stat.birds$DeathEvent)
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


both_birds <- stat.birds %>% filter(event == 0,
                                    Coverage >= 0.25 ,
                                    CovMum >= 0.25 ,
                                    CovDad >= 0.25,
                                    CovBrM >= 0.25)



#--------------------------lrs---------------

#test if maternal or paternal or social roh affects the lrs
#do this using 'n_off' variable, which is number of offspring that survived >1 year

hist(female_birds$n_off) #data is very 0 inflated


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

simulateResiduals(both.lrs.poisson, plot = T)  #looks pretty ok?

#check for collinearity and overdispersion
check_collinearity(both.lrs.poisson) #vif around one
check_overdispersion(both.lrs.poisson) #some  overdisp



#check nbinom1 family distribution

both.lrs.nbinom1 <- update(both.lrs.poisson, family = nbinom1) #convergance problems

simulateResiduals(both.lrs.nbinom1, plot = T) #looks good

#check for collinearity and overdispersion
check_collinearity(both.lrs.nbinom1) #vif around one
check_overdispersion(both.lrs.nbinom1) #no overdisp

summary(both.lrs.nbinom1)


#check nbinom2

both.lrs.nbinom2 <- update(both.lrs.poisson, family = nbinom2) #more convergance problems

simulateResiduals(both.lrs.nbinom2, plot = T) #looks good

#check for collinearity and overdispersion
check_collinearity(both.lrs.nbinom2) #vif around one
check_overdispersion(both.lrs.nbinom2) #no overdisp

summary(both.lrs.nbinom2)


#AIV comparison
AIC(both.lrs.nbinom1,both.lrs.nbinom2,both.lrs.poisson)
#both.lrs.nbinom1 31 1314.986
#both.lrs.nbinom2 31 1307.910
#both.lrs.poisson 30 1319.289

#go with poisson, aic problems in others.




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


#in summary, mumFROH AND SEX EFFECT IS ALL

#--------output---


#plot for eseb poster

poster_data <- select(both_birds ,n_off,Father = dadLargeFroh,Mother = mumLargeFroh)

long_poster <- pivot_longer(poster_data,!n_off, names_to = "Parent", values_to = "FROH") %>%
  mutate(n_off = as.numeric(n_off), )

long_poster

poster.inter <-  ggplot(long_poster,aes(x=FROH, y = n_off, colour = Parent))+
  geom_smooth(method = "lm")+
    geom_point()+
    theme_classic2()+
    theme(text = element_text(size = 40),
          plot.title = element_text(hjust = 0.2))+
    ggtitle("Intergenerational inbreeding depression")+
  labs( x = "Parental FROH > 3.3Mb", y= "Lifetime reproductive success")

poster.inter 

  
poster.int.plot <- ggplot(both_birds, aes(y = n_off)) +
  geom_smooth(aes(x = dadLargeFroh), method = "lm" , colour = "red", size = 2,fill = "red")+
  theme_classic2()+
  labs( x = "FROH > 3.3Mb", y= "Lifetime reproductive success")+
  theme(text = element_text(size = 40),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Individual inbreeding depression")



  
  #output mumfroh x sex interaction
  
summary(both.lrs.poisson.c)

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










#plot brm*rain variance

interact_plot(both.lrs.poisson.int.7,
              pred = BrMLargeFroh,
              modx = mean_rain_cv,
              modx.values = "terciles",
              int.type = "confidence",
              plot.points = T,
              data = plot_birds)+
  theme_classic2()+
  labs( x = "Social Father FROH",
        y= "LRS",
        fill = "Mean annual variance in rainfall over lifespan")+
  theme(text = element_text(size = 20)) +
  ggtitle("Fig. 7.1 - Predicted fit of GLMM LRS ~ Social Father FROH
              with Mean Annual Variance in rainfall interaction")


#plot mother froh Sex

ggplot(both_birds,aes(mumLargeFroh, n_off , group = Sex, fill = Sex))+
  geom_smooth(method = "lm")+
  geom_jitter(aes(color  = Sex))+
  theme_classic2()+
 labs( x = bquote("Mother's F"[ROH]),
        y= "Productivity" ) +
  theme(text = element_text(size = 20)) +
  ggtitle("")


#plot mother froh helper

ggplot(both_birds,aes(mumLargeFroh, n_off , group = helpF, fill = helpF))+
  geom_smooth(method = "lm")+
  geom_jitter(aes(color  = helpF))+
  theme_classic2()+
  labs( x = "Mother's FROH > 3.3Mb",
        y= "LRS",)+
  theme(text = element_text(size = 20)) +
  ggtitle("Fig. 7.3 Mother's FROH and LRS by helper presence")


#output this plot fora poster



LRS.parental.model <- both.lrs.poisson <- glmmTMB(n_off ~ LargeFROH + dadLargeFroh + mumLargeFroh + BrMLargeFroh +
                                                help + Sex + 
                                                mean_total_rain + mean_rain_cv + birth_year +
                                                (1 | mum) + (1 | dad),
                                              data = filter(both_birds, n_off < 16),
                                              ziformula=~1,
                                              family=poisson)


#for mothers

LRS.maternal.model.ggpred <-  ggpredict(LRS.parental.model,
                              terms = "mumLargeFroh[all]")


autoplot(LRS.maternal.model.ggpred)+
  geom_expected_line(linewidth = 1, color = "red")+
  geom_CI_ribbon(fill = "red")+
  theme_classic2()+
  labs( x = "          FROH > 3.3Mb", y= "Lifetime reproductive success")+
  theme(text = element_text(size = 50)) +
  ggtitle("     Intergenerational inbreeding")+
  layer_fit_data(color = "red")


#for social fathers

LRS.parental.model  <- glmmTMB(n_off ~ LargeFROH + dadLargeFroh + mumLargeFroh + BrMLargeFroh +
                                                    help + Sex + 
                                                    mean_total_rain + mean_rain_cv + birth_year +
                                                    (1 | mum) + (1 | dad),
                                                  data = filter(both_birds, n_off < 16, BrMLargeFroh < 0.45),
                                                  ziformula=~1,
                               family=poisson)

LRS.paternal.model.ggpred <-  ggpredict(LRS.parental.model,
                                        terms = "BrMLargeFroh[all]")


autoplot(LRS.paternal.model.ggpred )+
  geom_expected_line(linewidth = 1, color = "blue")+
  geom_CI_ribbon(fill = "blue")+
  theme_classic2()+
  labs( x = "Father's FROH > 3.3Mb", y= "Lifetime reproductive success")+
  theme(text = element_text(size = 40)) +
  ggtitle("         Paternal inbreeding")+
  layer_fit_data(alpha = 0.2,color = "blue")







#------------test lrs split by sex-------



female.lrs.poisson.1 <- glmmTMB(n_off ~ rescale(LargeFROH) + 
                                          rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                          rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                          sib_pres+
                                          rescale(birth_year) + 
                                          lifespan  ,
                                        data = female_birds,
                                        ziformula=~rescale(LargeFROH) + 
                                          rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                          rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                          sib_pres+
                                          rescale(birth_year) + 
                                          lifespan,
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



male.lrs.poisson.1 <- glmmTMB(n_off ~ rescale(LargeFROH) + 
                                  rescale(BrMLargeFroh) + rescale(dadLargeFroh) + rescale(mumLargeFroh) +
                                  rescale(mean_total_rain) + rescale(mean_rain_cv) +
                                  sib_pres+
                                  rescale(birth_year) + 
                                  lifespan  ,
                                data = female_birds,
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



#----------------epp split test----------
#lets split these into EPP and non EPP. Just for fun.


wp_birds <- filter(stat.birds, EPP == 0)
dim(wp_birds) #n = 1038

epp_birds <- filter(stat.birds, EPP == 1)
dim(epp_birds) #n = 356


wp.lrs.poisson.1 <- glmmTMB(n_off ~ LargeFROH + dadLargeFroh + mumLargeFroh + 
                              help + Sex + 
                              mean_total_rain + mean_rain_cv + birth_year +
                              (1 | mum) + (1 | dad),
                            data = both_birds,
                            ziformula=~1,
                            family=poisson)

summary(wp.lrs.poisson.1) #the interaction is significant


#and now extra pair!

#without int
epp.lrs.poisson.0 <- glmmTMB(n_off ~ LargeFROH + dadLargeFroh + mumLargeFroh + BrMLargeFroh +
                               h_countF + Sex + 
                               mean_total_rain + mean_rain_cv + birth_year +
                               (1 | mum) + (1 | dad),
                             data = both_birds,
                             ziformula=~1,
                             family=poisson)

summary(epp.lrs.poisson.0)
check_collinearity(epp.lrs.poisson.0)

#test brm int first

epp.lrs.poisson.1 <- glmmTMB(n_off ~ LargeFROH + dadLargeFroh + mumLargeFroh + BrMLargeFroh +
                              help + Sex + 
                              mean_total_rain + mean_rain_cv + birth_year +
                              BrMLargeFroh*mean_rain_cv +
                              (1 | mum) + (1 | dad),
                            data = both_birds,
                            ziformula=~1,
                            family=poisson)

summary(epp.lrs.poisson.1) #the interaction is significant

#and test gen dad int

epp.lrs.poisson.2 <- glmmTMB(n_off ~ LargeFROH + dadLargeFroh + mumLargeFroh + BrMLargeFroh +
                               help + Sex + 
                               mean_total_rain + mean_rain_cv + birth_year +
                               dadLargeFroh*mean_rain_cv +
                               (1 | mum) + (1 | dad),
                             data = both_birds,
                             ziformula=~1,
                             family=poisson)

summary(epp.lrs.poisson.2) #the interaction is also significant. What?


epp.lrs.poisson.3 <- glmmTMB(n_off ~ LargeFROH + dadLargeFroh + mumLargeFroh + BrMLargeFroh +
                               help + Sex + 
                               mean_total_rain + mean_rain_cv + birth_year +
                               dadLargeFroh*mean_rain_cv +
                               BrMLargeFroh*mean_rain_cv +
                               (1 | mum) + (1 | dad),
                             data = both_birds,
                             ziformula=~1,
                             family=poisson)

summary(epp.lrs.poisson.3) #neither are sig





#


































































#----------joint sex models


both_birds <- stat.birds %>% filter(event == 0, Coverage >= 1 , CovMum >= 1 , CovDad >= 1, CovBrM >= 1, !is.na(BirthInsect), !is.na(meanInsects))

chuen_birds <- both_birds %>% filter(first_catch > 0)

dim(both_birds) #534 samples

both.lrs.1 <- glmer(n_off ~ mumF + dadF + FROH + help + MumAge + Sex + meanInsects +(1 | birth_year) + (1 | mum) + (1 | dad),
                           data = both_birds,
                          family=poisson)

simulateResiduals(both.lrs.1, plot = T)
vif(both.adulthood.1)

summary(both.lrs.1)  #nothing significant when dataset is not split by EPP




#-------------------------tests for female birds------------
female.n_off.1 <- glmer.nb(n_off ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                            data = female_birds)

summary(female.n_off.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = female.n_off.1, plot = T)

#check for collinearirty with vif
vif(female.n_off.1)

##plot basic model with no interactions
#plot recruitment to n_off for different Fs

female.n_off.indF.plot <- ggplot(female_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.n_off.mumF.plot <- ggplot(female_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.n_off.dadF.plot <- ggplot(female_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to n_off", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.indF.plot <- ggplot(male_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.mumF.plot <- ggplot(male_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.dadF.plot <- ggplot(male_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to n_off", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


n_off.global <- (female.n_off.indF.plot + male.n_off.indF.plot) / (female.n_off.mumF.plot + male.n_off.mumF.plot) / (female.n_off.dadF.plot + male.n_off.dadF.plot)

ggsave("Figures/n_off_global.png", n_off.global)

#----------------female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
female.n_off.2 <- update(female.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.2) #interaction term not significant.
vif(female.n_off.2)
anova(female.n_off.1,female.n_off.2) #models not significantly different


#check dadF help interaction. 
female.n_off.3 <- update(female.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.3) #interaction term very almost significant p = 0.0503 .
vif(female.n_off.3)
anova(female.n_off.1,female.n_off.3)  #models are not significantly different

#visually asses the difference.
#It seems that having helpers means that birds are LESS likely to recruit to n_off when dads are inbred
female.n_off.dadF.helpYES.plot <- ggplot(female_birds[female_birds$help == 1,], aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  scale_y_continuous(limits = c(0, 0.25))+
  labs( x = "Recruitment to n_off", y = "Paternal F", title = "Helpers in Nest")+
  theme_bw()+
  theme(text = element_text(size = 20))

female.n_off.dadF.helpNO.plot <- ggplot(female_birds[female_birds$help == 0,], aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  scale_y_continuous(limits = c(0, 0.25))+
  labs( x = "Recruitment to n_off", y = "Paternal F", title = "No Helpers")+
  theme_bw()+
  theme(text = element_text(size = 20))

n_off.female.help <- female.n_off.dadF.helpYES.plot | female.n_off.dadF.helpNO.plot 

ggsave("Figures/n_off_female_help.png", n_off.female.help)

#check FROH help interaction. 
female.n_off.4 <- update(female.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.4) #interaction term not sig.
vif(female.n_off.3)
anova(female.n_off.1,female.n_off.3)  #models are significantly different

#----------------------------------female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
female.n_off.4 <- update(female.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.4) #interaction term not significant.
vif(female.n_off.2)
anova(female.n_off.1,female.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
female.n_off.5 <- update(female.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.5) #interaction term not sig .
vif(female.n_off.3)
anova(female.n_off.1,female.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
female.n_off.6 <- update(female.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.6) #interaction term not sig.
vif(female.n_off.3)
anova(female.n_off.1,female.n_off.6)  #models are not significantly different

#-------------------------tests for male birds------------
#basic model
male.n_off.1 <- glmer(n_off ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                          data = male_birds,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.1) #here FROH and mumage arte sig

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = male.n_off.1, plot = T)

#check for collinearirty with vif
vif(male.n_off.1)

#-----------------male helper interactions-----------

#check mumF help interaction. 
male.n_off.2 <- update(male.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.2) #interaction term not significant.
vif(male.n_off.2)
anova(male.n_off.1,male.n_off.2) #models not significantly different


#check dadF help interaction. 
male.n_off.3 <- update(male.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.3) #interaction term not significant
vif(male.n_off.3)
anova(male.n_off.1,male.n_off.3)  #models are not significantly different


#check FROH help interaction. 
male.n_off.4 <- update(male.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.4) #interaction term not sig. Nearly though
vif(male.n_off.3)
anova(male.n_off.1,male.n_off.3)  #models are not significantly different


#-------------------male mumage interaction-------

#check mumF MumAge interaction. 
male.n_off.4 <- update(male.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.4) #interaction term not significant.
vif(male.n_off.2)
anova(male.n_off.1,male.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
male.n_off.5 <- update(male.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.5) #interaction term not sig .
vif(male.n_off.3)
anova(male.n_off.1,male.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
male.n_off.6 <- update(male.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.6) #interaction term not sig.
vif(male.n_off.3)
anova(male.n_off.1,male.n_off.6)  #models are not significantly different


##--------------------------------EPP tests-------------------------------------


#select non cuck offspring
non_cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for non_cucks_female  birds------------
non_cucks_female.n_off.1 <- glmer(n_off ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                      data = non_cucks_female_birds,
                                      family=binomial,
                                      control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(non_cucks_female.n_off.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_female.n_off.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_female.n_off.1)

##plot basic model with no interactions
#plot recruitment to n_off for different Fs

non_cucks_female.n_off.indF.plot <- ggplot(non_cucks_female_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.n_off.mumF.plot <- ggplot(non_cucks_female_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.n_off.dadF.plot <- ggplot(non_cucks_female_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to n_off", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.indF.plot <- ggplot(male_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.mumF.plot <- ggplot(male_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.dadF.plot <- ggplot(male_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to n_off", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


n_off.global <- (non_cucks_female.n_off.indF.plot + male.n_off.indF.plot) / (non_cucks_female.n_off.mumF.plot + male.n_off.mumF.plot) / (non_cucks_female.n_off.dadF.plot + male.n_off.dadF.plot)

ggsave("Figures/n_off_global.png", n_off.global)

#----------------non_cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
non_cucks_female.n_off.2 <- update(non_cucks_female.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                         optCtrl = list(maxfun = 100000000)))

summary(non_cucks_female.n_off.2) #interaction term not significant.
vif(non_cucks_female.n_off.2)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.2) #models not significantly different


#check dadF help interaction. 
non_cucks_female.n_off.3 <- update(non_cucks_female.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.3) #interaction term very almost significant p = 0.0503 .
vif(non_cucks_female.n_off.3)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_female.n_off.4 <- update(non_cucks_female.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.4) #interaction term not sig.
vif(non_cucks_female.n_off.3)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.3)  #models are significantly different

#----------------------------------non_cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
non_cucks_female.n_off.4 <- update(non_cucks_female.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.4) #interaction term not significant.
vif(non_cucks_female.n_off.2)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_female.n_off.5 <- update(non_cucks_female.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.5) #interaction term not sig .
vif(non_cucks_female.n_off.3)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_female.n_off.6 <- update(non_cucks_female.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.6) #interaction term not sig.
vif(non_cucks_female.n_off.3)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select non bastard
non_cucks_male_birds <-filter(male_birds, EPP == 0) 

#basic model
non_cucks_male.n_off.1 <- glmer(n_off ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                                    data = non_cucks_male_birds,
                                    family=binomial,
                                    control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_male.n_off.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_male.n_off.1)


#-----------------non_cucks_male helper interactions-----------





#check mumF help interaction. 
non_cucks_male.n_off.2 <- update(non_cucks_male.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.2) #interaction term not significant.
vif(non_cucks_male.n_off.2)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.2) #models not significantly different


#check dadF help interaction. 
non_cucks_male.n_off.3 <- update(non_cucks_male.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.3) #interaction term not significant
vif(non_cucks_male.n_off.3)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_male.n_off.4 <- update(non_cucks_male.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.4) #interaction term not sig. Nearly though
vif(non_cucks_male.n_off.3)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.3)  #models are not significantly different



#-------------------non_cucks_male mumage interaction-------


#check mumF MumAge interaction. 
non_cucks_male.n_off.4 <- update(non_cucks_male.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.4) #interaction term not significant.
vif(non_cucks_male.n_off.2)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_male.n_off.5 <- update(non_cucks_male.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.5) #interaction term not sig .
vif(non_cucks_male.n_off.3)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_male.n_off.6 <- update(non_cucks_male.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.6) #interaction term not sig.
vif(non_cucks_male.n_off.3)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.6)  #models are not significantly different



##-------------------cuck tests--------------------------

#select cuck offspring
cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for cucks_female  birds------------
cucks_female.n_off.1 <- glmer(n_off ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                  data = cucks_female_birds,
                                  family=binomial,
                                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(cucks_female.n_off.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_female.n_off.1, plot = T)

#check for collinearirty with vif
vif(cucks_female.n_off.1)

##plot basic model with no interactions
#plot recruitment to n_off for different Fs

cucks_female.n_off.indF.plot <- ggplot(cucks_female_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.n_off.mumF.plot <- ggplot(cucks_female_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.n_off.dadF.plot <- ggplot(cucks_female_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to n_off", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.indF.plot <- ggplot(male_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.mumF.plot <- ggplot(male_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.dadF.plot <- ggplot(male_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to n_off", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


n_off.global <- (cucks_female.n_off.indF.plot + male.n_off.indF.plot) / (cucks_female.n_off.mumF.plot + male.n_off.mumF.plot) / (cucks_female.n_off.dadF.plot + male.n_off.dadF.plot)

ggsave("Figures/n_off_global.png", n_off.global)

#----------------cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
cucks_female.n_off.2 <- update(cucks_female.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                 optCtrl = list(maxfun = 100000000)))

summary(cucks_female.n_off.2) #interaction term not significant.
vif(cucks_female.n_off.2)
anova(cucks_female.n_off.1,cucks_female.n_off.2) #models not significantly different


#check dadF help interaction. 
cucks_female.n_off.3 <- update(cucks_female.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.3) #interaction term very almost significant p = 0.0503 .
vif(cucks_female.n_off.3)
anova(cucks_female.n_off.1,cucks_female.n_off.3)  #models are not significantly different


#check FROH help interaction. 
cucks_female.n_off.4 <- update(cucks_female.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.4) #interaction term not sig.
vif(cucks_female.n_off.3)
anova(cucks_female.n_off.1,cucks_female.n_off.3)  #models are significantly different

#----------------------------------cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
cucks_female.n_off.4 <- update(cucks_female.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.4) #interaction term not significant.
vif(cucks_female.n_off.2)
anova(cucks_female.n_off.1,cucks_female.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_female.n_off.5 <- update(cucks_female.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.5) #interaction term not sig .
vif(cucks_female.n_off.3)
anova(cucks_female.n_off.1,cucks_female.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_female.n_off.6 <- update(cucks_female.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.6) #interaction term not sig.
vif(cucks_female.n_off.3)
anova(cucks_female.n_off.1,cucks_female.n_off.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select bastard
cucks_male_birds <-filter(na.omit(male_birds, EPP == 0) )

#basic model
cucks_male.n_off.1 <- glmer(n_off ~ mumF + dadF +FROH + BrMF+  help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                                data = cucks_male_birds,
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_male.n_off.1, plot = T)

#check for collinearirty with vif
vif(cucks_male.n_off.1)


#-----------------cucks_male helper interactions-----------





#check mumF help interaction. 
cucks_male.n_off.2 <- update(cucks_male.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.2) #interaction term not significant.
vif(cucks_male.n_off.2)
anova(cucks_male.n_off.1,cucks_male.n_off.2) #models not significantly different


#check dadF help interaction. 
cucks_male.n_off.3 <- update(cucks_male.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.3) #interaction term not significant
vif(cucks_male.n_off.3)
anova(cucks_male.n_off.1,cucks_male.n_off.3)  #models are not significantly different


#check FROH help interaction. 
cucks_male.n_off.4 <- update(cucks_male.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.4) #interaction term not sig. Nearly though
vif(cucks_male.n_off.3)
anova(cucks_male.n_off.1,cucks_male.n_off.3)  #models are not significantly different


#-------------------cucks_male mumage interaction-------

#check mumF MumAge interaction. 
cucks_male.n_off.4 <- update(cucks_male.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.4) #interaction term not significant.
vif(cucks_male.n_off.2)
anova(cucks_male.n_off.1,cucks_male.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_male.n_off.5 <- update(cucks_male.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.5) #interaction term not sig .
vif(cucks_male.n_off.3)
anova(cucks_male.n_off.1,cucks_male.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_male.n_off.6 <- update(cucks_male.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.6) #interaction term not sig.
vif(cucks_male.n_off.3)
anova(cucks_male.n_off.1,cucks_male.n_off.6)  #models are not significantly different



##-----------------coverage subsample--------------


# select males and females, birds that are already dead and never translocated, and coverage above 2X

cov_female_birds <- stat.birds %>% filter(Sex == 0, event == 0, Coverage > 2, CovMum > 2, CovDad > 2)

cov_male_birds <- stat.birds %>% filter(Sex == 1, event == 0, Coverage > 2, CovMum > 2, CovDad > 2)

#--------------------------juvenile survival---------------

#test if maternal or paternal roh affects the juvenile survival
#do this using 'n_off' variable, which is 0 for bird died before 1 year old
#this is a success or failure, so should use binomial distribution

#-------------------------tests for cov_female birds------------
cov_female.n_off.1 <- glmer(n_off ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                data = cov_female_birds,
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

summary(cov_female.n_off.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cov_female.n_off.1, plot = T)

#check for collinearirty with vif
vif(cov_female.n_off.1)

##plot basic model with no interactions
#plot recruitment to n_off for different Fs

female.n_off.indF.plot <- ggplot(female_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.n_off.mumF.plot <- ggplot(female_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


female.n_off.dadF.plot <- ggplot(female_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to n_off", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.indF.plot <- ggplot(male_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.mumF.plot <- ggplot(male_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.dadF.plot <- ggplot(male_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to n_off", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


n_off.global <- (female.n_off.indF.plot + male.n_off.indF.plot) / (female.n_off.mumF.plot + male.n_off.mumF.plot) / (female.n_off.dadF.plot + male.n_off.dadF.plot)

ggsave("Figures/n_off_global.png", n_off.global)

#----------------female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
female.n_off.2 <- update(female.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.2) #interaction term not significant.
vif(female.n_off.2)
anova(female.n_off.1,female.n_off.2) #models not significantly different


#check dadF help interaction. 
female.n_off.3 <- update(female.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.3) #interaction term very almost significant p = 0.0503 .
vif(female.n_off.3)
anova(female.n_off.1,female.n_off.3)  #models are not significantly different

#visually asses the difference.
#It seems that having helpers means that birds are LESS likely to recruit to n_off when dads are inbred
female.n_off.dadF.helpYES.plot <- ggplot(female_birds[female_birds$help == 1,], aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  scale_y_continuous(limits = c(0, 0.25))+
  labs( x = "Recruitment to n_off", y = "Paternal F", title = "Helpers in Nest")+
  theme_bw()+
  theme(text = element_text(size = 20))

female.n_off.dadF.helpNO.plot <- ggplot(female_birds[female_birds$help == 0,], aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  scale_y_continuous(limits = c(0, 0.25))+
  labs( x = "Recruitment to n_off", y = "Paternal F", title = "No Helpers")+
  theme_bw()+
  theme(text = element_text(size = 20))

n_off.female.help <- female.n_off.dadF.helpYES.plot | female.n_off.dadF.helpNO.plot 

ggsave("Figures/n_off_female_help.png", n_off.female.help)

#check FROH help interaction. 
female.n_off.4 <- update(female.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.4) #interaction term not sig.
vif(female.n_off.3)
anova(female.n_off.1,female.n_off.3)  #models are significantly different

#----------------------------------female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
female.n_off.4 <- update(female.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.4) #interaction term not significant.
vif(female.n_off.2)
anova(female.n_off.1,female.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
female.n_off.5 <- update(female.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.5) #interaction term not sig .
vif(female.n_off.3)
anova(female.n_off.1,female.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
female.n_off.6 <- update(female.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(female.n_off.6) #interaction term not sig.
vif(female.n_off.3)
anova(female.n_off.1,female.n_off.6)  #models are not significantly different

#-------------------------tests for male birds------------
#basic model
male.n_off.1 <- glmer(n_off ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                          data = male_birds,
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.1) #here FROH and mumage arte sig

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = male.n_off.1, plot = T)

#check for collinearirty with vif
vif(male.n_off.1)

#-----------------male helper interactions-----------

#check mumF help interaction. 
male.n_off.2 <- update(male.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.2) #interaction term not significant.
vif(male.n_off.2)
anova(male.n_off.1,male.n_off.2) #models not significantly different


#check dadF help interaction. 
male.n_off.3 <- update(male.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.3) #interaction term not significant
vif(male.n_off.3)
anova(male.n_off.1,male.n_off.3)  #models are not significantly different


#check FROH help interaction. 
male.n_off.4 <- update(male.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.4) #interaction term not sig. Nearly though
vif(male.n_off.3)
anova(male.n_off.1,male.n_off.3)  #models are not significantly different


#-------------------male mumage interaction-------

#check mumF MumAge interaction. 
male.n_off.4 <- update(male.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.4) #interaction term not significant.
vif(male.n_off.2)
anova(male.n_off.1,male.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
male.n_off.5 <- update(male.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.5) #interaction term not sig .
vif(male.n_off.3)
anova(male.n_off.1,male.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
male.n_off.6 <- update(male.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(male.n_off.6) #interaction term not sig.
vif(male.n_off.3)
anova(male.n_off.1,male.n_off.6)  #models are not significantly different


##--------------------------------EPP tests-------------------------------------


#select non cuck offspring
non_cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for non_cucks_female  birds------------
non_cucks_female.n_off.1 <- glmer(n_off ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                      data = non_cucks_female_birds,
                                      family=binomial,
                                      control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(non_cucks_female.n_off.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_female.n_off.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_female.n_off.1)

##plot basic model with no interactions
#plot recruitment to n_off for different Fs

non_cucks_female.n_off.indF.plot <- ggplot(non_cucks_female_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.n_off.mumF.plot <- ggplot(non_cucks_female_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


non_cucks_female.n_off.dadF.plot <- ggplot(non_cucks_female_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to n_off", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.indF.plot <- ggplot(male_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.mumF.plot <- ggplot(male_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.dadF.plot <- ggplot(male_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to n_off", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


n_off.global <- (non_cucks_female.n_off.indF.plot + male.n_off.indF.plot) / (non_cucks_female.n_off.mumF.plot + male.n_off.mumF.plot) / (non_cucks_female.n_off.dadF.plot + male.n_off.dadF.plot)

ggsave("Figures/n_off_global.png", n_off.global)

#----------------non_cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
non_cucks_female.n_off.2 <- update(non_cucks_female.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                         optCtrl = list(maxfun = 100000000)))

summary(non_cucks_female.n_off.2) #interaction term not significant.
vif(non_cucks_female.n_off.2)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.2) #models not significantly different


#check dadF help interaction. 
non_cucks_female.n_off.3 <- update(non_cucks_female.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.3) #interaction term very almost significant p = 0.0503 .
vif(non_cucks_female.n_off.3)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_female.n_off.4 <- update(non_cucks_female.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.4) #interaction term not sig.
vif(non_cucks_female.n_off.3)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.3)  #models are significantly different

#----------------------------------non_cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
non_cucks_female.n_off.4 <- update(non_cucks_female.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.4) #interaction term not significant.
vif(non_cucks_female.n_off.2)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_female.n_off.5 <- update(non_cucks_female.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.5) #interaction term not sig .
vif(non_cucks_female.n_off.3)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_female.n_off.6 <- update(non_cucks_female.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(non_cucks_female.n_off.6) #interaction term not sig.
vif(non_cucks_female.n_off.3)
anova(non_cucks_female.n_off.1,non_cucks_female.n_off.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select non bastard
non_cucks_male_birds <-filter(male_birds, EPP == 0) 

#basic model
non_cucks_male.n_off.1 <- glmer(n_off ~ mumF + dadF +FROH +   help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                                    data = non_cucks_male_birds,
                                    family=binomial,
                                    control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = non_cucks_male.n_off.1, plot = T)

#check for collinearirty with vif
vif(non_cucks_male.n_off.1)


#-----------------non_cucks_male helper interactions-----------





#check mumF help interaction. 
non_cucks_male.n_off.2 <- update(non_cucks_male.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.2) #interaction term not significant.
vif(non_cucks_male.n_off.2)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.2) #models not significantly different


#check dadF help interaction. 
non_cucks_male.n_off.3 <- update(non_cucks_male.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.3) #interaction term not significant
vif(non_cucks_male.n_off.3)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.3)  #models are not significantly different


#check FROH help interaction. 
non_cucks_male.n_off.4 <- update(non_cucks_male.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.4) #interaction term not sig. Nearly though
vif(non_cucks_male.n_off.3)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.3)  #models are not significantly different



#-------------------non_cucks_male mumage interaction-------


#check mumF MumAge interaction. 
non_cucks_male.n_off.4 <- update(non_cucks_male.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.4) #interaction term not significant.
vif(non_cucks_male.n_off.2)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
non_cucks_male.n_off.5 <- update(non_cucks_male.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.5) #interaction term not sig .
vif(non_cucks_male.n_off.3)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
non_cucks_male.n_off.6 <- update(non_cucks_male.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(non_cucks_male.n_off.6) #interaction term not sig.
vif(non_cucks_male.n_off.3)
anova(non_cucks_male.n_off.1,non_cucks_male.n_off.6)  #models are not significantly different



##-------------------cuck tests--------------------------

#select cuck offspring
cucks_female_birds <-filter(female_birds, EPP == 0) 



#-------------------------tests for cucks_female  birds------------
cucks_female.n_off.1 <- glmer(n_off ~ mumF + dadF + FROH + help + MumAge + (1 | birth_year) + (1 | mum) + (1 | dad),
                                  data = cucks_female_birds,
                                  family=binomial,
                                  control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=2e5)))

summary(cucks_female.n_off.1)


#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_female.n_off.1, plot = T)

#check for collinearirty with vif
vif(cucks_female.n_off.1)

##plot basic model with no interactions
#plot recruitment to n_off for different Fs

cucks_female.n_off.indF.plot <- ggplot(cucks_female_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Individual F", title = "Females")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.n_off.mumF.plot <- ggplot(cucks_female_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "", y = "Maternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))


cucks_female.n_off.dadF.plot <- ggplot(cucks_female_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "pink", fill = "pink")+
  geom_boxplot(width = 0.1, fill = "pink")+
  labs( x = "Recruitment to n_off", y = "Paternal F")+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.indF.plot <- ggplot(male_birds, aes(n_off, FROH))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", title = "Males")+
  geom_signif(comparisons = list(c("0","1")))+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.mumF.plot <- ggplot(male_birds, aes(n_off,mumF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "", y = "", )+
  theme_bw()+
  theme(text = element_text(size = 20))

male.n_off.dadF.plot <- ggplot(male_birds, aes(n_off,dadF))+
  geom_violin()+
  geom_jitter(color = "blue", fill = "blue")+
  geom_boxplot(width = 0.1, fill = "blue")+
  labs( x = "Recruitment to n_off", y= "")+
  theme_bw()+
  theme(text = element_text(size = 20))


n_off.global <- (cucks_female.n_off.indF.plot + male.n_off.indF.plot) / (cucks_female.n_off.mumF.plot + male.n_off.mumF.plot) / (cucks_female.n_off.dadF.plot + male.n_off.dadF.plot)

ggsave("Figures/n_off_global.png", n_off.global)

#----------------cucks_female helper interactions-------------


#check for interaction effects.
#here i am interested if the presence of helpers or maternal age will affect the inbreeding depression.

#check mumF help interaction. 
cucks_female.n_off.2 <- update(cucks_female.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="nloptwrap",
                                                                                                 optCtrl = list(maxfun = 100000000)))

summary(cucks_female.n_off.2) #interaction term not significant.
vif(cucks_female.n_off.2)
anova(cucks_female.n_off.1,cucks_female.n_off.2) #models not significantly different


#check dadF help interaction. 
cucks_female.n_off.3 <- update(cucks_female.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.3) #interaction term very almost significant p = 0.0503 .
vif(cucks_female.n_off.3)
anova(cucks_female.n_off.1,cucks_female.n_off.3)  #models are not significantly different


#check FROH help interaction. 
cucks_female.n_off.4 <- update(cucks_female.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.4) #interaction term not sig.
vif(cucks_female.n_off.3)
anova(cucks_female.n_off.1,cucks_female.n_off.3)  #models are significantly different

#----------------------------------cucks_female mumage interactions--------------
##check for mother age interactions

#check mumF MumAge interaction. 
cucks_female.n_off.4 <- update(cucks_female.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.4) #interaction term not significant.
vif(cucks_female.n_off.2)
anova(cucks_female.n_off.1,cucks_female.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_female.n_off.5 <- update(cucks_female.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.5) #interaction term not sig .
vif(cucks_female.n_off.3)
anova(cucks_female.n_off.1,cucks_female.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_female.n_off.6 <- update(cucks_female.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="nloptwrap"))

summary(cucks_female.n_off.6) #interaction term not sig.
vif(cucks_female.n_off.3)
anova(cucks_female.n_off.1,cucks_female.n_off.6)  #models are not significantly different






#-------------------------tests for male birds------------

#select bastard
cucks_male_birds <-filter(na.omit(male_birds, EPP == 0) )

#basic model
cucks_male.n_off.1 <- glmer(n_off ~ mumF + dadF +FROH + BrMF+  help+  MumAge +(1 | birth_year) + (1 | mum) + (1 | dad),
                                data = cucks_male_birds,
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.1) #here dadF is significant ###################

#use DHARMA to check model vs simulated residuals. Seems ok
simulateResiduals(fittedModel = cucks_male.n_off.1, plot = T)

#check for collinearirty with vif
vif(cucks_male.n_off.1)


#-----------------cucks_male helper interactions-----------





#check mumF help interaction. 
cucks_male.n_off.2 <- update(cucks_male.n_off.1, ~. +mumF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.2) #interaction term not significant.
vif(cucks_male.n_off.2)
anova(cucks_male.n_off.1,cucks_male.n_off.2) #models not significantly different


#check dadF help interaction. 
cucks_male.n_off.3 <- update(cucks_male.n_off.1, ~. +dadF*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.3) #interaction term not significant
vif(cucks_male.n_off.3)
anova(cucks_male.n_off.1,cucks_male.n_off.3)  #models are not significantly different


#check FROH help interaction. 
cucks_male.n_off.4 <- update(cucks_male.n_off.1, ~. +FROH*help ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.4) #interaction term not sig. Nearly though
vif(cucks_male.n_off.3)
anova(cucks_male.n_off.1,cucks_male.n_off.3)  #models are not significantly different


#-------------------cucks_male mumage interaction-------

#check mumF MumAge interaction. 
cucks_male.n_off.4 <- update(cucks_male.n_off.1, ~. +mumF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.4) #interaction term not significant.
vif(cucks_male.n_off.2)
anova(cucks_male.n_off.1,cucks_male.n_off.2) #models not significantly different


#check dadF MumAge interaction. 
cucks_male.n_off.5 <- update(cucks_male.n_off.1, ~. +dadF*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.5) #interaction term not sig .
vif(cucks_male.n_off.3)
anova(cucks_male.n_off.1,cucks_male.n_off.5)  #models are not significantly different


#check FROH MumAge interaction. 
cucks_male.n_off.6 <- update(cucks_male.n_off.1, ~. +FROH*MumAge ,control=glmerControl(optimizer="bobyqa"))

summary(cucks_male.n_off.6) #interaction term not sig.
vif(cucks_male.n_off.3)
anova(cucks_male.n_off.1,cucks_male.n_off.6)  #models are not significantly different





