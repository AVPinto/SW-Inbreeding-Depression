#Inbreeding depression in the Seychelles warbler
#Alessandro V Pinto

install.packages("gtsummary")
install.packages("glue")


#packages
library(viridis)
library(stargazer)
library(see)
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
library(modelsummary)
library(glmer)
library(arm)
library(ggeffects)
library(ggggeffects)
library(gtsummary)
library(glue)
library(interactions)






#----------------------------------import data--------------

##Factorise
stat.birds <- read.csv("clean_data.csv")
stat.birds$helpF <- as.factor(stat.birds$help)
stat.birds$BreedSeason <- as.factor(stat.birds$BreedSeason)
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$mum <- as.factor(stat.birds$mum)
stat.birds$dad <- as.factor(stat.birds$dad)
stat.birds$birth_yearF <- as.factor(stat.birds$birth_year)
stat.birds$adulthood <- as.factor(stat.birds$adulthood)


dim(stat.birds)  #n = 1385

#check FROH against coverage

ggplot(stat.birds, aes(x = Coverage, y = LargeFROH))+
  geom_point()+
  geom_smooth()

#some decline in the low end.

ggplot((stat.birds %>% filter(Coverage > 0.25)), aes(x = Coverage, y = LargeFROH))+
  geom_point()+
  geom_smooth()

#nice flat line after filtering for coverage

#filter for birds that have died, and have coverage > 0.25

just_birds <- stat.birds %>% filter(event == 0, Coverage >= 0.25 )
dim(just_birds) #n = 1188

#---------------------lifespan models------------------

hist(just_birds$lifespan, breaks = 16)

#quite zero inflated

#check model distribution families and transformations for best fit

just.lifespan.poisson <- glmmTMB(lifespan ~ LargeFROH + Sex + help + mean_total_rain + mean_rain_cv +
                                    birth_year + (1 | mum) + (1 | dad) ,
                           data = just_birds,
                           ziformula=~1,
                           family=poisson
)


simulateResiduals(just.lifespan.poisson, plot = TRUE)
#deviation in the qq, residuals whack

#check for collinearity and overdispersion
check_collinearity(just.lifespan.poisson) #vif around one
check_overdispersion(just.lifespan.poisson) #no overdisp


#nbinom1 
just.lifespan.nbinom1 <- update(just.lifespan.poisson, family=nbinom1())

simulateResiduals(just.lifespan.nbinom1, plot = TRUE)
#qq is better, residuals still red but the points are not too bad

#check for collinearity and overdispersion
check_collinearity(just.lifespan.nbinom1) #vif around one
check_overdispersion(just.lifespan.nbinom1) #no overdisp


#nbinom2 
just.lifespan.nbinom2 <-  update(just.lifespan.poisson, family=nbinom2())

simulateResiduals(just.lifespan.nbinom2, plot = TRUE)
#qq same,  residuals looking bad



#AIC comparisons
AIC(just.lifespan.nbinom1, just.lifespan.nbinom2, just.lifespan.poisson)
#######################df      AIC
#just.lifespan.nbinom1 11 5776.268
#just.lifespan.nbinom2 11 5760.495
#just.lifespan.poisson 10 6090.430

summary(just.lifespan.nbinom2)
#nbinom2 has best AIC. However the data has linear observed vs residuals rather than quadratic
#so nbinom1 is more suitable

summary(just.lifespan.nbinom1)

#plot
ggplot(just_birds, aes(LargeFROH, lifespan)) +
  geom_point()+
  geom_smooth(method = "lm")


##now check for interactions

#FROH * sex interaction

just.lifespan.nbinom1.c <- update(just.lifespan.nbinom1, ~ . +  LargeFROH*Sex)

summary(just.lifespan.nbinom1.c) #model converges. int not sig

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.c) #anova says 0.8, not useful to include

#FROH * mean_rain interaction

just.lifespan.nbinom1.d <- update(just.lifespan.nbinom1, ~ . +  LargeFROH*mean_total_rain)

summary(just.lifespan.nbinom1.d) #model converges. int not sig

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.d) #anova says 0.55, not useful to include


#FROH * rain variance

just.lifespan.nbinom1.e <- update(just.lifespan.nbinom1, ~ . +  LargeFROH*mean_rain_cv)

summary(just.lifespan.nbinom1.e) #model converges. int not sig

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.e) #anova says 0.74, not useful to include


#FROH * Help interaction

just.lifespan.nbinom1.f <- update(just.lifespan.nbinom1, ~ . +  LargeFROH*help)

summary(just.lifespan.nbinom1.f) #model converges. int not sig

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.f) #anova says 0.51


##test for quadratic terms... perhaps the ID or weather is non-linear

#quadratic froh

just.lifespan.nbinom1.g <- update(just.lifespan.nbinom1, ~ . +  I(LargeFROH^2))

summary(just.lifespan.nbinom1.g) #model converges. quadratic effect not sig

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.g) #anova says  not useful to include

#quadratic rain

just.lifespan.nbinom1.h <- update(just.lifespan.nbinom1, ~ . +  I(mean_total_rain^2))

summary(just.lifespan.nbinom1.h) #model does not converge.


#drop lifespan

just.lifespan.nbinom1.x <- update(just.lifespan.nbinom1, ~ . - lifespan)

summary(just.lifespan.nbinom1.x) #model does not converge.




#so in summary model first model is best1
summary(just.lifespan.nbinom1)

#output table
tbl_regression(just.lifespan.nbinom1,
               intercept = T,
               show_single_row = "Sex",
               pvalue_fun = label_style_pvalue(digits = 3),
               estimate_fun = label_style_number(digits = 4),
               label = list(LargeFROH = "FROH > 3.3Mb", 
                            help = "Helper in natal territory",
                            mean_total_rain = "Mean annual rainfall over lifespan",
                            mean_rain_cv = "Mean variance in annual rainfall over lifespan",
                            birth_year = "Hatch year")
               ) %>%
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) 

####plots

##plot with ggpredict
#make ggpredict object
ggpred.just.lifespan.nbinom1 <-  ggpredict(just.lifespan.nbinom1,
                                           terms = "LargeFROH[all]")

plot(ggpred.just.lifespan.nbinom1)

#plot with ggggeffects for ggplot utility
autoplot(ggpred.just.lifespan.nbinom1)+
  geom_expected_line()+
  geom_CI_ribbon()+
  theme_classic2()+
  labs( x = "fROH > 3.3Mb", y= "Predicted Lifespan")+
  theme(text = element_text(size = 20)) +
  ggtitle("Fig. 2 - Predicted fit of GLMM Lifespan ~ FROH")+
  layer_fit_data(alpha = 0.2)
  


#ggplot
just.lifespan.largefroh.plot <- ggplot(just_birds, aes(LargeFROH, lifespan)) +
  geom_point()+
  geom_smooth(method = "lm" )+
  theme_classic2()+
  labs( x = "fROH > 3.3Mb", y= "Lifespan")+
  theme(text = element_text(size = 20))

  














#------------------------first year survival----------------------

#use binomial family as survival is binary

just.fys.1 <- glmer(adulthood ~ LargeFROH  + Sex + help + birth_total_rain + BirthRainCV + 
                      (1 | birth_year) + (1 | mum) + (1 | dad) ,
                    data = just_birds,
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))

)

summary(just.fys.1)

#model does not converge, try rescaling as suggested

just.fys.rescale <- glmer(adulthood ~ rescale(LargeFROH) + Sex + help  + rescale(birth_total_rain) + rescale(BirthRainCV)+
                      (1 | birth_year) + (1 | mum) + (1 | dad) ,
                    data = just_birds,
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))
                    
)

summary(just.fys.rescale)

simulateResiduals(just.fys.rescale, plot = T) #looks fine

#check for collinearity and overdispersion
check_collinearity(just.fys.rescale) #vif around one
check_overdispersion(just.fys.rescale) #no overdisp

###try interactions

#froh * rain

just.fys.rescale.b <- update(just.fys.rescale, ~ . + rescale(LargeFROH)*rescale(birth_total_rain))

summary(just.fys.rescale.b)# model converges. int not sig

anova(just.fys.rescale, just.fys.rescale.b) #anova says 0.4497

#froh * help

just.fys.rescale.c   <- update(just.fys.rescale, ~ . + rescale(LargeFROH)*help)

summary(just.fys.rescale.c)# model converges. int not sig

anova(just.fys.rescale, just.fys.rescale.c) #anova says 0.7322

#froh * sex interaction

just.fys.rescale.d <- update(just.fys.rescale, ~ . + rescale(LargeFROH)*Sex)

summary(just.fys.rescale.d)# model converges. int not sig

anova(just.fys.rescale, just.fys.rescale.d) #anova says 0.7115


#check quadratics 

#quadratic froh

just.fys.rescale.e <- update(just.fys.rescale, ~ . + I(rescale(LargeFROH^2)))

summary(just.fys.rescale.e)# model converges. int not sig

anova(just.fys.rescale, just.fys.rescale.e) #anova says 0.5

#quadratic rain

just.fys.rescale.f <- update(just.fys.rescale, ~ . + I(rescale(birth_total_rain^2)))

summary(just.fys.rescale.f)# model converges. quad IS sig

anova(just.fys.rescale, just.fys.rescale.f) #anova says 0.0007331 keep it!

#interaction with quadratic?

just.fys.rescale.g <- update(just.fys.rescale, ~ . + I(rescale(birth_total_rain^2))*LargeFROH  )

summary(just.fys.rescale.g)# model converges. interaction not significant

anova(just.fys.rescale, just.fys.rescale.g) 

####summary

#In summary we keep model f
summary(just.fys.rescale.f)


#output table
tbl_regression(just.fys.rescale.f, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            help = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall in hatch year",
                            "rescale(BirthRainCV)" = "Variance in rainfall in hatch year",
                            "I(rescale(birth_total_rain^2))" = "Quadratic of Annual rainfall"))%>% 

  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) 


#plot with ggpredict. doesn't quite work

##plots
#with ggpredict
##plot with ggpredict
#make ggpredict object
ggpred.just.fys.rescale.f <-  ggpredict(just.fys.rescale.f,
                                           terms = "LargeFROH[all]")

plot(ggpred.just.fys.rescale.f)

#plot with ggggeffects for ggplot utility
autoplot(ggpred.just.fys.rescale.f)+
  geom_expected_line()+
  geom_CI_ribbon()+
  theme_classic2()+
  labs( x = "fROH > 3.3Mb", y= "Predicted Lifespan")+
  theme(text = element_text(size = 20)) +
  ggtitle("Fig. 1 - Predicted fit of GLMM Lifespan ~ FROH")+
  layer_fit_data(alpha = 0.2)


predict.just.fys.rescale.f <- ggpredict(just.fys.rescale.f, terms="LargeFROH [all]")

plot(predict.just.fys.rescale.f) 

#plot with ggplot
  
  ggplot(just_birds, aes(LargeFROH, adulthood)) +
    geom_boxplot()+
    geom_jitter(alpha = 0.5)+
  theme_classic2()+
  labs( x = "fROH > 3.3Mb", y= "Recruitment to Adulthood")+
  theme(text = element_text(size = 20))+
  ggtitle("Fig. 2")+  
    geom_signif(comparisons = list(c("0","1")),
                                  map_signif_level = TRUE)
  
  
  ggplot(just_birds, aes(birth_total_rain, adulthood)) +
    geom_jitter()+
    theme_classic2()+
    labs( x = "Total annual rainfall", y= "Recruitment to Adulthood")+
    theme(text = element_text(size = 20))+
    ggtitle("Fig. 2")
  


  
  
  
  
  
  #----------------------------------lifetime reproductive success--------------  
  
  
  hist(just_birds$n_off, breaks = 16)
  
  #quite zero inflated
  
just.lrs.poisson <- glmmTMB(n_off ~ LargeFROH + help + Sex + 
                              mean_total_rain + mean_rain_cv + birth_year +
                      (1 | mum) + (1 | dad),
                      data = just_birds,
                      ziformula=~1,
                      family=poisson)


summary(just.lrs.poisson)

simulateResiduals(just.lrs.poisson, plot = T) #a bit wiggly there

#check for collinearity and overdispersion
check_collinearity(just.lrs.poisson) #vif around one
check_overdispersion(just.lrs.poisson) #no overdisp



#nbinom1 
just.lrs.nbinom1 <- update(just.lrs.poisson, family=nbinom1())

simulateResiduals(just.lrs.nbinom1, plot = TRUE)
#qq is better, but residuals are slightly bad

#check for collinearity and overdispersion
check_collinearity(just.lrs.nbinom1) #vif around one
check_overdispersion(just.lrs.nbinom1) # overdisp DETECTED
#some error messages

#nbinom2
just.lrs.nbinom2 <- update(just.lrs.poisson, family=nbinom2())
#some error messages. using different optimizer


simulateResiduals(just.lrs.nbinom2, plot = TRUE)
#qq is slightly worse, but residuals are slightly bad

#check for collinearity and overdispersion
check_collinearity(just.lrs.nbinom2) #vif around one
check_overdispersion(just.lrs.nbinom2) #no overdisp

#AIC comparisons
AIC(just.lrs.nbinom1, just.lrs.nbinom2, just.lrs.poisson)

#df      AIC
#ust.lrs.nbinom1 11       NA doesn't convert
#just.lrs.nbinom2 11 3986.743
#just.lrs.poisson 10 4162.900

summary(just.lrs.nbinom1)
summary(just.lrs.nbinom2)
summary(just.lrs.poisson)

#nbinoms have better AIC but are unstable. lets use poisson.




#try interactions

#froh*help

just.lrs.poisson.b <- update(just.lrs.poisson , ~ . + LargeFROH*help)

summary(just.lrs.poisson.b) #model converges, not sig

anova(just.lrs.poisson, just.lrs.poisson.b) #anova says not quite

#froh*sex

just.lrs.poisson.c <- update(just.lrs.poisson , ~ . + LargeFROH*Sex)

summary(just.lrs.poisson.c) #model converges, anova says not quite 

anova(just.lrs.poisson, just.lrs.poisson.c) 

#froh*mean_rain

just.lrs.poisson.d <- update(just.lrs.poisson , ~ . + LargeFROH*mean_total_rain)

summary(just.lrs.poisson.d) #not significant

anova(just.lrs.poisson, just.lrs.poisson.d) #anova says no!?


#froh*rain variance

just.lrs.poisson.e <- update(just.lrs.poisson , ~ . + LargeFROH*mean_rain_cv)

summary(just.lrs.poisson.e) #model converges, is not significant!

anova(just.lrs.poisson, just.lrs.poisson.e) #anova says no

#quadratics

just.lrs.poisson.f <- update(just.lrs.poisson , ~ . + I(mean_total_rain^2))

summary(just.lrs.poisson.f) #model does not converge


####summary
#in sumarry no interactions are significant so we go with base model

#we will output the base model first just in case
summary(just.lrs.poisson)

tbl_regression(just.lrs.poisson, intercept = T,
               show_single_row = "Sex",
               label = list("LargeFROH" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "mean_total_rain"	 = "Mean annual rainfall during lifespan",
                            "mean_rain_cv" = "Mean variance annual rainfall during lifespan",
                            "birth_year" = "Birth year",
                            "lifespan" ="Lifespan")
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  


##plots
#with ggpredict

#make ggpredict object
ggpred.just.lrs.poisson <-  ggpredict(just.lrs.poisson,
                                        terms = "LargeFROH[all]")

plot(ggpred.just.lrs.poisson)

#plot with ggggeffects for ggplot utility
autoplot(ggpred.just.lrs.poisson)+
  geom_expected_line()+
  geom_CI_ribbon()+
  theme_classic2()+
  labs( x = "fROH > 3.3Mb", y= "Lifetime reproductive success")+
  theme(text = element_text(size = 20)) +
  ggtitle("Fig. 3 - Predicted fit of ZI-GLMM LRS ~ FROH")+
  layer_fit_data(alpha = 0.2)



#output model with interaction - not real any more
summary(just.lrs.poisson.d)


#output table
tbl_regression(just.lrs.poisson.d, intercept = T,
               show_single_row = "Sex",
               label = list("LargeFROH" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "mean_total_rain"	 = "Mean annual rainfall during lifespan",
                            "mean_rain_cv" = "Mean variance annual rainfall during lifespan",
                            "birth_year" = "Birth year",
                            "lifespan" ="Lifespan")
)%>%
               
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))                

#plot with ggpredict. doesn't quite work

##plots

#try sliceplotting

interact_plot(just.lrs.poisson.d,
              pred = LargeFROH,
              modx = mean_total_rain,
              modx.values = "terciles",
              interval = T,
              plot.points = T)+
  theme_classic2()+
  labs( x = "fROH > 3.3Mb",
        y= "Lifetime reproductive success",
        fill = "Mean annual Rainfall over lifespan")+
  theme(text = element_text(size = 20)) +
  ggtitle("Fig. 3.2 - Predicted fit of GLMM LRS ~ FROH
          with Mean Annual Rainfall interaction")


#------------------------cox model---------------

#right censor still alive and translocated birds
#use stat.birds rather than just_birds as we can include censored birds


stat.cox.me <- coxme(Surv(lifespan, event) ~  rescale(LargeFROH) +Sex + help + 
                       mean_total_rain + mean_rain_cv +
                    (1 | mum) + (1 | dad) ,
                  data = stat.birds , 
                  control = coxme.control(iter.max  = 200)
                  ) 

#model converges ok

#some tests

#proportional hazards tests
cox.test <- cox.zph(stat.cox.me)
print(cox.test) #GLOBAL and mean_rain_cv violates the assumptions of proportional hazards

plot(cox.test, var = "mean_rain_cv")# seems to increase over time.
#maybe this is related to weather getting wackier recently

summary(stat.cox.me) #it's not significant


#we will remove it
stat.cox.me.a <- coxme(Surv(lifespan, event) ~  rescale(LargeFROH) +Sex + help + 
                         mean_total_rain +  
                         (1 | mum) + (1 | dad) ,
                       data = stat.birds , 
                       control = coxme.control(iter.max  = 200)
) 

cox.zph(stat.cox.me.a) #GLOBAL now conforms to ph assumptions

#now test interactions

#froh * sex

stat.cox.me.b <- coxme(Surv(lifespan, event) ~  rescale(LargeFROH) +Sex + help + 
                            mean_total_rain + rescale(LargeFROH)*Sex +
                            (1 | mum) + (1 | dad) ,
                          data = stat.birds , 
                          control = coxme.control(iter.max  = 200)
                       )

summary(stat.cox.me.b) #model converges. not significant

cox.zph(stat.cox.me.b) #fine

anova(stat.cox.me.a, stat.cox.me.b) #no sig diff

#froh * help

stat.cox.me.c <- coxme(Surv(lifespan, event) ~  rescale(LargeFROH) +Sex + help + 
                         mean_total_rain + rescale(LargeFROH)*help +
                         (1 | mum) + (1 | dad) ,
                       data = stat.birds , 
                       control = coxme.control(iter.max  = 200)
)

summary(stat.cox.me.c) #model converges. not significant

cox.zph(stat.cox.me.c) #fine

anova(stat.cox.me.a, stat.cox.me.c) #no sig diff

#froh * MumAge

#stat.cox.me.d <- coxme(Surv(lifespan, event) ~  rescale(LargeFROH) +Sex + help + MumAge +
#                         mean_total_rain + rescale(LargeFROH)*MumAge +
#                         (1 | mum) + (1 | dad) ,
#                       data = stat.birds , 
#                       control = coxme.control(iter.max  = 200)
#)
#
#summary(stat.cox.me.d) #model converges. not significant
#
#cox.zph(stat.cox.me.d) #fine
#
#anova(stat.cox.me.a, stat.cox.me.d) #no sig diff

#froh * rain

stat.cox.me.e <- coxme(Surv(lifespan, event) ~  rescale(LargeFROH) +Sex + help + 
                         mean_total_rain + rescale(LargeFROH)*mean_total_rain +
                         (1 | mum) + (1 | dad) ,
                       data = stat.birds , 
                       control = coxme.control(iter.max  = 200)
)

summary(stat.cox.me.e) #model converges. not significant

cox.zph(stat.cox.me.e) #fine

anova(stat.cox.me.a, stat.cox.me.e) #no sig diff

#quadratic rain

stat.cox.me.f <- coxme(Surv(lifespan, event) ~  rescale(LargeFROH) +Sex + help + 
                         mean_total_rain + I(mean_total_rain^2) +
                         (1 | mum) + (1 | dad) ,
                       data = stat.birds , 
                       control = coxme.control(iter.max  = 200)
)

summary(stat.cox.me.f) #model converges. quadratic rain IS significant

cox.zph(stat.cox.me.f) #GLOBAL is not fine but all variables are?
plot(cox.zph(stat.cox.me.f))

#adding the quadratic term improves the fit but makes GLOBAL not fit assumptions

#check for interactions with the term

stat.cox.me.g <- coxme(Surv(lifespan, event) ~  rescale(LargeFROH) +Sex + help + MumAge +
                         mean_total_rain + I(mean_total_rain^2) + I(mean_total_rain^2)*rescale(LargeFROH) +
                         (1 | mum) + (1 | dad) ,
                       data = stat.birds , 
                       control = coxme.control(iter.max  = 200)
)

summary(stat.cox.me.g) #model converges. interaction term not significant.

#in summary, we will go with a. Better to remove the quadratic rain term and be consistent.
  #Especially as there is no significant interactions with it

summary(stat.cox.me.a)
#table output
tbl_regression(stat.cox.me.a, intercept = T, show_single_row = "Sex")

tbl_regression(stat.cox.me.a, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            "help" = "Helper in natal territory",
                            "mean_total_rain"	 = "Mean annual rainfall during lifespan")
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  




#plot

#first bin the rohs

#sort out sensible bin sizes

cox.plot.birds <- stat.birds
cox.plot.birds$FROH_class <- 0
cox.plot.birds <- cox.plot.birds %>% select(-DeathEventL)


hist(cox.plot.birds$LargeFROH)

quantile(cox.plot.birds$LargeFROH, probs = c(0.1,0.9))


for (n in seq(length(cox.plot.birds$LargeFROH))) {
  
  if(cox.plot.birds$LargeFROH[n] < 0.16 ) {
    cox.plot.birds$FROH_class [n] <- "< 10%"
    
  } 
  
  if (cox.plot.birds$LargeFROH[n] > 0.16 && cox.plot.birds$LargeFROH[n] <= 0.3){
    cox.plot.birds$FROH_class [n] <- " 10%-90%"
  }

  if (cox.plot.birds$LargeFROH[n] > 0.3 ){
    cox.plot.birds$FROH_class [n] <- ">90%"
  }
}



cox.surv.froh <- survfit(Surv(lifespan,event)~ FROH_class, data = cox.plot.birds)


cox.surv.froh.plot <- ggsurvplot(cox.surv.froh, data = cox.plot.birds,
                                  xlab="Age",
                                  legend.title = "FROH >33MB quantile",
                                 legend.labs = c("10%-90%","<10%", "> 90%"),
                                  pval = FALSE,
                                  conf.int = FALSE,
                                  font.y = 20,
                                  font.x = 20,
                                  font.legend  = c(20, "blue"))



cox.surv.froh.plot+ggtitle("Fig. 4.1 - Cox Proportional Hazards Model with FROH Quantiles")

#ggpredict plot

cox.plot.ph.ggpedict.data <- ggpredict(cox.plot.ph , terms = "LargeFROH[all]" )

cox.plot.ph <- coxph(Surv(lifespan, event) ~  rescale(LargeFROH), 
                     data = cox.plot.birds )

plot(cox.plot.ph)


plot(cox.plot.ph.ggpedict.data, show_data = TRUE)  


predict.cox.me.a <- ggpredict(stat.cox.me.a, terms = "LargeFROH[all]")




