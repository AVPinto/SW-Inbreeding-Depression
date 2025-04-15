#Inbreeding depression in the Seychelles warbler
#Alessandro V Pinto

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
library(gtsummary)
library(glue)
install.packages("gtsummary")
install.packages("glue")

#----------------------------------import data--------------

##Factorise
stat.birds <- read.csv("clean_data.csv")
stat.birds$helpF <- as.factor(stat.birds$help)
stat.birds$BreedSeason <- as.factor(stat.birds$BreedSeason)
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$DeathEventL <- as.logical(stat.birds$DeathEvent)
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

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.f) #anova says 0.51, almmost useful to include!



##test for quadratic terms... perhaps the ID or weather is non-linear

#quadratic froh

just.lifespan.nbinom1.g <- update(just.lifespan.nbinom1, ~ . +  I(LargeFROH^2))

summary(just.lifespan.nbinom1.g) #model converges. quadratic effect not sig

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.g) #anova says  not useful to include

#quadratic rain

just.lifespan.nbinom1.h <- update(just.lifespan.nbinom1, ~ . +  I(mean_total_rain^2))

summary(just.lifespan.nbinom1.h) #model does not converge.



#so in summary model first model is best1
summary(just.lifespan.nbinom1)

#output table
tbl_regression(just.lifespan.nbinom1, intercept = T, show_single_row = "Sex")

#plot
plot(ggpredict(just.lifespan.nbinom1), show_data = TRUE, jitter = T) 

#ggplot
just.lifespan.largefroh.plot <- ggplot(just_birds, aes(LargeFROH, lifespan)) +
  geom_point()+
  geom_smooth(method = "lm" )+
  theme_classic2()+
  labs( x = "fROH > 3.3Mb", y= "Lifespan")+
  theme(text = element_text(size = 20))

  














#------------------------first year survival----------------------

#use binomial family as survival is binary

just.fys.1 <- glmer(adulthood ~ LargeFROH + + Sex + help + birth_total_rain + BirthRainCV + 
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

anova(just.fys.rescale, just.fys.rescale.g) #anova says 0.001826  keep it!


#In summary we keep model f
summary(just.fys.rescale.f)



#output table
tbl_regression(just.fys.rescale.f, intercept = T, show_single_row = "Sex")

#plot
predict.just.fys.rescale.f <- ggpredict(just.fys.rescale.f)



plot(predict.just.fys.rescale.f$LargeFROH, show_data = TRUE) 



#plot
  #need to plot this better
#just.fys.largefroh.plot <- 
  
  ggplot(just_birds, aes(LargeFROH, adulthood)) +
    geom_boxplot()+
  theme_classic2()+
  labs( x = "fROH > 3.3Mb", y= "Recruitment to Adulthood")+
  theme(text = element_text(size = 20))
  


  
  
  
  
  
  
  
#----------------------lifetime reproductive success----------
  #unfinished#
  
  hist(just_birds$n_off, breaks = 16)
  
  #quite zero inflated
  
just.lrs.poisson <- glmmTMB(n_off ~ LargeFROH + lifespan + help + MumAge + Sex + 
                              mean_total_rain + mean_rain_cv + birth_year +
                      + (1 | mum) + (1 | dad),
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
check_overdispersion(just.lrs.nbinom1) #no overdisp


#nbinom2
just.lrs.nbinom2 <- update(just.lrs.poisson, family=nbinom2())
#some error messages. using different optimizer


simulateResiduals(just.lrs.nbinom2, plot = TRUE)
#qq is better, but residuals are slightly bad

#check for collinearity and overdispersion
check_collinearity(just.lrs.nbinom2) #vif around one
check_overdispersion(just.lrs.nbinom2) #no overdisp

#AIC comparisons
AIC(just.lrs.nbinom1, just.lrs.nbinom2, just.lrs.poisson)

#df      AIC
#just.lrs.nbinom1 13 3230.470
#just.lrs.nbinom2 13 3218.410
#just.lrs.poisson 12 3386.051

#lets go with nbinom1 again, it seems good

summary(just.lrs.nbinom1)

#try interactions

#froh*help

just.lrs.nbinom1.b <- update(just.lrs.nbinom1 , ~ . + LargeFROH*help)

summary(just.lrs.nbinom1.b) #model converges, not sig

anova(just.lrs.nbinom1, just.lrs.nbinom1.b) #anova says no

#froh*sex

just.lrs.nbinom1.c <- update(just.lrs.nbinom1 , ~ . + LargeFROH*Sex - (1|MumID) - (1|DadID))

summary(just.lrs.nbinom1.c) #model does not converge! 

anova(just.lrs.nbinom1, just.lrs.nbinom1.c) 

#froh*mean_rain

just.lrs.nbinom1.d <- update(just.lrs.nbinom1 , ~ . + LargeFROH*mean_total_rain)

summary(just.lrs.nbinom1.d) #model converges, not sig

anova(just.lrs.nbinom1, just.lrs.nbinom1.c) #anova says no

