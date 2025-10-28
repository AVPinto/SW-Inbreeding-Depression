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
library(arm)
library(ggeffects)
library(ggggeffects)
library(gtsummary)
library(glue)
library(interactions)
library(ggplot2)
library(marginaleffects)

#----------------------------------import data--------------

##Read in data
stat.birds <- read.csv("data/clean_data.csv")
str(stat.birds)

#Factorise some stuff
stat.birds$h_countF <- as.factor(stat.birds$h_count)
stat.birds$BreedSeason <- as.factor(stat.birds$BreedSeason)
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$mum <- as.factor(stat.birds$mum)
stat.birds$dad <- as.factor(stat.birds$dad)
stat.birds$birth_yearF <- as.factor(stat.birds$birth_year)
stat.birds$adulthood <- as.factor(stat.birds$adulthood)
stat.birds$sib_pres <- as.factor(stat.birds$sib_pres)


dim(stat.birds)  #n = 1210

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

just_birds <- stat.birds %>% filter(event == 0,
                                    Coverage >= 0.25,
                                         # remove some NAs
                                      !is.na(mean_total_rain),
                                      !is.na(mean_rain_cv),
                                      !is.na(mum),
                                      !is.na(dad)
                                     )
dim(just_birds) #n = 1019


#--------------exploratory analysis--------------------

#are males and females similar?

with(just_birds, t.test(lifespan[Sex == 0], lifespan[Sex == 1]))
wilcox.test(as.numeric(just_birds$lifespan) ~ just_birds$Sex)


with(just_birds, t.test(n_off[Sex == 0], n_off[Sex == 1]))
wilcox.test(as.numeric(just_birds$n_off) ~ just_birds$Sex)


#---------------------lifespan models------------------

hist(just_birds$lifespan, breaks = 29)


#check model distribution families and transformations for best fit


just.lifespan.poisson <- glmmTMB(
  lifespan ~ rescale(LargeFROH) + Sex + h_countF +
    rescale(mean_total_rain) + rescale(mean_rain_cv) +
    birth_year + sib_pres + rescale(natal_group_size) +
    (1 | mum) + (1 | dad) ,
  data = just_birds,
  ziformula =  ~ 0,
  family = poisson
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
check_zeroinflation(just.lifespan.nbinom1)
testZeroInflation(just.lifespan.nbinom1)


#nbinom2 
just.lifespan.nbinom2 <-  update(just.lifespan.poisson, family=nbinom2())

simulateResiduals(just.lifespan.nbinom2, plot = TRUE)
#qq same,  significant outliers
check_overdispersion(just.lifespan.nbinom2)
check_zeroinflation(just.lifespan.nbinom2)


#AIC comparisons
AIC(just.lifespan.nbinom1, just.lifespan.nbinom2, just.lifespan.poisson)
#######################df      AIC
#just.lifespan.nbinom1 13 5375.993
#just.lifespan.nbinom2 13 5326.276
#just.lifespan.poisson 12 5668.534

#nbinoms are similar but nbinom2 residuals are a bad, so go for nbinom1

summary(just.lifespan.nbinom1)

#quick plot
ggplot(just_birds, aes(LargeFROH, lifespan)) +
  geom_point()+
  geom_smooth(method = "lm")


##now check for interactions

#FROH * sex interaction

just.lifespan.nbinom1.c <- update(just.lifespan.nbinom1, ~ . +  rescale(LargeFROH)*Sex)

summary(just.lifespan.nbinom1.c) #model converges. int not sig

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.c) #anova says 0.8, not useful to include

#FROH * mean_rain interaction

just.lifespan.nbinom1.d <- update(just.lifespan.nbinom1, ~ . +  rescale(LargeFROH)*rescale(mean_total_rain))

summary(just.lifespan.nbinom1.d) #model converges. int not sig

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.d) #anova says 0.55, not useful to include


#FROH * rain variance

just.lifespan.nbinom1.e <- update(just.lifespan.nbinom1, ~ . +  rescale(LargeFROH)*rescale(mean_rain_cv))

summary(just.lifespan.nbinom1.e) #model converges. int not sig

anova(just.lifespan.nbinom1, just.lifespan.nbinom1.e) #anova says 0.74, not useful to include


#FROH * Help interaction

just.lifespan.nbinom1.f <- update(just.lifespan.nbinom1, ~ . +  rescale(LargeFROH)*h_countF)

summary(just.lifespan.nbinom1.f) #model converges. F1 almost significant

print(anova(just.lifespan.nbinom1, just.lifespan.nbinom1.f)) #anova says no


#FROH * sib_pres interaction

just.lifespan.nbinom1.g <- update(just.lifespan.nbinom1.f, ~ . +  rescale(LargeFROH)*sib_pres)

summary(just.lifespan.nbinom1.g) #model converges. not sig

anova(just.lifespan.nbinom1.f, just.lifespan.nbinom1.g) #anova says no


#FROH * natal group size interaction

just.lifespan.nbinom1.h <- update(just.lifespan.nbinom1, ~ . +  rescale(LargeFROH)*rescale(natal_group_size))

summary(just.lifespan.nbinom1.h) #converges, not sig?? 


#so in summary base model is best
summary(just.lifespan.nbinom1)


#output table
tbl_regression(just.lifespan.nbinom1,
               intercept = T,
               show_single_row = "Sex",
               pvalue_fun = label_style_pvalue(digits = 3),
               estimate_fun = label_style_number(digits = 4),
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            h_countF = "Helper in natal territory",
                            "rescale(mean_total_rain)" = "Mean annual rainfall over lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance in annual rainfall over lifespan",
                            birth_year = "Hatch year",
                            sib_pres = "Sibling presence in nest",
                            "rescale(natal_group_size)" = "Natal group size")
               ) %>%
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) 

#output expected predictions

# Create a dataset with LargeFROH values differing by 0.1
# Hold other variables at means (numeric) or modes (factors)
newdata_low <- data.frame(
  LargeFROH = mean(just_birds$LargeFROH, na.rm = TRUE),
  Sex = 1,  # Mode for factors
  h_countF = 0,
  mean_total_rain = mean(just_birds$mean_total_rain, na.rm = TRUE),
  mean_rain_cv = mean(just_birds$mean_rain_cv, na.rm = TRUE),
  birth_year = mean(just_birds$birth_year, na.rm = TRUE),
  sib_pres = 0,
  natal_group_size = mean(just_birds$natal_group_size, na.rm = TRUE),
  mum = NA,  # Set random effects to NA for population-level predictions
  dad = NA
)

newdata_high <- newdata_low
newdata_high$LargeFROH <- newdata_high$LargeFROH + 0.1

# Predict expected lifespan for both scenarios
pred_low <- predictions(just.lifespan.nbinom1,
                        newdata = newdata_low,
                        component        = "conditional", 
                        allow.new.levels = TRUE
)

pred_high <- predictions(just.lifespan.nbinom1,
                         newdata = newdata_high,
                         type = "response", 
                         allow.new.levels = TRUE)


pred_low <- glmmTMB::predict.glmmTMB(
  object     = just.lifespan.nbinom1,
  newdata    = newdata_low,
  type       = "response",
  component  = "conditional",
  re.form    = NA
)


# Compute the discrete change
change <- pred_high - pred_low
change


# Create a combined dataset
newdata_both <- rbind(newdata_low, newdata_high)
newdata_both$LargeFROH_scenario <- c("low", "high")

# Get predictions with confidence intervals
preds <- predictions(just.lifespan.nbinom1, 
                     newdata = newdata_both, 
                     type = "response", 
                     allow.new.levels = TRUE,
                     component = "dispersion")

# Compute the difference
diff <- diff(preds, hypothesis = "revpairwise")
print(diff)


####plots

##plot with ggpredict
#make ggpredict object
ggpred.just.lifespan.nbinom1.f <-  ggpredict(just.lifespan.nbinom1.f,
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
  geom_point(alpha = 1)+
  geom_smooth(method = "lm",
              colour = "black",
              linewidth = 2)+
  theme_classic2()+
  labs( x=expression(F[ROH]),
        y= "Lifespan")+
  theme(text = element_text(size = 20))+
  ggtitle("")

  
plot(just.lifespan.largefroh.plot)


#------------------------first year survival----------------------

#use binomial family as survival is binary

just.fys.1 <- glmer(unix_fys ~ LargeFROH  + Sex + h_count + birth_total_rain + BirthRainCV + 
                      (1 | birth_year) + (1 | mum) + (1 | dad) ,
                    data = just_birds,
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))

)

summary(just.fys.1)

#model does not converge, try rescaling as suggested

just.fys.rescale <- glmer(unix_fys ~ rescale(LargeFROH) + Sex + h_countF  + 
                            rescale(birth_total_rain) + rescale(BirthRainCV)+
                            sib_pres + rescale(natal_group_size) +
                      (1 | birth_year) + (1 | mum) + (1 | dad) ,
                    data = just_birds,
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))
                    
)

summary(just.fys.rescale)
isSingular(just.fys.rescale)
#fits is singular

just.fys.rescale.a <- update(just.fys.rescale, ~ . - (1| mum) - (1| dad))
isSingular(just.fys.rescale.a) #no longer singular

simulateResiduals(just.fys.rescale.a, plot = T) #looks fine

#check for collinearity and overdispersion
check_collinearity(just.fys.rescale.a) #vif around one
check_overdispersion(just.fys.rescale.a) #no overdisp

summary(just.fys.rescale.a)

###try interactions

#froh * rain

just.fys.rescale.b <- update(just.fys.rescale.a, ~ . + rescale(LargeFROH)*rescale(birth_total_rain))

summary(just.fys.rescale.b)# model converges. int not sig

anova(just.fys.rescale, just.fys.rescale.b) #anova says no

#froh * help

just.fys.rescale.c   <- update(just.fys.rescale.a, ~ . + rescale(LargeFROH)*h_countF)

summary(just.fys.rescale.c)# model converges. int not sig

anova(just.fys.rescale, just.fys.rescale.c) #anova says no

#froh * sex interaction

just.fys.rescale.d <- update(just.fys.rescale.a, ~ . + rescale(LargeFROH)*Sex)

summary(just.fys.rescale.d)# model converges. int not sig

anova(just.fys.rescale, just.fys.rescale.d) #anova says no

#froh * sib_pres interaction

just.fys.rescale.e <- update(just.fys.rescale.a, ~ . + rescale(LargeFROH)*sib_pres)

summary(just.fys.rescale.e)# model converges. int not sig

anova(just.fys.rescale, just.fys.rescale.d) #anova says no

#froh * natal_group_size

#froh * sex interaction

just.fys.rescale.g <- update(just.fys.rescale.a, ~ . + rescale(LargeFROH)*rescale(natal_group_size))

summary(just.fys.rescale.g)# model converges. int not sig

anova(just.fys.rescale, just.fys.rescale.d) #anova says no


####summary

#In summary we keep model a
summary(just.fys.rescale.a)



#output table
tbl_regression(just.fys.rescale.a, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            h_countF = "Helper in natal territory",
                            "rescale(birth_total_rain)"	 = "Annual rainfall in hatch year",
                            "rescale(BirthRainCV)" = "Variance in rainfall in hatch year",
                            "sib_pres" = "Presence of sibling in nest",
                            "rescale(natal_group_size)" = "Natal group size"),
                estimate_fun = label_style_sigfig(digits = 3),
                pvalue_fun = label_style_pvalue(digits = 3)
                            )%>% 

  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) %>%
  
  modify_fmt_fun(
    std.error  = function(x) style_sigfig(x, digits = 3),
    statistic  = function(x) style_sigfig(x, digits = 3)
  )



#plot with ggpredict. doesn't quite work


data(efc, package = "ggeffects")

str(efc)

fit <- lm(barthtot ~ c12hour + neg_c_7 + c161sex + c172code, data = efc)

predict_response(fit, terms = "c12hour")





##plots
#with ggpredict
##plot with ggpredict
#make ggeffect object

just.fys.glm <- glm(unix_fys ~ LargeFROH  + Sex + help + birth_total_rain + BirthRainCV  ,
                    data = just_birds,
                    family=binomial
)
                    
glmobj <- predict_response(just.fys.glm,
                           terms = "LargeFROH[all]")


plot(glmobj)



ggpred.just.fys.rescale.a <-  predict_response(just.fys.rescale.a,
                                               terms = "LargeFROH[all]",
                                               condition=c(natal_group_size =1, sib_pres=1)) 


ggpred.just.fys.rescale.a %>% as.data.frame()

ggplot(ggpred.just.fys.rescale.a, aes(x,predicted)) +geom_line()

plot(ggpred.just.fys.rescale.a)

ggplot(ggpred.just.fys.rescale.a, aes(x = x, y = y)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Rescaled LargeFROH",
    y = "Predicted Probability of unix_fys",
    title = "Logistic Regression Prediction by LargeFROH and sib_pres"
  ) +
  theme_minimal()


#plot line graph with ggplot

ggplot(just_birds, aes(LargeFROH, unix_fys))+
  geom_smooth(method = "lm",
              col = "black",
              size = 2)+
  geom_point(shape = 1,
             size = 2)+
  scale_y_continuous(breaks = c(0,1))+
  labs(title = "Fig. 1",
       y = "First year survival",
       x=expression(F[ROH])) +
  theme_classic()+
  theme(text = element_text(size = 20)) 




#we will now try this with different ROH classes to see what is happening
  
  
just.fys.small <- glmer(unix_fys ~ rescale(RegFROH) + Sex + help  + 
                              rescale(birth_total_rain) + rescale(BirthRainCV)+
                              sib_pres + rescale(natal_group_size) +
                              (1 | birth_year) ,
                            data = just_birds,
                            family=binomial,
                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))
                            
  )

summary(just.fys.small) #still not sig with ROH > 1.6MB


just.fys.F10 <- glmer(unix_fys ~ rescale(F10) + Sex + help  + 
                          rescale(birth_total_rain) + rescale(BirthRainCV)+
                          sib_pres + rescale(natal_group_size) +
                          (1 | birth_year) ,
                        data = just_birds,
                        family=binomial,
                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))
                        
)

summary(just.fys.F10) #still not significant with RZooRoH derived FROH
  
  
  
#----------------------------------lifetime reproductive success--------------  

hist(just_birds$n_off, breaks = 25)
  
just.lrs.poisson <- glmmTMB(n_off ~  rescale(LargeFROH) + Sex + h_countF  + 
                              rescale(mean_total_rain) + rescale(mean_rain_cv)+
                              sib_pres + rescale(natal_group_size) + birth_year +
                              lifespan +
                              (1 | birth_year) + (1 | mum) + (1 | dad),
                      data = just_birds,
                      ziformula=~0,
                      family=poisson)


summary(just.lrs.poisson)

simulateResiduals(just.lrs.poisson, plot = T) #a bit wiggly there

#check for collinearity and overdispersion
check_collinearity(just.lrs.poisson) #vif around one
check_overdispersion(just.lrs.poisson) #overdispersion not bad
check_zeroinflation(just.lrs.poisson)
testZeroInflation(just.lrs.poisson)

##try adding zero inflation

just.lrs.poisson <- glmmTMB(n_off ~ rescale(LargeFROH) + Sex + h_countF  + 
                              rescale(mean_total_rain) + rescale(mean_rain_cv)+
                              sib_pres + rescale(natal_group_size) + birth_year +
                              lifespan +
                              (1 | birth_year) +
                              (1 | mum) + (1 | dad),
                            data = just_birds,
                            ziformula=~rescale(LargeFROH) + Sex + h_countF  + 
                              rescale(mean_total_rain) + rescale(mean_rain_cv)+
                              sib_pres + rescale(natal_group_size) + birth_year +
                              lifespan,
                            family=poisson)

check_overdispersion(just.lrs.poisson)
simulateResiduals(just.lrs.poisson, plot = T) #much better
check_zeroinflation(just.lrs.poisson) #still underfitting, try nbinom
check_collinearity(just.lifespan.poisson)

#nbinom1 
just.lrs.nbinom1 <- update(just.lrs.poisson, family=nbinom1())
simulateResiduals(just.lrs.nbinom1, plot = TRUE) #qq is better, but residuals are slightly bad

#check for collinearity and overdispersion
check_collinearity(just.lrs.nbinom1) #vif around one
check_overdispersion(just.lrs.nbinom1) # no overdispersion
check_zeroinflation(just.lrs.nbinom1)

#nbinom2
just.lrs.nbinom2 <- update(just.lrs.poisson, family=nbinom2())

simulateResiduals(just.lrs.nbinom2, plot = TRUE)
#qq is slightly worse, but residuals are slightly bad

#check for collinearity and overdispersion
check_collinearity(just.lrs.nbinom2) #vif around one
check_overdispersion(just.lrs.nbinom2) #no overdisp
check_zeroinflation(just.lrs.nbinom2)


#AIC comparisons
AIC(just.lrs.nbinom1, just.lrs.nbinom2, just.lrs.poisson)

#df      AIC
#                 df      AIC
#just.lrs.nbinom1 13 4206.997
#just.lrs.nbinom2 13 4214.194
#just.lrs.poisson 12 4354.147

summary(just.lrs.nbinom1)
summary(just.lrs.nbinom2)
summary(just.lrs.poisson)

#nbinoms have better AIC but are unstable. lets use poisson.



#try interactions

#froh*help

just.lrs.poisson.b <- update(just.lrs.poisson , ~ . + rescale(LargeFROH)*h_countF)

summary(just.lrs.poisson.b) #model converges,  sig

anova(just.lrs.poisson, just.lrs.poisson.b) #anova says not better fit

#froh*sex

just.lrs.poisson.c <- update(just.lrs.poisson , ~ . + rescale(LargeFROH)*Sex)

summary(just.lrs.poisson.c) #model has trouble. interaction marginal

anova(just.lrs.poisson, just.lrs.poisson.c) #anova says no

#froh*mean_rain

just.lrs.poisson.d <- update(just.lrs.poisson , ~ . + rescale(LargeFROH)*rescale(mean_total_rain))

summary(just.lrs.poisson.d) #not significant

anova(just.lrs.poisson, just.lrs.poisson.d) #anova says no!?


#froh*rain variance

just.lrs.poisson.e <- update(just.lrs.poisson , ~ . + rescale(LargeFROH)*rescale(mean_rain_cv))

summary(just.lrs.poisson.e) #lots of failures of convergence here

anova(just.lrs.poisson, just.lrs.poisson.e) #anova says no!


#froh*sib_pres

just.lrs.poisson.f <- update(just.lrs.poisson , ~ . + rescale(LargeFROH)*sib_pres)

summary(just.lrs.poisson.f) #not sig

anova(just.lrs.poisson, just.lrs.poisson.f) #anova says no!!!


#froh*natal_group_size

just.lrs.poisson.g <- update(just.lrs.poisson , ~ . + rescale(LargeFROH)*rescale(natal_group_size))

summary(just.lrs.poisson.g) #nope

anova(just.lrs.poisson, just.lrs.poisson.g) #anova says no



#lets check these hold true with other roh types

just.lrs.poisson.reg <- glmmTMB(n_off ~ rescale(RegFROH) + Sex + h_countF  + 
                              rescale(mean_total_rain) + rescale(mean_rain_cv)+
                              sib_pres + rescale(natal_group_size) +
                              (1 | birth_year) +
                              (1 | mum) + (1 | dad),
                            data = just_birds,
                            ziformula=~1,
                            family=poisson)

summary(just.lrs.poisson.reg) #same same

#with RZooROH
just.lrs.poisson.F10 <- glmmTMB(n_off ~  rescale(F10) + Sex + h_countF  + 
                                  rescale(birth_total_rain) + rescale(BirthRainCV)+
                                  sib_pres + rescale(natal_group_size) +
                                  lifespan+
                                  (1 | birth_year) +
                                  (1 | mum) + (1 | dad),
                                data = just_birds,
                                ziformula=~0,
                                family=poisson)


summary(just.lrs.poisson.F10) #same same



####summary
#in summary no interactions are significantly better so we go with base model

summary(just.lrs.poisson)

just.lrs.poisson$explained_variance$rescale(LargeFROH)

explained_variance(just.lrs.poisson$rescale(LargeFROH))

tbl_regression(just.lrs.poisson, intercept = T,
               show_single_row = "Sex",
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            h_countF = "Helper in natal territory",
                            "rescale(mean_total_rain)" = "Mean annual rainfall over lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance in annual rainfall over lifespan",
                            birth_year = "Hatch year",
                            sib_pres = "Sibling presence in nest",
                            "rescale(natal_group_size)" = "Natal group size",
                            lifespan = "Lifespan"),
               estimate_fun = label_style_sigfig(digits = 3),
               pvalue_fun = label_style_pvalue(digits = 3)
               
)%>%
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic))  %>%
  
  modify_fmt_fun(
         std.error  = function(x) style_sigfig(x, digits = 3),
        statistic  = function(x) style_sigfig(x, digits = 3)
  )


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


#do this plot again with a slightly trimmed dataset to make it look nice


LRS.plot.model <- glmmTMB(n_off ~ LargeFROH + help + Sex + 
                              mean_total_rain + mean_rain_cv + birth_year +
                              (1 | mum) + (1 | dad),
                            data = filter(just_birds, n_off < 16),
                            ziformula=~1,
                            family=poisson)


LRS.plot.ggpred <-  ggpredict(LRS.plot.model,
                                      terms = "LargeFROH[all]")


autoplot(LRS.plot.ggpred)+
  geom_expected_line(linewidth = 1, color = "red")+
  geom_CI_ribbon(fill = "red")+
  theme_classic2()+
  labs( x = "FROH > 3.3Mb", y= "Lifetime reproductive success")+
  theme(text = element_text(size = 40)) +
  ggtitle("         Individual inbreeding")+
  layer_fit_data(alpha = 0.2,color = "red")


#plot with ggplot

#ggplot
just.lrs.largefroh.plot <- ggplot(just_birds, aes(LargeFROH, n_off)) +
  geom_point(alpha = 1)+
  geom_smooth(method = "lm" , 
              colour = "black", 
              linewidth = 2)+
  theme_classic2()+
  labs( x = expression(F[ROH]),
        y= "Productivity")+
  theme(text = element_text(size = 20))+
  ggtitle("")


plot(just.lrs.largefroh.plot)

#plot for poster design

poster.plot <- ggplot(just_birds, aes(LargeFROH, lifespan)) +
  geom_point(colour = "red")+
  geom_smooth(method = "lm" , colour = "red", size = 2,fill = "red")+
  theme_classic2()+
  labs( x = "FROH > 3.3Mb", y= "Lifetime reproductive success")+
  theme(text = element_text(size = 40),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Individual inbreeding depression")

poster.plot

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
stat.cox.me.a <- coxme(Surv(lifespan, event) ~  rescale(LargeFROH) + Sex + help  +
                       sib_pres + natal_group_size + MumAge + mean_total_rain +  
                         (1 | mum) + (1 | dad) + (1 | birth_year),
                       data = stat.birds , 
                       control = coxme.control(iter.max  = 200)
) 

cox.zph(stat.cox.me.a) #GLOBAL now conforms to ph assumptions

summary(stat.cox.me.a)

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




#----------------------EPP distribution---------------

#check if inbred social fathers are more likely to be cucked.
#i.e. do individuals that are EPP have higher social father FROH

hist(just_birds$BrMLargeFroh)

boxplot(BrMLargeFroh ~ EPP, data = just_birds) #do epp offspring have higher social father froh? 

wilcox.test(BrMLargeFroh ~ EPP, data = just_birds) #no

#do epp offspring have higher social father froh than genetic father froh?
wilcox.test(filter(just_birds, EPP == 1)$BrMLargeFroh, filter(just_birds, EPP == 1)$dadLargeFroh )
#no


  
dadrohs <- c(just_birds$BrMLargeFroh,just_birds$dadLargeFroh)

dadrohlabels <- c(rep("Social", length(just_birds$BrMLargeFroh)), rep("Genetic", length(just_birds$BrMLargeFroh)))

boxplot_dadrohs <- data.frame(dadrohs, dadrohlabels)

boxplot(dadrohs ~ dadrohlabels, data = boxplot_dadrohs)

#---inbreeding load------------

#regression of the natural logarithm of *trait* ~ inbreeding gives 'inbreeding load'

#as per doi: 10.1111/eva.12713, use a GLM with link = log


#genetic load for whole lifespan

load_model <- glm(lifespan ~ LargeFROH + birth_year + mean_total_rain,
                  family = poisson(link = "log"),
                  data = just_birds)


summary(load_model)

library(sandwich) # library to get robust estimators
library(lmtest) # library to test coefficients

#get robust SE
( inb.se <- coeftest(load_model, vcov = sandwich) )

#find 95% confidence limits using 1.96 * SE
c( -(inb.se["LargeFROH", "Estimate"] + 1.96*inb.se["LargeFROH", "Std. Error"]),
   -inb.se["LargeFROH", "Estimate"],
   -(inb.se["LargeFROH", "Estimate"] - 1.96*inb.se["LargeFROH", "Std. Error"]) )



just_birds$m



#---------------------------Animal model tests---------------------

#use bayesian animal models to verify if adding a grm does anything

#install packages

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(MCMCglmm)
library(lqmm)
library(Matrix)
library(readxl)
library(pedigree)


#script taken from kiran lee

#Load pedigree of combined masterbayes and sequioa Pedigree

doubleped <- read.csv("combined_pedigree.csv") %>% select(BirdID, mum, dad)
str(doubleped)

# Step 1: Identify parents not in BirdID
unique_parents <- unique(c(doubleped$mum, doubleped$dad))
 
str(unique_parents)

missing_parents <- unique_parents[!unique_parents %in% doubleped$BirdID]
str(missing_parents)


# Step 2: Create new rows for missing parents
new_rows <- data.frame(
  BirdID = missing_parents,
  mum = NA,
  dad = NA
)

# Step 3: Combine with original data
expanded_ped <- rbind(doubleped, new_rows)

# Step 4: Verify
str(expanded_ped)
tail(expanded_ped)  # Check the tail for added parents

ord<- orderPed(expanded_ped) #So offspring follow parents

ordered_ped <- expanded_ped[order(ord),]
head(ordered_ped)
ordered_ped <- ordered_ped[-1,]

str(ordered_ped)

inverseAmatrix <- inverseA(pedigree = ordered_ped)$Ainv # convert pedigree to invert relatedness matrix
dim(inverseAmatrix)

#add to main datatable

animal_birds <- just_birds %>%
  mutate(animal = factor(BirdID, levels = rownames(inverseAmatrix))) %>%
  filter(                                          # remove some NAs
    !is.na(mean_total_rain),
    !is.na(mean_rain_cv),
    !is.na(mum),
    !is.na(dad),
    !is.na(animal),
  ) 

str(animal_birds)

#do the model

animal_lifespan.1 <- MCMCglmm(
                      fixed =  lifespan ~ LargeFROH + Sex  + mean_total_rain + mean_rain_cv , #Response and Fixed effect formula
                     random = ~ mum + dad + birth_year + animal, # Random effect formula
                     ginverse = list(animal = inverseAmatrix), # correlations among random effect levels (here breeding values)
                     data = animal_birds)# data set


summary(animal_lifespan.1)

plot(animal_lifespan.1$VCV)


var_factor_lifespan <- animal_lifespan.1$VCV[, "animal"]
var_residual_lifespan <- animal_lifespan.1$VCV[, "units"]

# Compute proportion of variance
prop_variance_lifespan <- var_factor_lifespan / (var_factor_lifespan + var_residual_lifespan)

# Summarise
mean(prop_variance_lifespan)

summarise(data.frame(c = prop_variance_lifespan))




#do the same for lifetime reproductive success

animal_LRS.1 <- MCMCglmm(
  fixed =  n_off ~ LargeFROH + Sex + h_countF + mean_total_rain + mean_rain_cv , #Response and Fixed effect formula
  random = ~ mum + dad + birth_year + animal, # Random effect formula
  ginverse = list(animal = inverseAmatrix), # correlations among random effect levels (here breeding values)
  data = animal_birds)# data set

summary(animal_LRS.1)


var_factor_LRS <- animal_LRS.1$VCV[, "animal"]
var_residual_LRS <- animal_LRS.1$VCV[, "units"]

# Compute proportion of variance
prop_variance_LRS <- var_factor_LRS / (var_factor_LRS + var_residual_LRS)

# Summarise
t.test(prop_variance_LRS)

summarise(data.frame(c = prop_variance_LRS))



#animal model FYS

animal.FYS.1 <- MCMCglmm(
  fixed =  unix_fys ~ LargeFROH + Sex +  h_countF+ mean_total_rain + mean_rain_cv , #Response and Fixed effect formula
  random = ~ mum + dad +animal, # Random effect formula
  ginverse = list(animal = inverseAmatrix), # correlations among random effect levels (here breeding values)
  data = animal_birds)# data set

summary(animal.FYS.1)



var_factor_FYS <- animal.FYS.1$VCV[, "animal"]
var_residual_FYS <- animal.FYS.1$VCV[, "units"]

# Compute proportion of variance
prop_variance_LRS <- var_factor_LRS / (var_factor_LRS + var_residual_LRS)

# Summarise
t.test(prop_variance_LRS)

