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

#---------------------import data--------------

##Read in data
stat.birds <- read.csv("data/clean_data.csv")
str(stat.birds)

#Factorise 
stat.birds$h_countF <- as.factor(stat.birds$h_count)
stat.birds$BreedSeason <- as.factor(stat.birds$BreedSeason)
stat.birds$Sex <- as.factor(stat.birds$Sex)
stat.birds$mum <- as.factor(stat.birds$mum)
stat.birds$dad <- as.factor(stat.birds$dad)
stat.birds$birth_yearF <- as.factor(stat.birds$birth_year)
stat.birds$adulthood <- as.factor(stat.birds$adulthood)
stat.birds$sib_pres <- as.factor(stat.birds$sib_pres)


dim(stat.birds)  #n = 1210


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
#simplify random effects structure to avoid singularity

just.fys.rescale.a <- update(just.fys.rescale, ~ . - (1| mum) - (1| dad))
isSingular(just.fys.rescale.a) #no longer singular

simulateResiduals(just.fys.rescale.a, plot = T) #looks fine

#check for collinearity and overdispersion
check_collinearity(just.fys.rescale.a) #vif around one
check_overdispersion(just.fys.rescale.a) #no overdispersion

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

#correlation r squared

r2(just.fys.rescale.a)

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



##plots

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


just.fys.small <- glmer(unix_fys ~ rescale(RegFROH) + Sex   + 
                          rescale(birth_total_rain) + rescale(BirthRainCV)+
                          sib_pres + rescale(natal_group_size) +
                          (1 | birth_year) ,
                        data = just_birds,
                        family=binomial,
                        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))
                        
)

summary(just.fys.small) #still not sig with ROH > 1.6MB


just.fys.F10 <- glmer(unix_fys ~ rescale(F10) + Sex +  
                        rescale(birth_total_rain) + rescale(BirthRainCV)+
                        sib_pres + rescale(natal_group_size) +
                        (1 | birth_year) ,
                      data = just_birds,
                      family=binomial,
                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e8))
                      
)

summary(just.fys.F10) #still not significant with RZooRoH derived FROH



#---------------------lifespan models------------------

hist(just_birds$lifespan, breaks = 29)


just.lifespan.nbinom2 <- glmmTMB(
  lifespan ~ rescale(LargeFROH) + Sex + h_countF +
    rescale(mean_total_rain) + rescale(mean_rain_cv) +
    birth_year + sib_pres + rescale(natal_group_size) ,
  data = just_birds %>% filter(lifespan > 0),
  ziformula =  ~ 0,
  family = nbinom2
)

simulateResiduals(just.lifespan.nbinom2, plot = TRUE)
#qq plot looks ok
check_overdispersion(just.lifespan.nbinom2) #no overdispersion
check_singularity(just.lifespan.nbinom2)    #is NOT singular
check_convergence(just.lifespan.nbinom2) #it converges!

summary(just.lifespan.nbinom2)

#quick plot
ggplot(just_birds, aes(LargeFROH, lifespan)) +
  geom_point()+
  geom_smooth(method = "lm")

##now check for interactions

#FROH * sex interaction

just.lifespan.nbinom2.c <- update(just.lifespan.nbinom2, ~ . +  rescale(LargeFROH)*Sex)

summary(just.lifespan.nbinom2.c) #model converges. int not sig

anova(just.lifespan.nbinom2, just.lifespan.nbinom2.c) #anova says 0.8, not useful to include

#FROH * mean_rain interaction

just.lifespan.nbinom2.d <- update(just.lifespan.nbinom2, ~ . +  rescale(LargeFROH)*rescale(mean_total_rain))

summary(just.lifespan.nbinom2.d) #model converges. int not sig

anova(just.lifespan.nbinom2, just.lifespan.nbinom2.d) #anova says 0.55, not useful to include


#FROH * rain variance

just.lifespan.nbinom2.e <- update(just.lifespan.nbinom2, ~ . +  rescale(LargeFROH)*rescale(mean_rain_cv))

summary(just.lifespan.nbinom2.e) #model converges. int not sig

anova(just.lifespan.nbinom2, just.lifespan.nbinom2.e) #anova says 0.74, not useful to include


#FROH * Help interaction

just.lifespan.nbinom2.f <- update(just.lifespan.nbinom2, ~ . +  rescale(LargeFROH)*h_countF)

summary(just.lifespan.nbinom2.f) #model converges. F1 almost significant

print(anova(just.lifespan.nbinom2, just.lifespan.nbinom2.f)) #anova says no


#FROH * sib_pres interaction

just.lifespan.nbinom2.g <- update(just.lifespan.nbinom2.f, ~ . +  rescale(LargeFROH)*sib_pres)

summary(just.lifespan.nbinom2.g) #model converges. not sig

anova(just.lifespan.nbinom2.f, just.lifespan.nbinom2.g) #anova says no


#FROH * natal group size interaction

just.lifespan.nbinom2.h <- update(just.lifespan.nbinom2, ~ . +  rescale(LargeFROH)*rescale(natal_group_size))

summary(just.lifespan.nbinom2.h) #converges, not sig?? 


#so in summary base model is best
summary(just.lifespan.nbinom2)


#get r2 

r.squaredGLMM(just.lifespan.nbinom2)



#output table
tbl_regression(just.lifespan.nbinom2,
               intercept = T,
               show_single_row = "Sex",
               pvalue_fun = label_style_pvalue(digits = 3),
               estimate_fun = label_style_number(digits = 3),
               label = list("rescale(LargeFROH)" = "FROH > 3.3Mb", 
                            h_countF = "Helper in natal territory",
                            "rescale(mean_total_rain)" = "Mean annual rainfall over lifespan",
                            "rescale(mean_rain_cv)" = "Mean variance in annual rainfall over lifespan",
                            birth_year = "Hatch year",
                            sib_pres = "Sibling presence in nest",
                            "rescale(natal_group_size)" = "Natal group size"),

)%>% 
  
  modify_column_hide(column = conf.low) %>%
  
  modify_column_unhide(column = c(std.error,statistic)) %>%
  
  modify_fmt_fun(
    std.error  = function(x) style_sigfig(x, digits = 3),
    statistic  = function(x) style_sigfig(x, digits = 3)
  )


#plot


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

simulateResiduals(just.lrs.poisson, plot = T) #slightly off centre

#check for collinearity and overdispersion
check_collinearity(just.lrs.poisson) #vif around one
check_overdispersion(just.lrs.poisson) #overdispersion not bad
check_zeroinflation(just.lrs.poisson)
testZeroInflation(just.lrs.poisson) #model is underfitting zeros

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
check_zeroinflation(just.lrs.poisson) #good obs/expected zero ratio
check_collinearity(just.lifespan.poisson) #vifs low


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

summary(just.lrs.poisson.reg) #similar results

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


summary(just.lrs.poisson.F10) #similar



####summary
#in summary no interactions are significantly better so we go with base model

summary(just.lrs.poisson)


r2(just.lrs.poisson)


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



#do this using the actual model

# get the slope of the scaled variable
beta_scaled <-  tidy(just.lifespan.poisson)[2,4:5]

# compute SD of un-rescaled FROH
sd_FROH  <- sd(just_birds$LargeFROH, na.rm=TRUE)


#  Convert the slope back to the unscaled FROH
beta_original <- beta_scaled / (2 * sd_FROH)

#times 2 for LE per diploid genome

LE <- beta_original * 2

LE




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

