##################################################################################################################################
# 1.Descriptive statistics
#~~~~~~~~~~~~~~~~~~~~~~~~~
summary(AAD)

mean(AAD$D1.Shannon.diversity,na.rm=TRUE)
median(AAD$D1.Shannon.diversity,na.rm=TRUE)
min(AAD$D1.Shannon.diversity,na.rm=TRUE)
max(AAD$D1.Shannon.diversity,na.rm=TRUE)
quantile(AAD$D1.Shannon.diversity,na.rm=TRUE,c(0.25,0.75))

mean(AAD$D6.Shannon.diversity,na.rm=TRUE)
median(AAD$D6.Shannon.diversity,na.rm=TRUE)
min(AAD$D6.Shannon.diversity,na.rm=TRUE)
max(AAD$D6.Shannon.diversity,na.rm=TRUE)
quantile(AAD$D6.Shannon.diversity,na.rm=TRUE,c(0.25,0.75))

mean(AAD$D1.Chao1.diversity,na.rm=TRUE)
median(AAD$D1.Chao1.diversity,na.rm=TRUE)
min(AAD$D1.Chao1.diversity,na.rm=TRUE)
max(AAD$D1.Chao1.diversity,na.rm=TRUE)
quantile(AAD$D1.Chao1.diversity,na.rm=TRUE,c(0.25,0.75))

mean(AAD$D6.Chao1.diversity,na.rm=TRUE)
median(AAD$D6.Chao1.diversity,na.rm=TRUE)
min(AAD$D6.Chao1.diversity,na.rm=TRUE)
max(AAD$D6.Chao1.diversity,na.rm=TRUE)
quantile(AAD$D6.Chao1.diversity,na.rm=TRUE,c(0.25,0.75))

mean(AAD$D1.D6.Jaccard.distance,na.rm=TRUE)
median(AAD$D1.D6.Jaccard.distance,na.rm=TRUE)
min(AAD$D1.D6.Jaccard.distance,na.rm=TRUE)
max(AAD$D1.D6.Jaccard.distance,na.rm=TRUE)
quantile(AAD$D1.D6.Jaccard.distance,na.rm=TRUE,c(0.25,0.75))

table(AAD$Antibiotic.class)
table(AAD$Outcome)

cor(AAD$D1.Shannon.diversity,AAD$D6.Shannon.diversity, use="complete.obs")
cor(AAD$D1.Chao1.diversity,AAD$D6.Chao1.diversity, use="complete.obs")

##################################################################################################################################
# 2.Graphics
#~~~~~~~~~~~~

#1st: Bar chart
barplot(table(AAD$Outcome), xlab="Outcome",ylab="Frequency", col = c("cadetblue1","cadetblue2","cadetblue3"),main="Outcome barchart")
        #ND is the most frequent in data

#2nd: Bar chart
barplot(tapply(AAD$D1.D6.Jaccard.distance, list(Antibiotic=AAD$Antibiotic.class),mean,na.rm=T), xlab="Antibiotics",ylab="Mean jaccard.distance",
        main="Mean_Jaccard.distance barchart", col = c("lightpink","lightpink3","#8B5F65"))

barplot(tapply(AAD$D6.Shannon.diversity, list(Antibiotic=AAD$Antibiotic.class),mean,na.rm=T), xlab="Antibiotics",ylab="Mean D6.Shannon.diversity",
        main="Mean_D6.Shannon.diversity", col = c("purple","purple3","purple4"))

#3rd: Histogram
hist(AAD$D1.Shannon.diversity,xlab="D1 Shannon",main="Distribution of D1 Shannon",col = "aquamarine2")
hist(AAD$D6.Shannon.diversity,xlab="D6 Shannon",main="Distribution of D6 Shannon",col = "yellow2")
        #both of them is left slewed 

#4th: Scattar plot
plot(D6.Shannon.diversity[Antibiotic.class=="FQN"]~D1.Shannon.diversity[Antibiotic.class=="FQN"],data = AAD,xlab="D1 Shannon",ylab="D6 Shannon",col="red",main="Scatterplot of D1_Shannon & D6_Shannon")
points(D6.Shannon.diversity[Antibiotic.class=="OBL"]~ D1.Shannon.diversity[Antibiotic.class=="OBL"],data = AAD,xlab="D1 Shannon",ylab="D6 Shannon",col="darkblue")
points(D6.Shannon.diversity[Antibiotic.class=="PBL"]~ D1.Shannon.diversity[Antibiotic.class=="PBL"],data = AAD,xlab="D1 Shannon",ylab="D6 Shannon",col="chartreuse3")
legend("topleft", legend = c("FQN","OBL","PBL"), fill=c("red","darkblue","chartreuse3"))
#regression lines
abline(lm(AAD$D6.Shannon.diversity[AAD$Antibiotic.class=="FQN"]~AAD$D1.Shannon.diversity[AAD$Antibiotic.class=="FQN"]),col="red")
abline(lm(AAD$D6.Shannon.diversity[AAD$Antibiotic.class=="OBL"]~AAD$D1.Shannon.diversity[AAD$Antibiotic.class=="OBL"]),col="darkblue")
abline(lm(AAD$D6.Shannon.diversity[AAD$Antibiotic.class=="PBL"]~AAD$D1.Shannon.diversity[AAD$Antibiotic.class=="PBL"]),col="chartreuse3")
        #the relashionship between them is not linear
        
#5th: Box plot
boxplot(D1.D6.Jaccard.distance~Antibiotic.class,main="Boxplot for Antibiotics & jaccard",col = c("darkolivegreen2","darkolivegreen4","darkolivegreen"),xlab="Antibiotics",data = AAD)
        #PBL differ from others 

##################################################################################################################################
#3. OUTLIERS
#~~~~~~~~~~~~

boxplot(AAD$D1.Shannon.diversity)
boxplot(AAD$D1.Shannon.diversity)$out

boxplot(AAD$D6.Shannon.diversity)
boxplot(AAD$D6.Shannon.diversity)$out       #The most one contain outliers

boxplot(AAD$D1.Chao1.diversity)
boxplot(AAD$D1.Chao1.diversity)$out     

boxplot(AAD$D6.Chao1.diversity)
boxplot(AAD$D6.Chao1.diversity)$out        

boxplot(AAD$D1.D6.Jaccard.distance)
boxplot(AAD$D1.D6.Jaccard.distance)$out     #Has no outliers


##################################################################################################################################
#4. normality & homoscedasticity 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#4.1.testing for normality
#~~~~~~~~~~~~~~~~~~~~~~~~~~
#method 1: QQ plots & method 2: shapiro test , we applied normality test on all numeric variables in dataset

#1- D1.Shannon.diversity
qqnorm(AAD$D1.Shannon.diversity, ylab="D1 Shannon Diversity", main="Q-Q of D1 Shannon")
qqline(AAD$D1.Shannon.diversity,ylab="D1 Shannon Diversity")
#from the qq plot we can see the d1 shannon diversity data is not normal as the points are further away from the reference line

shapiro.test(AAD$D1.Shannon.diversity) #p-value = 4.134e-13 < 0.05 there is significant difference (not normal)
#from the shapiro test, we can see the p-value is less than 0.05 so we can reject the assumption of normality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2- D6.Shannon.diversity
qqnorm(AAD$D6.Shannon.diversity,ylab="D6 Shannon Diversity", main="Q-Q of D6 Shannon")
qqline(AAD$D6.Shannon.diversity, ylab="D6 Shannon Diversity")   #from the qq plot we can see points did not fit lines so not normal

shapiro.test(AAD$D6.Shannon.diversity) #(p-value = 8.795e-13) < 0.05 
#from shapiro the p value is less than 0.05 so we can reject the assumption of normality

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3- D1.Chao1.diversity
hist(AAD$D1.Chao1.diversity)  #histogram is skewed not normal
qqnorm(AAD$D1.Chao1.diversity,ylab="D1 Chao1 Diversity",main="Q-Q of D1 Chao")    #from the qq plot we can see points did not fit lines so not normal
qqline(AAD$D1.Chao1.diversity)

shapiro.test(AAD$D1.Chao1.diversity) #(p-value = 2.489e-09) < 0.05
#from shapiro the p value < 0.05 so we can reject the assumption of normality for d1 chao1 diversity 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4- D6.Chao1.diversity
qqnorm(AAD$D6.Chao1.diversity,ylab="D6 Chao Diversity",main="Q-Q of D6 Chao")
qqline(AAD$D6.Chao1.diversity)                           #from the qq plot we can see points did not fit lines so not normal

shapiro.test(AAD$D6.Chao1.diversity)  #(p-value = 1.776e-06) < 0.05
#from shapiro the p value < 0.05 so we can reject assumption of normality 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#5- D1.D6.Jaccard.distance
qqnorm(AAD$D1.D6.Jaccard.distance, ylab="Jaccard Distance D1 & D6",main="Jaccard Distance D1 & D6")
qqline(AAD$D1.D6.Jaccard.distance)
shapiro.test(AAD$D1.D6.Jaccard.distance) #(p-value = 8.022e-05) < 0.05
#from shapiro the p value < 0.05 so we can reject assumption of normality 

#All Data variables are not normal according to tests

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.2.testing for homoscedasticity
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Graphical visualizing using box plots 
#method 1: levene test & method 2: barlett test , we applied normality test on all numeric variables in dataset
library(car)
#Note: levene test as is robust to non-normality 

#1- D1.Shannon.diversity

boxplot(AAD$D1.Shannon.diversity,ylab="D1 Shannon values",main="D1 Shannon Data") #

leveneTest(D1.Shannon.diversity ~ Antibiotic.class, data=AAD)  #p value 0.857 > 0.05 (homoscedastic)

bartlett.test(D1.Shannon.diversity ~ Antibiotic.class, data = AAD)  #p-value = 0.9746  > 0.05 so assume homoscedasticity
#From 2 methods D1.Shannon.diversity so we have enough evidence to accept assumption that it is homoscedastic

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2- D6.Shannon.diversity
boxplot(AAD$D6.Shannon.diversity,ylab="D6 Shannon values",main="D6 Shannon Data")

leveneTest(D6.Shannon.diversity~Antibiotic.class,data=AAD) #p value= 0.4289 > 0.05 assume homoscedasticity

bartlett.test(D6.Shannon.diversity ~ Antibiotic.class, data = AAD)  #p-value = 0.1757 > 0.05 so we have enough evidence to accept assumption of homoscedasticity
#From 3 methods D6.Shannon.diversity so we have enough evidence to accept assumption that it is homoscedastic
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3- D1.Chao1.diversity
boxplot(AAD$D1.Chao1.diversity,ylab="D1 Chao values",main="D1 Chao Data")

leveneTest(D1.Chao1.diversity~Antibiotic.class,data=AAD) # p value 0.9278 > 0.05 assume homoscedasticity

bartlett.test(D1.Chao1.diversity ~ Antibiotic.class, data = AAD) #p-value = 0.5621 > 0.05
#since p-value is greater than 0.05 so we fail to reject the assumption of homoscedasticity
#From 2 methods D1.Chao1.diversity so we have enough evidence to accept assumption that it is homoscedastic

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4- D6.Chao1.diversity

boxplot(AAD$D6.Chao1.diversity,ylab="D6 Chao values",main="D6 Chao Data")

leveneTest(D6.Chao1.diversity~Antibiotic.class,data=AAD) # p value 0.295 > 0.05 assume homoscedasticity

bartlett.test(D6.Chao1.diversity ~ Antibiotic.class, data = AAD) # p-value = 0.206 > 0.05
#since p-value is greater than 0.05 so we fail to reject the assumption of homoscedasticity
#From 2 methods D6.Chao1.diversity so we have enough evidence to accept assumption that it is homoscedastic

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#5- D1.D6.Jaccard.distance

boxplot(AAD$D1.D6.Jaccard.distance,ylab="D1 D6 Jaccard values",main="D1 D6 Jaccard") #no outliers

leveneTest(D1.D6.Jaccard.distance ~ Antibiotic.class,data=AAD) #p value 0.02891 < 0.05 is significant so we fail to assume homoscedasticity

bartlett.test(D1.D6.Jaccard.distance ~ Antibiotic.class, data=AAD) #p-value = 0.1448 > 0.05 assume homoscedasticity
#but levene's test is more reliable as its robust against non-normality in the data
#so we still fail to assume homoscedatsicity

##################################################################################################################################
#5. statistical inference
#~~~~~~~~~~~~~~~~~~~~~~~~~
#note -> lm: Calculate the mean and standard error

antibiotics = unique(AAD$Antibiotic.class)
for (antibiotic in antibiotics){
  level_90 = confint(lm(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class == antibiotic] ~ 1, AAD), level=0.90)
  print(antibiotic)
  cat("level 0.90 of ", antibiotic, level_90,"\n")
}

for (antibiotic in antibiotics){
  level_95 = confint(lm(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class == antibiotic] ~ 1, AAD), level=0.95)
  print(antibiotic)
  cat("level 0.95 of ", antibiotic, level_95,"\n")
}

for (antibiotic in antibiotics){
  level_99 = confint(lm(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class == antibiotic] ~ 1, AAD), level=0.99)
  print(antibiotic)
  cat("level 0.99 of ", antibiotic, level_99,"\n")
}
#as level increased from 0.90 to 0.99 confidence intervals got much greater 
#Note: The wider the confidence level, the wider the interval becomes, reflecting a higher level of confidence in capturing the true value.

##################################################################################################################################
#6.	Hypothesis testing
#~~~~~~~~~

#1st point

#We hypothesis that Chao/Shannon at day 6 different between CDI vs ND

#null hypothesis: there is no difference between D6.chao/shannon between CDI and ND 
#alternative: there is a difference between D6.chao/shannon between CDI and ND 

#Assuming normality and homoscedasticity we use t-test
t.test(AAD$D6.Chao1.diversity[AAD$Outcome=="CDI"],AAD$D6.Chao1.diversity[AAD$Outcome=="ND"], var.equal= TRUE)
# p-value = 0.3362 (not significant so we fail to reject the null ) so we can say there is no difference 
t.test(AAD$D6.Shannon.diversity[AAD$Outcome=="CDI"],AAD$D6.Shannon.diversity[AAD$Outcome=="ND"], var.equal= TRUE)
# p-value = 0.7164 (not significant so we fail to reject the null ) so we can say there is no difference


# assessing normality assumption
shapiro.test(AAD$D6.Chao1.diversity[AAD$Outcome=="CDI"])   #p-value = 0.04323 assumption of normality is rejected
shapiro.test(AAD$D6.Chao1.diversity[AAD$Outcome=="ND"])    #p-value = 1.831e-05 assumption of normality is rejected

shapiro.test(AAD$D6.Shannon.diversity[AAD$Outcome=="CDI"]) #p-value = 0.1407 #normal
shapiro.test(AAD$D6.Shannon.diversity[AAD$Outcome=="ND"])  #p-value = 6.086e-12


#assessing homoscedasticity 
boxplot(AAD$D6.Chao1.diversity[AAD$Outcome=="CDI"],AAD$D6.Chao1.diversity[AAD$Outcome=="ND"],names = c("Chao6 CDI", "Chao6 ND"))     #variances differ in the two groups
boxplot(AAD$D6.Shannon.diversity[AAD$Outcome=="CDI"],AAD$D6.Shannon.diversity[AAD$Outcome=="ND"],names = c("Shannon6 CDI", "Shannon6 ND")) #variances differ in the two groups


# since the two assumptions are not met we can use a non-parametric test
wilcox.test(AAD$D6.Chao1.diversity[AAD$Outcome=="CDI"],AAD$D6.Chao1.diversity[AAD$Outcome=="ND"])
# p-value = 0.2786 (not significant so we fail to reject the null) so we can say there is no difference 
wilcox.test(AAD$D6.Shannon.diversity[AAD$Outcome=="CDI"],AAD$D6.Shannon.diversity[AAD$Outcome=="ND"])
#p-value = 0.8014 (not significant so we fail to reject the null) so we can say there is no difference 

#~~~~~~~~~~~~~~~~~~~
#2nd point: We hypothesis that Jacard distance “different” in the group receiving OBL Antibiotics compared to the FQN antibiotics

#We hypothesis that Jacard distance “different” in the group receiving OBL Antibiotics  compared to the FQN
#null hypothesis: there is no difference in the jaccard distance groups receiving obl and fqn 
#alternative there is difference between the two groups

#using welch's t test assuming heteroscedasticity
t.test(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="OBL"],AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="FQN"],var.equal=FALSE )
#p-value = 0.4538 (greater than 0.05 so not significant)
#so we fail to reject the null 

#assessing the assumption
var.test(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="OBL"],AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="FQN"])
# p-value = 0.4596 (not significant) so data is homoscedastic and assumption is not met 

shapiro.test(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="OBL"]) #p-value = 0.05651
shapiro.test(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="FQN"]) #p-value = 0.573

#~~~~~~~~~~~~~~~~~~~
#3rd point:	We hypothesis that Jacard distance is different between the different Antibiotics

#H0: there is no difference in jaccard distance between different antibiotics 
#HA: there is at least one antibiotic different than one other antibiotic


#assuming normality and homoscedasticity
ANOVAModel = aov(D1.D6.Jaccard.distance ~ Antibiotic.class,data=AAD)
summary(ANOVAModel)
coef(ANOVAModel)
#p value=0.287 not signifcant (so we fail to reject the null)

TukeyHSD(ANOVAModel) #there is no significant differenences between the groups
plot(TukeyHSD(ANOVAModel))

#assessing normality and homoscedasticity
shapiro.test(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="OBL"]) #p-value = 0.05651 not significant assume normality
qqnorm(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="OBL"])
qqline(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="OBL"])
hist(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="OBL"])
# we assume normality

shapiro.test(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="PBL"]) #p-value = 0.001216 we can't assume normality
qqnorm(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="PBL"])
qqline(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="PBL"])
hist(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="PBL"])

shapiro.test(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="FQN"]) #p-value = 0.573 assume normality 
qqnorm(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="FQN"])
qqline(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="FQN"])
hist(AAD$D1.D6.Jaccard.distance[AAD$Antibiotic.class=="FQN"])
# since one group is not normal we can use a non-parametric test to be more sure

leveneTest(D1.D6.Jaccard.distance~Antibiotic.class,data=AAD) #p value=0.02891 #heteroscedastic

kruskal.test(D1.D6.Jaccard.distance ~ Antibiotic.class,data=AAD) #p-value = 0.3324 


##################################################################################################################################
#7.	Linear Model
#~~~~~~~~~~~~~~~~
#7.1. simple regression: we used two models first: D1.Shannon.diversity & D6.Shannon.diversity, second: D1.Chao1.diversity & D6.Chao1.diversity and say which model is better.

# STEP 1: Draw a graph of the data to make sure the relationship between D1 & D6 Shannon.diversity make sense
plot(AAD$D1.Shannon.diversity, AAD$D6.Shannon.diversity)

# STEP 2: Do the regression
simple.regression1 <- lm(AAD$D6.Shannon.diversity ~ AAD$D1.Shannon.diversity , data=AAD)

# STEP 3: Look at the R^2, F-value and p-value
summary(simple.regression1)
#In R^2 only 4.875% of the variance in D1.Shannon.diversity is explained by D6.Shannon.diversity. 
#associated p-value is 4.569e-05 (very small). This indicates that the regression model as a whole is statistically significant, suggesting that D6.Shannon.diversity is a significant predictor of D1.Shannon.diversity.

plot(AAD$D1.Shannon.diversity, AAD$D6.Shannon.diversity)
abline(simple.regression1, lwd=5, col="red")
#the results suggest that there is a statistically significant relationship between D6.Shannon.diversity and D1.Shannon.diversity.
#However, the low R-squared value suggests that D6.Shannon.diversity explains only a small portion of the variance in D1.Shannon.diversity.
#Based on these considerations, it appears that the model may not be particularly strong in explaining the variance in the dependent variable. 

#Interpret the regression coefficient
intercept <- coef(simple.regression1)[1]
coefficient <- coef(simple.regression1)[2]

cat("Intercept:", intercept, "\n")
cat("Coefficient:", coefficient, "\n")

qqnorm(simple.regression1 $residuals)
qqline(simple.regression1 $residuals)

#~~~~~~~~~~~~~~ Apply Step 1 2 3 on the Chao Columns and see which model is better 

simple.regression2 <- lm( AAD$D6.Chao1.diversity, data=AAD ~ AAD$D1.Chao1.diversity)
summary(simple.regression2)
#In R^2, only 9.157% of the variance in D1.Chao1.diversity is explained by D6.Chao1.diversity. 
#the associated p-value is 1.597e-08 (very small). This indicates that the regression model as a whole is statistically significant, suggesting that D6.Chao1.diversity is a significant predictor of D1.Chao1.diversity.

plot(AAD$D1.Chao1.diversity, AAD$D6.Chao1.diversity)
abline(simple.regression2, lwd=5, col="red")
#the results suggest that there is a statistically significant relationship between D6.Chao1.diversity and D1.Chao1.diversity. 
#However, the low R-squared value suggests that D6.Chao1.diversity explains only a small portion of the variance in D1.Chao1.diversity.

#Interpret the regression coefficient
intercept <- coef(simple.regression2)[1]
coefficient <- coef(simple.regression2)[2]
cat("Intercept:", intercept, "\n")
cat("Coefficient:", coefficient, "\n")

qqnorm(simple.regression2 $residuals)
qqline(simple.regression2 $residuals)

# In summary, the Chao1 model seems to have a slightly better fit than the Shannon model, as indicated by the higher R-squared (9.157% ). However, it is important to note that both models have relatively low R-squared values, 
#indicating that the chosen predictors (D6.Chao1.diversity and D6.Shannon.diversity) explain only a small portion of the variance in the respective dependent variables.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#7.2. MULTIPLE regression: we used D6.Shannon.diversity with D1 & Antibiotic.class

plot(data.frame(AAD$D6.Shannon.diversity, AAD$D1.Shannon.diversity, AAD$Antibiotic.class))
multiple.regression <- lm(D6.Shannon.diversity ~ D1.Shannon.diversity + Antibiotic.class, data=AAD)
summary(multiple.regression)

if(!require(car)){
  install.packages("car")
}
# Load car package
library(car)

# Produce added variable plots
avPlots(multiple.regression)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BONUS
#1st
confint(simple.regression1, level=0.95)
confint(simple.regression2, level=0.95)

#2nd
Estimate <- predict(lm(AAD$D1.D6.Jaccard.distance ~ AAD$Antibiotic.class), AAD)
Change <- Estimate[2] - Estimate[1]

#Since the value is negative, it indicates a decrease in the average Jacard distance when the Antibiotics value is changed from the first level to the second level.


