library(gRim)
library(tidyverse)
library(caret)
library(ggpubr)
library(dplyr)
library(heatmaply)
theme_set(theme_bw())

#Import dataset
df_original = read_csv("data_logreg.csv")
summary(df_original)

df_diet = read.csv("pesi.csv", encoding = "UTF-8 BOM")
colnames(df_diet)[1] <- "X1"

#union of the two datasets
df = cbind(df_original, df_diet)

df_V2 = read.csv("breast_cancer_dataset_dse.csv", sep = ",")
df_V2 = df_V2["V2"]
summary(df_V2)

df = cbind(df, df_V2)
head(df)

#Cancer = 1/ no cancer = 0
df$V2[df$V2 == 2] = 0
summary(df)

#Scale data
df$X1 = scale(df$X1)
df$X2 = scale(df$X2)
df$X3 = scale(df$X3)
df$X4 = scale(df$X4)

#Check correlation in order to obtain confounding by logistic and linear model:

#Age is not a confounding: it is not associated with cancer but it is associated with diets
summary(glm(V2~V8, data = df, family = binomial))

summary(lm(X1~V8, data = df))
summary(lm(X2~V8, data = df))
summary(lm(X3~V8, data = df))
summary(lm(X4~V8, data = df))

#Edu is a confounding: it is associated with cancer and with diets
summary(glm(V2~V12, data = df, family = binomial))

summary(lm(X1~V12, data = df))
summary(lm(X2~V12, data = df))
summary(lm(X3~V12, data = df))
summary(lm(X4~V12, data = df))

#Menopausal status is a confounding: it is associated both with cancer and all diets
summary(glm(V2~GIN4, data = df, family = binomial))

summary(lm(X1~GIN4, data = df))
summary(lm(X2~GIN4, data = df))
summary(lm(X3~GIN4, data = df))
summary(lm(X4~GIN4, data = df))

#BMI is not a confounding: it is associated with diets but not with Cancer (X2)
summary(glm(V2~BMI, data = df, family = binomial))

summary(lm(X1~BMI, data = df))
summary(lm(X2~BMI, data = df))
summary(lm(X3~BMI, data = df))
summary(lm(X4~BMI, data = df))

#Number of children is a confounding: there is association between the number of children and cancer
summary(glm(V2~V11, data = df, family = binomial))

summary(lm(X1~V11, data = df))
summary(lm(X2~V11, data = df))
summary(lm(X3~V11, data = df))
summary(lm(X4~V11, data = df))

#Smoking is a confounding: there is association between the smoking status and cancer and diet (3 and)
summary(glm(V2~FUM1, data = df, family = binomial))

summary(lm(X1~FUM1, data = df))
summary(lm(X2~FUM1, data = df))
summary(lm(X3~FUM1, data = df))
summary(lm(X4~FUM1, data = df))

#Alcohol is a confounding: there is association between the alcohol status and cancer
summary(glm(V2~ALC1, data = df, family = binomial))

summary(lm(X1~ALC1, data = df))
summary(lm(X2~ALC1, data = df))
summary(lm(X3~ALC1, data = df))
summary(lm(X4~ALC1, data = df))

#Physical activity: a sedentary life between 15-19 is not associated with cancer but is associated with diets
summary(glm(V2~FIS4, data = df, family = binomial))

summary(lm(X1~ FIS4, data = df))
summary(lm(X2~FIS4, data = df))
summary(lm(X3~FIS4, data = df))
summary(lm(X4~FIS4, data = df))


#Logistic model with Cancer as dependent variable and the 4 types of diet as independent variables
model_unadj <- glm( V2 ~ X1 + X2 + X3 + X4, 
              data = df, family = binomial)


summary(model_unadj)


#Logistic model with Cancer as dependent variable and the 4 types of diet as independent variables,
#the model includes also confounding that can be found in the literature
model_adj <- glm( V2 ~ X1 + X2 + X3 + X4 + V12 + GIN4 + V11 + FUM1 + ALC1,
              data = df, family = binomial)

summary(model_adj)

#Odds ratio:
exp(coef(model_unadj))
exp(coef(model_adj))
