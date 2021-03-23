library(gRim)
library(tidyverse)
library(caret)
library(ggpubr)
library(dplyr)
theme_set(theme_bw())


#Import dataset
df_original = read_csv("breast_cancer_dataset_dse.csv")
colnames(df_original)[1] <- "CENTRO"

summary(df_original)

df_diet = read.csv("pesi.csv", encoding = "UTF-8 BOM")
colnames(df_diet)[1] <- "X1"

#union of the two datasets
df = cbind(df_original, df_diet)

#Cancer = 1/ no cancer = 0
df$V2[df$V2 == 2] = 0
summary(df)


#Logistic model with Cancer as dependent variable and the 4 types of diet as independent variables
model <- glm( V2 ~ X1 + X2 + X3 + X4, 
              data = df, family = binomial)


summary(model)

#Logistic model with Cancer as dependent variable and the 4 types of diet as independent variables,
#the model includes also confounding that can be found in the literature
model <- glm( V2 ~ X1 + X2 + X3 + X4 + V8 + V12 + GIN4 + ANTR1 + ANTR2,
              data = df, family = binomial)

summary(model)

