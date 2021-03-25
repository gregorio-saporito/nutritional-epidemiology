# load libraries ----------------------------------------------------------
library(tidyverse)
library(bnlearn)
library(Rgraphviz)
library(gRain)

# load data ---------------------------------------------------------------
data <- read_csv("data/breast_cancer_dataset_dse.csv")

# select features
features <- c('NUT2','NUT3','NUT22','NUT23','NUT37','NUT27','NUT9',
              'NUT11','NUT13','NUT16','NUT15','NUT17','NUT31','NUT14','NUT32','NUT18',
              'NUT19','NUT21','NUT33','NUT34','NUT35','NUT28','NUT29','NUT46','NUT36',
              'NUT30','NUT62')

# their respective labels
new_labels <- c('animal_protein', 'vegetable_protein', 'cholesterol', 'saturated_fatty_acids', 'monosaturated_fatty_acids',
                'polyunsaturated_fatty_acids', 'soluble_carbohydrates', 'starch', 'alcohol', 'sodium', 'calcium', 'potassium',
                'phosphorus','iron', 'zinc', 'vitamin_B1', 'vitamin_B2', 'vitamin_C', 'vitamin_B6', 'total_folate',
                'niacin', 'retinol', 'beta_carotene_equivalents', 'lycopene', 'vitamin_D', 'vitamin_E', 'total_fiber')

independent <- select(data,features)
colnames(independent) <- new_labels
dependent <- ifelse(data$V2==2,0,data$V2)

bndf = discretize(independent,method='quantile',breaks=4) %>%
  mutate(cancer=factor(dependent))

# assign quartiles labels with 1 the lowest and 4 the highest
for(i in new_labels){
  levels(bndf[,i]) = c(1:4)
}

# Bayesian Networks -------------------------------------------------------
# learn with hill climbing and AIC criterion algorithm
dag <- hc(bndf,score='aic')
# plot it
graphviz.plot(dag)
# fit it
model <- bn.fit(dag,bndf)

# query the marginal: probability of having cancer
querygrain(as.grain(model), nodes=c("cancer"), type="marginal")

# high starch and low retinol check dimension 3
set.seed(1)
cpquery(model, event = (cancer == "1"),
        evidence =  ((starch == "4") & (retinol == "1")),n=10000000)

# high sodium low beta-carotene check dimension 4
set.seed(1)
cpquery(model, event = (cancer == "1"),
        evidence =  ((sodium == "4") & (beta_carotene_equivalents == "1")),n=10000000)












