library(tidyverse)

data <- read_csv("data/breast_cancer_dataset_dse.csv")

# select features
features <- c('NUT2','NUT3','NUT22','NUT23','NUT37','NUT27','NUT9',
              'NUT11','NUT13','NUT16','NUT15','NUT17','NUT31','NUT14','NUT32','NUT18',
              'NUT19','NUT21','NUT33','NUT34','NUT35','NUT28','NUT29','NUT46','NUT36',
              'NUT30','NUT62')

# their respective labels
new_labels <- c('animal protein', 'vegetable protein', 'cholesterol', 'saturated fatty acids', 'monosaturated fatty acids',
               'polyunsaturated fatty acids', 'soluble carbohydrates', 'starch', 'alcohol', 'sodium', 'calcium', 'potassium',
               'phosphorus','iron', 'zinc', 'vitamin B1', 'vitamin B2', 'vitamin C', 'vitamin B6', 'total folate',
               'niacin', 'retinol', 'beta-carotene equivalents', 'lycopene', 'vitamin D', 'vitamin E', 'total fiber')

independent <- select(data,features)
colnames(independent) <- new_labels
dependent <- ifelse(data$V2==2,0,data$V2)

cleaned = data.frame(cancer=factor(dependent)) %>% bind_cols(independent)

mylogit <- glm(cancer ~ ., data = cleaned, family = "binomial")

summarylogit = summary(mylogit)

summarylogit$coefficients %>% as.data.frame() %>% 
  filter(`Pr(>|z|)`< 0.1)

# odds ratios and 95% CI
odds_ratios = exp(cbind(OR = coef(mylogit), confint(mylogit)))
odds_ratios %>% as.data.frame() %>%
  filter(`2.5 %` > 1 | `97.5 %`<1)

# ROC curve
prob=predict(mylogit,type=c("response"))
cleaned$prob=prob
library(pROC)
g <- roc(cancer ~ prob, data = cleaned)
plot(g)    
