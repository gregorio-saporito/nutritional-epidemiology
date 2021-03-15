library(tidyverse)

data <- read_csv("data/breast_cancer_dataset_dse.csv")

selected_nutrients <- c('NUT12','NUT2','NUT3','NUT5','NUT6','NUT22','NUT23','NUT37','NUT27','NUT9',
              'NUT11','NUT13','NUT16','NUT15','NUT17','NUT31','NUT14','NUT32','NUT18',
              'NUT19','NUT21','NUT33','NUT34','NUT35','NUT28','NUT29','NUT46','NUT36',
              'NUT30','NUT62')

returnvars <- function(pattern,number){
  map(1:number,function(i){paste0(pattern,i)}) %>% unlist()
}


