### Data Manipulation
library(ForImp)

# Importing data
path <- '../data/breast_cancer_dataset_dse.csv'
data <- read.csv(path)

# Selected variables
vars <- c('CENTRO', 'V2', 'V8', 'V9', 'V11', 'V12', 'FUM1', 'FUM7', 'FIS2', 
          'FIS4', 'ANTR1', 'ANTR2', 'ANTR12', 'ALC4', 'ALC8', 'ALC12',
          'ALC16', 'ALC20', 'ALC25', 'ANAM1', 'ANAM2', 'ANAM3', 'ANAM4', 
          'ANAM5', 'ANAM6', 'ANAM7', 'ANAM8', 'ANAM9', 'ANAM10', 'ICDAN11', 
          'ANAM11', 'ANAM12', 'ANAM13', 'ANAM14', 'ANAM15', 'ANAM16',
          'ANAM17', 'ANAM18', 'ANAM19', 'ANAM20', 'ANAM21', 'ANAM22', 'ANAM23', 
          'ANAM24', 'ANAM25', 'ANAM26', 'ANAM27', 'FAM3', 'FAM8', 'FAM13', 
          'FAM18', 'FAM23', 'FAM28', 'GIN1', 'GIN3', 'GIN4', 'GIN7', 'GIN8', 'GIN9')

df <- data[vars]

# For alcohol variables, 998 (occasional) -> 0
alc_var <- c('ALC4', 'ALC8', 'ALC12', 'ALC16', 'ALC20', 'ALC25')
df[alc_var][df[alc_var]==998] <- 0

# Missing value code for each variable of interest
na_99 <- c('V12', 'FUM7', 'ANAM2', 'ANAM3', 'ANAM5', 'ANAM6', 'ANAM7', 'ANAM8', 'ANAM12',
           'ANAM14', 'ANAM18', 'ANAM21', 'ANAM22', 'ANAM24', 'ANAM25', 'ANAM27', 'GIN1')
na_98 <- c('ANAM4', 'ANAM6', 'ANAM15', 'ANAM16')
na_9 <- c('FUM1', 'FIS2', 'FIS4', 'ANTR12', 'GIN3', 'GIN4', 'GIN7', 'GIN8', 'GIN9')
na_999 <- c('ANTR1', 'ANTR2', 'ALC4', 'ALC8', 'ALC12', 'ALC16', 'ALC20')
na_0 <- c('ANTR12', 'GIN3')

# Replace code with NA
df[na_99][df[na_99] == 99] <- NA
df[na_98][df[na_98] == 98] <- NA
df[na_9][df[na_9] == 9] <- NA
df[na_999][df[na_999] == 999] <- NA
df[na_0][df[na_0] == 0] <- NA

# Forward Imputation algorithm
df <- ForImp(df)
df <- as.data.frame(df)
colnames(df) <- vars

# Generate one alcohol var by summing the others and drop them
df$ALC0 <- rowSums(df[alc_var])
df <- df[, !(names(df) %in% alc_var)]

# Transform medical history variables in binary
med_var = c('ANAM1', 'ANAM2', 'ANAM3', 'ANAM4', 'ANAM5', 'ANAM6', 'ANAM7', 'ANAM8',
            'ANAM9', 'ANAM10', 'ICDAN11', 'ANAM11', 'ANAM12', 'ANAM13', 'ANAM14', 'ANAM15',
            'ANAM16', 'ANAM17', 'ANAM18', 'ANAM19', 'ANAM20', 'ANAM21', 'ANAM22', 'ANAM23',
            'ANAM24', 'ANAM25', 'ANAM26', 'ANAM27')
df[med_var][df[med_var] > 0] <- 1

# Convert continuous variables to categorical
to_cat <- c('V8', 'V11', 'V12', 'FUM7', 'ANTR1', 'ANTR2', 'ALC0', 'GIN1', 'GIN7', 'GIN8', 'GIN9')

