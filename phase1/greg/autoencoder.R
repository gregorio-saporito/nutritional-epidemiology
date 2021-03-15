library(keras)
library(tidyverse)
library(readr)
library(factoextra)

data <- read_csv("data/breast_cancer_dataset_dse.csv")

features <- c('NUT12','NUT2','NUT3','NUT5','NUT6','NUT22','NUT23','NUT37','NUT27','NUT9',
              'NUT11','NUT13','NUT16','NUT15','NUT17','NUT31','NUT14','NUT32','NUT18',
              'NUT19','NUT21','NUT33','NUT34','NUT35','NUT28','NUT29','NUT46','NUT36',
              'NUT30','NUT62')

# set training data
minmax <- function(x) (x - min(x))/(max(x) - min(x))
# standardise the data
x_train <- apply(select(data,features), 2, minmax)

# PCA
pca <- prcomp(x_train)

summary(pca)

screeplot(pca)

ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, col = factor(data$V2))) + geom_point()

# autoencoder
x_train <- as.matrix(x_train)

### some parameters
epochs = 50
verbose = 1
# how many dimensions to check later
setdim=5

# set model
model <- keras_model_sequential()
model %>%
  layer_dense(units = 20, activation = "tanh", input_shape = ncol(x_train)) %>%
  layer_dense(units = 4, activation = "tanh", name = "bottleneck") %>%
  layer_dense(units = 20, activation = "tanh") %>%
  layer_dense(units = ncol(x_train))

# view model layers
summary(model)

# compile model
model %>% compile(
  loss = "mean_squared_error", 
  optimizer = "adam"
)

# fit model
model %>% fit(
  x = x_train, 
  y = x_train, 
  epochs = epochs,
  verbose = verbose
)

# evaluate the performance of the model
mse.ae2 <- evaluate(model, x_train, x_train)
mse.ae2

# extract the bottleneck layer
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, x_train)

# print dietary patterns to csv
write_csv(data.frame(intermediate_output),'outputs/autoencoder.csv')

ggplot(data.frame(PC1 = intermediate_output[,1], PC2 = intermediate_output[,2]), aes(x = PC1, y = PC2, col = factor(data$V2))) + geom_point()

# pCA reconstruction
pca.recon <- function(pca, x, k){
  mu <- matrix(rep(pca$center, nrow(pca$x)), nrow = nrow(pca$x), byrow = T)
  recon <- pca$x[,1:k] %*% t(pca$rotation[,1:k]) + mu
  mse <- mean((recon - x)^2)
  return(list(x = recon, mse = mse))
}

xhat <- rep(NA, setdim)
for(k in 1:setdim){
  xhat[k] <- pca.recon(pca, x_train, k)$mse
}

ae.mse <- rep(NA, 5)
for(k in 1:setdim){
  modelk <- keras_model_sequential()
  
  modelk %>%
    layer_dense(units = 20, activation = "tanh", input_shape = ncol(x_train)) %>%
    layer_dense(units = k, activation = "tanh", name = "bottleneck") %>%
    layer_dense(units = 20, activation = "tanh") %>%
    layer_dense(units = ncol(x_train))
  
  modelk %>% compile(
    loss = "mean_squared_error", 
    optimizer = "adam"
  )
  
  modelk %>% fit(
    x = x_train, 
    y = x_train, 
    epochs = epochs,
    verbose = verbose
  )
  
  print(paste0('number of latent variables estimated: ',k))
  
  ae.mse[k] <- unname(evaluate(modelk, x_train, x_train))
}

df <- data.frame(k = c(1:setdim, 1:setdim), mse = c(xhat, ae.mse), method = c(rep("pca", setdim), rep("autoencoder", setdim)))
ggplot(df, aes(x = k, y = mse, col = method)) + geom_line()

