# load libraries ----------------------------------------------------------
library(tidyverse)
library(readr)
library(factoextra)
set.seed(123)
library(keras)
library(tensorflow)

# set seeds ---------------------------------------------------------------
use_session_with_seed(123)

# TensorFlow session configuration that uses only a single thread. Multiple threads are a 
# potential source of non-reproducible results, see: https://stackoverflow.com/questions/42022950/which-seeds-have-to-be-set-where-to-realize-100-reproducibility-of-training-res
session_conf <- tf$ConfigProto(intra_op_parallelism_threads = 1L, 
                               inter_op_parallelism_threads = 1L)

# Set TF random seed
tf$set_random_seed(123)

# Create the session using the custom configuration
sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)

# Instruct Keras to use this session
K <- backend()
K$set_session(sess)

# load the data -----------------------------------------------------------
data <- read_csv("data/breast_cancer_dataset_dse.csv")

# select features
features <- c('NUT2','NUT3','NUT22','NUT23','NUT37','NUT27','NUT9',
              'NUT11','NUT13','NUT16','NUT15','NUT17','NUT31','NUT14','NUT32','NUT18',
              'NUT19','NUT21','NUT33','NUT34','NUT35','NUT28','NUT29','NUT46','NUT36',
              'NUT30','NUT62')

# their respective labels
new_labels = c('animal protein', 'vegetable protein', 'cholesterol', 'saturated fatty acids', 'monosaturated fatty acids',
               'polyunsaturated fatty acids', 'soluble carbohydrates', 'starch', 'alcohol', 'sodium', 'calcium', 'potassium',
               'phosphorus','iron', 'zinc', 'vitamin B1', 'vitamin B2', 'vitamin C', 'vitamin B6', 'total folate',
               'niacin', 'retinol', 'beta-carotene equivalents', 'lycopene', 'vitamin D', 'vitamin E', 'total fiber')

# analysis ----------------------------------------------------------------
# scaler
minmax <- function(x) (x - min(x))/(max(x) - min(x))
# standardise the data
x_train <- apply(select(data,features), 2, minmax)

# PCA ---------------------------------------------------------------------
pca <- prcomp(x_train)

summary(pca)
screeplot(pca)
ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, col = factor(data$V2))) + geom_point()

# autoencoder -------------------------------------------------------------
x_train <- as.matrix(x_train)

### some parameters
epochs = 50
verbose = 1

# set model
model <- keras_model_sequential()
model %>%
  layer_dense(units = 4, activation = "tanh", input_shape = ncol(x_train), name = "bottleneck", kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 123)) %>%
  layer_dense(units = ncol(x_train),kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 123))

# view model layers
summary(model)

# compile model
model %>% compile(loss = "mean_squared_error", optimizer = "adam")

# fit model
model %>% fit(x = x_train, y = x_train, 
              epochs = epochs, verbose = verbose)

# evaluate the performance of the model
mse.ae2 <- evaluate(model, x_train, x_train)
mse.ae2

# extract the bottleneck layer
intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, x_train)

# get weights
weights = get_weights(model)[[1]]

# plot the table of weights
weights %>% as.data.frame() %>% 
  mutate(Var2=new_labels) %>%
  gather(key='Var1',value='value',-Var2) %>%
  mutate(Var2=factor(Var2,levels=rev(new_labels))) %>%
  ggplot(aes(x=Var1, y=Var2, fill=value, label=round(value,2))) +
  geom_tile() +
  geom_text(size=2) +
  scale_fill_gradient2() +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  ) +
  ggtitle('weights')

# save dietary patterns extracted with autoencoder to csv
write_csv(data.frame(intermediate_output),'outputs/one_hidden_layer.csv')

# exploring the dimensions
ggplot(data.frame(PC1 = intermediate_output[,1], PC2 = intermediate_output[,2]), aes(x = PC1, y = PC2, col = factor(data$V2))) + geom_point()
ggplot(data.frame(PC1 = intermediate_output[,1], PC2 = intermediate_output[,3]), aes(x = PC1, y = PC2, col = factor(data$V2))) + geom_point()
ggplot(data.frame(PC1 = intermediate_output[,1], PC2 = intermediate_output[,4]), aes(x = PC1, y = PC2, col = factor(data$V2))) + geom_point()
ggplot(data.frame(PC1 = intermediate_output[,2], PC2 = intermediate_output[,3]), aes(x = PC1, y = PC2, col = factor(data$V2))) + geom_point()
ggplot(data.frame(PC1 = intermediate_output[,2], PC2 = intermediate_output[,4]), aes(x = PC1, y = PC2, col = factor(data$V2))) + geom_point()
ggplot(data.frame(PC1 = intermediate_output[,3], PC2 = intermediate_output[,4]), aes(x = PC1, y = PC2, col = factor(data$V2))) + geom_point()

# interpret the dimensions with correlations
corcheck = data.frame(matrix(nrow = 27, ncol = 4))
colnames(corcheck) = c('DIM1','DIM2','DIM3','DIM4')
rownames(corcheck) = new_labels

for(i in 1:27){
  for(d in 1:4){
    corcheck[i,d] = cor(x_train[,i],intermediate_output[,d])
  }
}

# exploring the potential presence of nonlinearities
plot(x_train[,27],intermediate_output[,1])
plot(x_train[,20],intermediate_output[,3])
plot(x_train[,15],intermediate_output[,4])
plot(x_train[,10],intermediate_output[,1])
plot(x_train[,5],intermediate_output[,2])
plot(x_train[,1],intermediate_output[,3])

# plot table to compare the correlations
corcheck %>%
  mutate(Var2=new_labels) %>%
  gather(key='Var1',value='value',-Var2) %>%
  mutate(Var2 = factor(Var2, levels=rev(new_labels))) %>% 
  ggplot(aes(x=Var1, y=Var2, fill=value, label=round(value,2))) +
  geom_tile() +
  geom_text(size=2) +
  scale_fill_gradient2() +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  ) +
  ggtitle('correlations')

# compare corcheck and weights to see if they are similar
library(CCA)
X = as.matrix(corcheck)
colnames(X) = c('C1','C2','C3','C4')
Y = as.matrix(weights)
rownames(Y) = new_labels
colnames(Y) = c('W1','W2','W3','W4')
res.cc <- cc(X,Y)
plt.cc(res.cc, var.label = T)

# compare the weights with the loadings of the pca
# they look indeed different from the pca so we are capturing different dietary patterns
Xpca = pca$rotation[,1:4]
colnames(Xpca) = c('PC1','PC2','PC3','PC4')
rownames(Xpca) = new_labels
auto.pca <- cc(Xpca,Y)
plt.cc(auto.pca, var.label = T)

# cor between PCA components
cor(pca$x[,1:4])
# cor between variables in the latent space with autoencoder, larger than PCA!
cor(intermediate_output)


# compare the performances ------------------------------------------------
# how many dimensions to check
setdim=3
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

ae.mse <- rep(NA, setdim)
for(k in 1:setdim){
  modelk <- keras_model_sequential()
  
  modelk %>%
    layer_dense(units = k, activation = "tanh", input_shape = ncol(x_train), name = "bottleneck", kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 123)) %>%
    layer_dense(units = ncol(x_train),kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 123))
  
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

