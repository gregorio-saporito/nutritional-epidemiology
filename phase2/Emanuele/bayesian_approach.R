library(gRbase)
library(Rgraphviz)
library(igraph)
library(bnlearn)
library(dplyr)
library(heatmaply)
theme_set(theme_bw())


#Import dataset:
#df_cleaned is the dataframe with the variables categorized.
df_cleaned = read.csv("data_cleaned.csv", sep = ";", encoding = "UTF-8 BOM")
colnames(df_cleaned)[1] <- "CENTRO"

df_diet = read.csv("pesi.csv", encoding = "UTF-8 BOM")
colnames(df_diet)[1] <- "X1"

#union of the two datasets
df = cbind(df_cleaned, df_diet)

#Cancer = 1/ no cancer = 0
df$V2[df$V2 == 2] = 0
summary(df)

df$X1 = scale(normalize(df$X1))
df$X2 = scale(normalize(df$X2))
df$X3 = scale(normalize(df$X3))
df$X4 = scale(normalize(df$X4))

hist(df$X4)


#categorize the diets
df[,c("X1", "X2", "X3", "X4")] = as.data.frame(sapply(df[,c("X1", "X2", "X3", "X4")], function (x) {
  qx <- quantile(x)
  cut(x, qx, include.lowest = TRUE,
      labels = 1:4)
}))


#Factorize the categorical variables
names <- c(1:ncol(df))
df[,names] <- lapply(df[,names] , factor)
str(df)


#Diet_tot

diet = c("V2", "EDU", "AGE", "GIN4", "ANTR0", "X4", "X2", "X3")
diet.df = df[,diet]

wl = data.frame(from =c("X4", "X2", "X3", "EDU", "AGE", "GIN4", "ANTR0"), to = c("V2", "V2", "V2", "V2", "V2", "V2", "V2"))
bl = data.frame(from = c("V2", "V2"), to = c("EDU", "GIN4"))

#Models:

#PC stable structure learning
diet.model_pc = pc.stable(diet.df, whitelist = wl, blacklist = bl)
plot(diet.model_pc)
pc_mat = amat(diet.model_pc)
pc_mat

#Grow-Shrink structure learning
diet.model_gs = gs(diet.df, whitelist = wl, blacklist = bl)
plot(diet.model_gs)
gs_mat =amat(diet.model_gs)
gs_mat

#Hill-climbing structure learning
diet.model_hc = hc(diet.df, score = "bic", whitelist = wl, blacklist = bl)
plot(diet.model_hc)
hc_mat = amat(diet.model_hc)

#Tabu greedy Search
diet.model_tabu = tabu(diet.df, score = "bic", whitelist = wl, blacklist = bl)
plot(diet.model_tabu)
tabu_mat = amat(diet.model_tabu)

#Max-Min Hill Climbing
diet.model_MM = rsmax2(diet.df, whitelist = wl, blacklist = bl)
plot(diet.model_MM)
MM_mat = amat(diet.model_MM)

mat_sum = hc_mat + gs_mat + pc_mat + tabu_mat + MM_mat
mat_sum

#create a bayesian network based on the results of the structure learning algorithms:

adj = mat_sum
adj[which(mat_sum < 2)] = 0
adj[which(mat_sum >= 2)] = 1

adj

model = empty.graph(colnames(mat_sum))
title(main="Bayesian Network ", col.main="red")

adj["X4","X3"] = 0L
adj["X2", "X4"] = 0L
adj["X3", "X4"] = 0L
adj["X3", "X2"] = 0L

#creation of the bayesian network based on the adjacency matrix 
#retrieved by the combination of the learning structure models applied

amat(model) = adj
model = pdag2dag(model, ordering =colnames(mat_sum))

graphviz.plot(model)

fit = bn.fit(model, diet.df)


#Calculate the probability of cancer by varying the confounding variables for each diet:

age_values = unique(df$AGE)

#X2_____________________________________________________________________________
#Probability of cancer by varying X2
probs_X2 = c()

for (i in 1:4){
  value = toString(i)
  prob = cpquery(fit, (V2 == "1") , (X2 == value))
  probs_X2[i] = prob
  }

probs_X2

#X2: Probability of cancer by varying X2 with AGE as confounder:

list_probs_X2 = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & AGE == age_values[j]))
    probs[i] = prob
  }
  
  list_probs_X2[[j]] = probs
}
list_probs_X2

#X2: Probability with confounding X4:

list_probs_X2_conf = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & X4 == toString(j)))
    probs[i] = prob
  }
  
  list_probs_X2_conf[[j]] = probs
}
list_probs_X2_conf

##X3_____________________________________________________________________________
#prob of X3:
probs_X3 = c()
for (i in 1:4){
  value = toString(i)
  prob = cpquery(fit, (V2 == "1") , (X3 == value))
  probs_X3[i] = prob
}

probs_X3

#prob of X3 with AGE as confounding:
list_probs_X3 = list()
for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X3 == value & AGE == age_values[j]))
    probs[i] = prob
  }
  
  list_probs_X3[[j]] = probs
}
list_probs_X3

#prob of X3 with X2 as confounding:

list_probs_X3_conf = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X3 == value & X2 == toString(j)))
    probs[i] = prob
  }
  
  list_probs_X3_conf[[j]] = probs
}
list_probs_X3_conf

##X4_____________________________________________________________________________
#prob of X4

probs_X4 = c()

for (i in 1:4){
  value = toString(i)
  prob = cpquery(fit, (V2 == "1") , (X4 == value))
  probs_X4[i] = prob
}
probs_X4

#prob of X4 with AGE as confounding

list_probs_X4 = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X4 == value & AGE == age_values[j]))
    probs[i] = prob
  }
  
  list_probs_X4[[j]] = probs
}

list_probs_X4

prob_cancer = cpquery(fit, (V2 == "1"), TRUE)
prob_cancer

#Plots
#X2
g_range = range(0.3, (probs_X2+0.2),(list_probs_X2_conf[[1]]+0.2), (list_probs_X2_conf[[2]]+0.2))

plot(probs_X2, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
abline(h = prob_cancer, col = "grey")
lines(list_probs_X2_conf[[1]], type = "o", col = "green",lwd = 1)
lines(list_probs_X2_conf[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(list_probs_X2_conf[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
lines(list_probs_X2_conf[[4]], type = "o", col = "blue", pch = 19, lwd = 1)
title(main="Probabilities of cancer increasing value of diet X2", col.main="red")
title(xlab="Quartiles of diet X2")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("1","2", "3", "4"), cex=0.8, 
       col=c("green","red", "orange", "blue"), pch = 21, lty = 1, title = "Quartiles of diet X4;", horiz = FALSE);


#X3
g_range = range(0.2, (probs_X3+0.2),(list_probs_X3[[1]]+0.2), (list_probs_X3[[2]]+0.2))

plot(probs_X3, type = "o", col = "black",lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
abline(h = prob_cancer, col = "grey")
lines(list_probs_X3_conf[[1]], type = "o", col = "green")
lines(list_probs_X3_conf[[2]], type = "o", col = "red")
lines(list_probs_X3_conf[[3]], type = "o", col = "orange", pch = 19)
lines(list_probs_X3_conf[[4]], type = "o", col = "blue", pch = 19)
title(main="Probabilities of cancer increasing value of diet X3", col.main="red")
title(xlab="Quartiles of diet X3")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("1","2", "3", "4"), cex=0.8, 
       col=c("green","red", "orange", "blue"), pch = 21, lty = 1, title = "Quartiles of diet X2;", horiz = FALSE);

#X4
g_range = range(0.1, (probs_X4+0.2),(list_probs_X4[[1]]+0.2), (list_probs_X4[[2]]+0.2))

plot(probs_X4, type = "o", col = "black", ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
abline(h = prob_cancer, col = "grey")
lines(list_probs_X4[[1]], type = "o", col = "red", pch = 19)
lines(list_probs_X4[[2]], type = "o", col = "green", pch = 19)
lines(list_probs_X4[[3]], type = "o", col = "purple")
title(main="Probabilities of cancer increasing value of diet X4", col.main="red")
title(xlab="Quartiles of diet X4")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("<=40",">40 & < 60", ">=60", "all"), cex=0.8, 
       col=c("purple","red", "green", "black"), pch = 21, lty = 1, title = "Age of the patient:")

