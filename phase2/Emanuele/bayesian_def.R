library(gRbase)
library(Rgraphviz)
library(igraph)
library(bnlearn)
library(dplyr)
library(heatmaply)

#Import dataset:
#df_cleaned is the dataframe with the variables categorized.
df_cleaned = read.csv("data_cleaned.csv", sep = ";", encoding = "UTF-8 BOM")
colnames(df_cleaned)[1] <- "CENTRO"

df_diet = read.csv("pesi.csv", encoding = "UTF-8 BOM")
colnames(df_diet)[1] <- "X1"

#union of the two datasets
df = cbind(df_cleaned, df_diet)
head(df)

#Cancer = 1/ no cancer = 0
df$V2[df$V2 == 2] = 0
summary(df)

#Scale data
df$X1 = scale(df$X1)
df$X2 = scale(df$X2)
df$X3 = scale(df$X3)
df$X4 = scale(df$X4)


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

colnames(df)

#Diet_tot
diet = c("AGE", "EDU", "GIN4", "ANTR0", "CHILD", "SMOKE", "ALCOHOL", "FIS4", "X2","X3","X4", "V2")
diet.df = df[,diet]

wl = data.frame(from =c("AGE", "AGE", "AGE","EDU","EDU","EDU","GIN4","GIN4","GIN4","GIN4", "ANTR0", "ANTR0", "CHILD","CHILD", "SMOKE",
                        "SMOKE","SMOKE","SMOKE","ALCOHOL","ALCOHOL","ALCOHOL","ALCOHOL","FIS4","FIS4","FIS4", "X2", "X3", "X4"),
                
                to = c("X2", "X3", "X4","V2","X2","X3", "V2", "X2", "X3", "X4","X2","X4", "V2","X2", "V2", "X2", "X3", "X4", 
                       "V2", "X2", "X3", "X4", "X2", "X3", "X4", "V2", "V2", "V2"))
bl = data.frame(from = c("V2"), to = c("EDU"))


#PC stable structure learning
diet.model_pc = pc.stable(diet.df, blacklist = bl, whitelist = wl)
graphviz.plot(diet.model_pc)
pc_mat = amat(diet.model_pc)

#Grow-Shrink structure learning
diet.model_gs = gs(diet.df, whitelist = wl, blacklist = bl)
graphviz.plot(diet.model_gs)
gs_mat =amat(diet.model_gs)

#Hill-climbing structure learning
diet.model_hc = hc(diet.df, score = "bic", whitelist = wl, blacklist = bl)
graphviz.plot(diet.model_hc)
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
adj[which(mat_sum < 3)] = 0
adj[which(mat_sum >= 3)] = 1

adj

model = empty.graph(colnames(adj))

#creation of the bayesian network based on the adjacency matrix 
#retrieved by the combination of the learning structure models applied
amat(model) = adj

gR = graphviz.plot(model, layout = "dot", shape = "rectangle", 
              highlight = list(nodes = c("V2"), col = c("tomato"), fill = "orange"))

node.attrs = nodeRenderInfo(gR)
node.attrs$textCol[c("X3","X2","X4")] ="tomato"
node.attrs$shape[c("X3","X2","X4")] ="ellipse"
nodeRenderInfo(gR) = node.attrs
renderGraph((gR))

arc.attrs = edgeRenderInfo(gR)
arc.attrs$col[c("V2~ANTR0", "ANTR0~X2", "X2~V2")] = "red"
arc.attrs$col[c("AGE~X4", "X4~V2")] = "orange"
arc.attrs$col[c("GIN4~V2", "GIN4~X3","X3~V2")] = "green"
edgeRenderInfo(gR) = arc.attrs
renderGraph(gR)

fit = bn.fit(model, diet.df)

#Values of the variables:
age_values = unique(df$AGE)
edu_values = unique(df$EDU)
smoke_values = unique(df$SMOKE)
alcohol_values = unique(df$ALCOHOL)
bmi_values = unique(df$ANTR0)
child_values = unique(df$CHILD)

#X2_____________________________________________________________________________
#Probability of cancer by varying X2
probs_X2 = c()

for (i in 1:4){
  value = toString(i)
  prob = cpquery(fit, (V2 == "1") , (X2 == value))
  probs_X2[i] = prob
}
probs_X2

#Probability of cancer by varying X2 with AGE as confounder:

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


#X2: Probability with confounding EDU:
str(df$EDU)

X2_EDU = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & EDU == edu_values[j]))
    probs[i] = prob
  }
  
  X2_EDU[[j]] = probs
}
X2_EDU

#X2: Probability with confounding GIN4:
str(df$GIN4)

X2_GIN = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & GIN4 == toString(j)))
    probs[i] = prob
  }
  
  X2_GIN[[j]] = probs
}
X2_GIN

#X2: Probability with confounding SMOKE:
str(df$SMOKE)

X2_SMOKE = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & SMOKE == smoke_values[j]))
    probs[i] = prob
  }
  
  X2_SMOKE[[j]] = probs
}
X2_SMOKE

#X2: Probability with confounding ALCOHOL:
str(df$ALCOHOL)

X2_ALCOHOL = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & ALCOHOL == alcohol_values[j]))
    probs[i] = prob
  }
  
  X2_ALCOHOL[[j]] = probs
}
X2_ALCOHOL

#X2: Probability with confounding FIS4:
str(df$FIS4)

X2_FIS4 = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & FIS4 == toString(j)))
    probs[i] = prob
  }
  
  X2_FIS4[[j]] = probs
}
X2_FIS4


#X2: Probability with confounding BMI:
str(df$ANTR0)

X2_BMI = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & ANTR0 == bmi_values[i]))
    probs[i] = prob
  }
  
  X2_BMI[[j]] = probs
}
X2_BMI

#X2: Probability with confounding CHILD:
str(df$CHILD)

X2_CHILD = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & CHILD == child_values[j]))  
    probs[i] = prob
  }
  
  X2_CHILD[[j]] = probs
}
X2_CHILD

#X2 PLOT:
g_range = range(0.25, 0.5)

plot(probs_X2, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(X2_BMI[[1]], type = "o", col = "green",lwd = 1)
lines(X2_BMI[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(X2_BMI[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
lines(X2_BMI[[4]], type = "o", col = "blue", pch = 19, lwd = 1)
title(main="Probabilities of cancer given X2 and BMI", col.main="red")
title(xlab="Quartiles of diet X2")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("1","2", "3", "4", "no strat."), cex=0.8, 
       col=c("green","red", "orange", "blue", "black"), pch = 21, lty = 1, title = "Stratified by BMI:", horiz = FALSE);


#X3_____________________________________________________________________________
#Probability of cancer by varying X3
probs_X3 = c()

for (i in 1:4){
  value = toString(i)
  prob = cpquery(fit, (V2 == "1") , (X3 == value))
  probs_X3[i] = prob
}
probs_X3

#Probability of cancer by varying X3 with AGE as confounding:

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

#X3: Probability with confounding EDU:
str(df$EDU)

X3_EDU = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X3 == value & EDU == edu_values[j]))
    probs[i] = prob
  }
  
  X3_EDU[[j]] = probs
}
X3_EDU

#X3: Probability with confounding GIN4:
str(df$GIN4)

X3_GIN = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X3 == value & GIN4 == toString(j)))
    probs[i] = prob
  }
  
  X3_GIN[[j]] = probs
}
X3_GIN

#X3: Probability with confounding SMOKE:
str(df$SMOKE)

X3_SMOKE = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X3 == value & SMOKE == smoke_values[j]))
    probs[i] = prob
  }
  
  X3_SMOKE[[j]] = probs
}
X3_SMOKE

#X3: Probability with confounding ALCOHOL:
str(df$ALCOHOL)

X3_ALCOHOL = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X3 == value & ALCOHOL == alcohol_values[j]))
    probs[i] = prob
  }
  
  X3_ALCOHOL[[j]] = probs
}
X3_ALCOHOL

alcohol_values

#X3: Probability with confounding FIS4:
str(df$FIS4)

X3_FIS4 = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X3 == value & FIS4 == toString(j)))
    probs[i] = prob
  }
  
  X3_FIS4[[j]] = probs
}
X3_FIS4

#X3 PLOT:
g_range = range(0, 0.6)

#X3_EDU
plot(probs_X3, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(X3_EDU[[1]], type = "o", col = "green",lwd = 1)
lines(X3_EDU[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(X3_EDU[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
title(main="Probabilities of cancer given X3 and EDU", col.main="red")
title(xlab="Quartiles of diet X3")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("elementary","middle/high", "university", "P(cancer|X3)"), cex=0.8, 
       col=c("green","red", "orange", "black"), pch = 21, lty = 1, title = "Level of education;", horiz = FALSE);

#X3_GIN
g_range = range(0.025, 0.5)
plot(probs_X3, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(X3_GIN[[1]], type = "o", col = "green",lwd = 1)
lines(X3_GIN[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(X3_GIN[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
title(main="Probabilities of cancer given X3 and Menopausal Status", col.main="red")
title(xlab="Quartiles of diet X3")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("pre-menopause","peri-menopause", "post-menopause", "No strat."), cex=0.8, 
       col=c("green","red", "orange", "black"), pch = 21, lty = 1, title = "Stratified by Menopausal Status:", horiz = FALSE);


#X4_____________________________________________________________________________
#Probability of cancer by varying X4
probs_X4 = c()

for (i in 1:4){
  value = toString(i)
  prob = cpquery(fit, (V2 == "1") , (X4 == value))
  probs_X4[i] = prob
}
probs_X4

#Probability of cancer by varying X4 with AGE as confounding:

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


#X4: Probability with confounding GIN4:
str(df$GIN4)

X4_GIN = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X4 == value & GIN4 == toString(j)))
    probs[i] = prob
  }
  
  X4_GIN[[j]] = probs
}
X4_GIN

#X4: Probability with confounding BMI:
str(df$ANTR0)

X4_BMI = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X4 == value & ANTR0 == bmi_values[i]))
    probs[i] = prob
  }
  
  X4_BMI[[j]] = probs
}
X4_BMI

#X4: Probability with confounding SMOKE:
str(df$SMOKE)

X4_SMOKE = list()

for (j in 1:3){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X4 == value & SMOKE == smoke_values[j]))
    probs[i] = prob
  }
  
  X4_SMOKE[[j]] = probs
}
X4_SMOKE

#X4: Probability with confounding ALCOHOL:
str(df$ALCOHOL)

X4_ALCOHOL = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X4 == value & ALCOHOL == alcohol_values[j]))
    probs[i] = prob
  }
  
  X4_ALCOHOL[[j]] = probs
}
X4_ALCOHOL

#X4: Probability with confounding FIS4:
str(df$FIS4)

X4_FIS4 = list()

for (j in 1:4){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X4 == value & FIS4 == toString(j)))
    probs[i] = prob
  }
  
  X4_FIS4[[j]] = probs
}
X4_FIS4

#X4 Plots
#X4_AGE
g_range = range(0.1, 0.5)
plot(probs_X4, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(list_probs_X4[[1]], type = "o", col = "green",lwd = 1)
lines(list_probs_X4[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(list_probs_X4[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
title(main="Probabilities of cancer given X4 and AGE", col.main="red")
title(xlab="Quartiles of diet X4")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("<40",">=40 & < 60", ">=60", "No strat."), cex=0.8, 
       col=c("green","red", "orange", "black"), pch = 21, lty = 1, title = "Stratified by Age", horiz = FALSE)
