library(edgeR)
library(glmnet)
library(pROC)
library(pheatmap)
library("survminer")
library("ggplot2")
require("survival")
library(RColorBrewer)
library(tidyverse)
library(lmtest)

clinInfo = read.csv("Gide_ClinInfo.csv", header = T)
names(clinInfo)[11] = "patient"

data_S100B = read.csv("pseudo_S100B_mat_Gide_decon.csv", header = TRUE, row.names = 1)
data_CD68 = read.csv("pseudo_CD68_mat_Gide_decon.csv", header = TRUE, row.names = 1)
data_CD45 = read.csv("pseudo_CD45_mat_Gide_decon.csv", header = TRUE, row.names = 1)
data_Stroma = data_CD68 + data_CD45

data_S100B[is.na(data_S100B)] = 0
data_CD68[is.na(data_CD68)] = 0
data_CD45[is.na(data_CD45)] = 0
data_Stroma[is.na(data_Stroma)] = 0

selRows = which(rownames(data_S100B) %in% clinInfo$patient)
data_S100B = data_S100B[selRows,]
selRows = which(rownames(data_CD68) %in% clinInfo$patient)
data_CD68 = data_CD68[selRows,]
selRows = which(rownames(data_CD45) %in% clinInfo$patient)
data_CD45 = data_CD45[selRows,]
selRows = which(rownames(data_Stroma) %in% clinInfo$patient)
data_Stroma = data_Stroma[selRows,]

patIds = rownames(data_S100B)
nPat = length(patIds)
resps = rep(0,nPat)
for (i in 1:nPat){
  idx = which(clinInfo$patient == patIds[i])
  idx
  if (clinInfo$response[idx]=="R"){
    resps[i] = 1
  }
}

nSig = 4
sigScores = matrix(0,nPat,nSig)

for (cSig in 1:nSig){
  # if (cSig==1){
  #   coeffs = read.csv("final_coeffs_pseudoBulk_sm.csv")
  #   c_data = data
  # }
  if (cSig==1){
    coeffs = read.csv("final_coeffs_Stroma_df_6_genes.csv")
    c_data = data_Stroma
  }
  #
  if (cSig==2){
    coeffs = read.csv("final_coeffs_S100B_8_genes.csv")
    c_data = data_S100B
  }
  if (cSig==3){
    coeffs = read.csv("final_coeffs_CD45_4_genes.csv")
    c_data = data_CD45
  }
  if (cSig==4){
    coeffs = read.csv("final_coeffs_CD68_8_genes.csv")
    c_data = data_CD68
  }
  
  genes = coeffs$X[2:dim(coeffs)[1]]
  nGene = length(genes)
  cfs = coeffs$s0[2:dim(coeffs)[1]]
  
  for (cPat in 1:nPat){
    vals = rep(0,nGene)
    for (i in 1:nGene){
      idx = which(colnames(c_data)==genes[i])
      if (length(idx)>0){
        if (is.finite(c_data[cPat,idx])){
          vals[i] = c_data[cPat,idx]
        }
      }
    }
    sigScores[cPat,cSig] = sum((as.matrix(vals))*(as.matrix(cfs)))
  } 
  print(sum(vals!=0))  
}

sigScores2 = scale(sigScores)
sigScores2[is.na(sigScores2)] = 0

# # quantile norm
# for (i in 1:nSig){
#   idx = order(sigScores2[,i])
#   for (j in 1:nPat){
#     sigScores2[idx[j],i] = j
#   }
#   sigScores2[,i] = sigScores2[,i] / nPat
# }

sigScores2 = as.data.frame(sigScores2)
sigScores2$resps = t(t(resps))
colnames(sigScores2) = c('pseudoStroma', 'S100B','CD45','CD68','resps')

model0 = glm(resps ~ 1, data = sigScores2, family = "binomial")
# model1 = glm(resps ~ pseudoBulk, data = sigScores2, family = "binomial")
model1 = glm(resps ~ pseudoStroma, data = sigScores2, family = "binomial")
model2 = glm(resps ~ S100B, data = sigScores2, family = "binomial")
model3 = glm(resps ~ CD45, data = sigScores2, family = "binomial")
model4 = glm(resps ~ CD68, data = sigScores2, family = "binomial")

lrtest(model1, model0)
lrtest(model2, model0)
lrtest(model3, model0)
lrtest(model4, model0)
# lrtest(model5, model0)

# model12 = glm(resps ~ pseudoBulk + S100B, data = sigScores2, family = "binomial")
# lrtest(model12, model1)
# 
# model13 = glm(resps ~ pseudoBulk + CD45, data = sigScores2, family = "binomial")
# lrtest(model13, model1)
# 
# model14 = glm(resps ~ pseudoBulk + CD68, data = sigScores2, family = "binomial")
# lrtest(model14, model1)

# model35 = glm(resps ~ S100B + CD68, data = sigScores2, family = "binomial")
# lrtest(model35, model0)
# 
# model91011 = glm(resps ~ S100B_c + CD68_c, data = sigScores2, family = "binomial")
# lrtest(model91011, model0)

### AUC tests

for (cSig in 1:nSig){
  
  # if (cSig==1){
  #   lasso_scores = sigScores2$pseudoBulk
  # }
  if (cSig==1){
    lasso_scores = sigScores2$pseudoStroma
  }
  if (cSig==2){
    lasso_scores = sigScores2$S100B
  }
  if (cSig==3){
    lasso_scores = sigScores2$CD45
  }
  if (cSig==4){
    lasso_scores = sigScores2$CD68
  }
  
  pROC_obj <- roc(sigScores2$resps,lasso_scores,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=FALSE, max.auc.polygon=FALSE, grid=FALSE,
                  print.auc=TRUE, show.thres=TRUE, direction = '<')
  
  pROC_obj
  
}




