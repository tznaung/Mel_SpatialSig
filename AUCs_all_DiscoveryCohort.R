library("survminer")
library("ggplot2")
require("survival")
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(edgeR)
library(pROC)


pBulk = read.csv("PseudoBulk_df.csv", row.names = 1) # pseudoBulk_av.csv
pStroma = read.csv("PseudoStroma_df.csv", row.names = 1) # pseudoStroma_av.csv
S100B = read.csv("S100B_df.csv", row.names = 1) # S100B_av.csv
CD68 = read.csv("CD68_df.csv", row.names = 1) # CD68_av.csv
CD45 = read.csv("CD45_df.csv", row.names = 1) # CD45_av.csv

data_pBulk = pBulk[, 1:7320]
data_pStroma = pStroma[, 2:2336]
data_S100B = S100B[, 3:9286]
data_CD45 = CD45[, 3:3738]
data_CD68 = CD68[, 3:1713]


patIds = rownames(data_pBulk)
nPat = length(patIds)
resps = rep(0,nPat)
for (i in 1:nPat){
  idx = which(pBulk$SpotID == patIds[i])
  idx
  if (pBulk$response[idx]=="yes"){
    resps[i] = 1
  }
}

maxN_coeff = 10
nSig = 5

aucs_all = matrix(0,maxN_coeff,nSig)

for (nLargestCoeffs in 1:maxN_coeff){
  
  print(nLargestCoeffs)
  
  sigScores = matrix(0,nPat,nSig)
  
  for (cSig in 1:nSig){
    
    
    if (cSig==1){
      # coeffs = read.csv("final_coeffs_S100B_sm.csv")
      coeffs = read.csv("final_coeffs_S100B_sm.csv")
      c_data = data_S100B
    }
    if (cSig==2){
      coeffs = read.csv("final_coeffs_CD45_sm.csv")
      c_data = data_CD45
    }
    if (cSig==3){
      coeffs = read.csv("final_coeffs_CD68_sm.csv")
      c_data = data_CD68
    }
    if (cSig==4){
      coeffs = read.csv("final_coeffs_pseudoStroma_sm.csv")
      c_data = data_pStroma
    }
    if (cSig==5){
      coeffs = read.csv("final_coeffs_pseudoBulk_sm.csv")
      c_data = data_pBulk
    }
    
    genes = coeffs$X[2:dim(coeffs)[1]]
    nGene = length(genes)
    cfs = coeffs$s0[2:dim(coeffs)[1]]
    
    cfs_sort = sort(abs(cfs),decreasing = TRUE)
    if ((min(nLargestCoeffs,nGene)+1) <= nGene){
      cfs_sort = cfs_sort[(min(nLargestCoeffs,nGene)+1):nGene]  
      cf_idx = which(abs(cfs) %in% cfs_sort)
      cfs[cf_idx] = 0    
    }
    
    for (cPat in 1:nPat){
      vals = rep(0,nGene)
      pat_idx = which(rownames(c_data)==patIds[cPat])
      if(length(pat_idx)>0){
        for (i in 1:nGene){
          idx = which(colnames(c_data)==genes[i])
          if (length(idx)>0){
            if (is.finite(c_data[pat_idx,idx])){
              vals[i] = c_data[pat_idx,idx]
            }
          }else{
            if(cPat==1){
              #print(genes)
              print(genes[i])
            }
          }
        }
        sigScores[cPat,cSig] = sum((as.matrix(vals))*(as.matrix(cfs)))
      }else{
        sigScores[cPat,cSig] = NA
      }
      if(cPat==35){
        print(c(sum(vals!=0),nGene)) 
      }
    } 
    #print(c(sum(vals!=0),nGene))  
  }
  
  #sigScores0 = sigScores
  #sigScores2 = scale(sigScores)
  ##sigScores2[is.na(sigScores2)] = 0
  #sigScores = sigScores2
  
  idx = which(!is.na(sigScores[,1]))
  pROC_obj <- roc(resps[idx],sigScores[idx,1],
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)
  pROC_obj
  aucs_all[nLargestCoeffs,1] = pROC_obj$auc
  
  idx = which(!is.na(sigScores[,2]))
  pROC_obj <- roc(resps[idx],sigScores[idx,2],
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)
  pROC_obj
  aucs_all[nLargestCoeffs,2] = pROC_obj$auc
  
  idx = which(!is.na(sigScores[,3]))
  pROC_obj <- roc(resps[idx],sigScores[idx,3],
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)
  pROC_obj
  aucs_all[nLargestCoeffs,3] = pROC_obj$auc
  
  idx = which(!is.na(sigScores[,4]))
  pROC_obj <- roc(resps[idx],sigScores[idx,4],
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)
  pROC_obj
  aucs_all[nLargestCoeffs,4] = pROC_obj$auc
  
  
  idx = which(!is.na(sigScores[,5]))
  pROC_obj <- roc(resps[idx],sigScores[idx,5],
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)
  pROC_obj
  aucs_all[nLargestCoeffs,5] = pROC_obj$auc
  
  
}

print(aucs_all)
write.csv(aucs_all, "aucs_all_DiscoveryCohort_smSig.csv")

