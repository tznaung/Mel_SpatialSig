library("survminer")
library("ggplot2")
require("survival")
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(edgeR)
library(pROC)

# get data for all compartments 

CD68 = read.csv("Q3_Norm_10%_CD68.csv", header = T)
CD45 = read.csv("Q3_Norm_10%_CD45.csv", header = T)
S100B = read.csv("Q3_Norm_10%_S100B.csv", header = T)

#get the average counts of both blocks and combine CD45 and CD68 to get bulk compartment
dim1 = dim(CD45)[1]
dim2 = dim(CD45)[2]
ids_short = rep(0,dim1)
for(i in 1:dim1){
  nm = CD45[i,2]
  ids_short[i] = nm
}
unq = unique(ids_short)
CD45_av <- data.frame(matrix(ncol = dim2-2, nrow = 0))
for(i in 1:length(unq)){
  vec = colMeans(CD45[which(ids_short==unq[i]),3:dim2])
  CD45_av = rbind(CD45_av,vec)
}
rownames(CD45_av) = unq
colnames(CD45_av) = colnames(CD45)[3:dim2]

dim1 = dim(CD68)[1]
dim2 = dim(CD68)[2]
ids_short = rep(0,dim1)
for(i in 1:dim1){
  nm = CD68[i,2]
  ids_short[i] = nm
}
unq = unique(ids_short)
CD68_av <- data.frame(matrix(ncol = dim2-2, nrow = 0))
for(i in 1:length(unq)){
  vec = colMeans(CD68[which(ids_short==unq[i]),3:dim2])
  CD68_av = rbind(CD68_av,vec)
}
rownames(CD68_av) = unq
colnames(CD68_av) = colnames(CD68)[3:dim2]

dim1 = dim(S100B)[1]
dim2 = dim(S100B)[2]
ids_short = rep(0,dim1)
for(i in 1:dim1){
  nm = S100B[i,2]
  ids_short[i] = nm
}
unq = unique(ids_short)
S100B_av <- data.frame(matrix(ncol = dim2-2, nrow = 0))
for(i in 1:length(unq)){
  vec = colMeans(S100B[which(ids_short==unq[i]),3:dim2])
  S100B_av = rbind(S100B_av,vec)
}
rownames(S100B_av) = unq
colnames(S100B_av) = colnames(S100B)[3:dim2]

# for(i in 1:dim(CD45_av)[1]){  
#   CD45_av[i,] = 1000000 * CD45_av[i,] / sum(CD45_av[i,])
# }
# for(i in 1:dim(CD68_av)[1]){  
#   CD68_av[i,] = 1000000 * CD68_av[i,] / sum(CD68_av[i,])
# }
# for(i in 1:dim(S100B_av)[1]){  
#   S100B_av[i,] = 1000000 * S100B_av[i,] / sum(S100B_av[i,])
# }

## combine 3 compartment (CD45+CD68+S100B) tables 
allnames = union(colnames(S100B_av), colnames(CD45_av))
allnames = union(allnames, colnames(CD68_av))
CD45_av2 = matrix(0,dim(CD45_av)[1],length(allnames))
f <- function(x) which(allnames==x)
idxs = lapply(colnames(CD45_av),f)
for(i in 1:length(idxs)){
  CD45_av2[,unlist(idxs)[i]] = CD45_av[,i]
}
CD68_av2 = matrix(0,dim(CD68_av)[1],length(allnames))
f <- function(x) which(allnames==x)
idxs = lapply(colnames(CD68_av),f)
for(i in 1:length(idxs)){
  CD68_av2[,unlist(idxs)[i]] = CD68_av[,i]
}
S100B_av2 = matrix(0,dim(S100B_av)[1],length(allnames))
f <- function(x) which(allnames==x)
idxs = lapply(colnames(S100B_av),f)
for(i in 1:length(idxs)){
  S100B_av2[,unlist(idxs)[i]] = S100B_av[,i]
}

allnames2 = intersect(colnames(CD68_av), colnames(CD45_av))
allnames2 = intersect(allnames2, colnames(S100B_av))
idx = which(allnames %in% allnames2)
CD45_av2 = CD45_av2[,idx]
CD68_av2 = CD68_av2[,idx]
S100B_av2 = S100B_av2[,idx]
allnames =  allnames[idx]

pseudo_bulk = CD45_av2
subj_names = rownames(CD45_av)
nGenes_total = length(allnames)
for(i in 1:dim(CD68_av2)[1]){
  idx = which(rownames(CD68_av)[i]==subj_names)
  if(length(idx)>0){
    for(j in 1:nGenes_total){
      pseudo_bulk[idx,j] = pseudo_bulk[idx,j] + CD68_av2[i,j]      
    }
  }else{
    pseudo_bulk = rbind(pseudo_bulk,CD68_av2[i,])
    subj_names = append(subj_names,rownames(CD68_av)[i])
  }
}
for(i in 1:dim(S100B_av2)[1]){
  idx = which(rownames(S100B_av)[i]==subj_names)
  if(length(idx)>0){
    for(j in 1:nGenes_total){
      pseudo_bulk[idx,j] = pseudo_bulk[idx,j] + S100B_av2[i,j]      
    }
  }else{
    pseudo_bulk = rbind(pseudo_bulk,S100B_av2[i,])
    subj_names = append(subj_names,rownames(S100B_av)[i])
  }
}

# for(i in 1:dim(pseudo_bulk)[1]){  
#   pseudo_bulk[i,] = 1000000 * pseudo_bulk[i,] / sum(pseudo_bulk[i,])
# }
# nSubj = length(subj_names)
# 
# pseudo_bulk_norm = log2(pseudo_bulk + 1)

pseudo_bulk_df = as.data.frame(pseudo_bulk)
colnames(pseudo_bulk_df) = allnames
rownames(pseudo_bulk_df) = subj_names
pseudo_bulk_df$SpotID = rownames(pseudo_bulk_df)

# matching patient's ID with the response.

Res_Bulk = read.csv("pseudoBulk_response_table.csv", header = T)
Res_Bulk$SpotID <- Res_Bulk$ROI_ID_2
Bulk_df = merge(pseudo_bulk_df, Res_Bulk, by = "SpotID", all = FALSE)

Bulk_df$group <- Bulk_df$response # define the group you want to compare
data = Bulk_df[,2:1631]
prop_test = 0.20
idx0 = which(Bulk_df$group=="no")
idx1 = which(Bulk_df$group=="yes")

#number_of_splits = 300

aucs = rep(-1,300)

for (cSplit in 1:300){
  
  print(cSplit)
  set.seed(cSplit)
  idx0_shuffled = sample(idx0, length(idx0), replace = FALSE)
  idx1_shuffled = sample(idx1, length(idx1), replace = FALSE)
  
  test_idx = c(idx0_shuffled[1:(length(idx0)*prop_test)],idx1_shuffled[1:(length(idx1)*prop_test)])
  train_idx = setdiff(1:dim(data)[1],test_idx)
  train_data = Bulk_df[train_idx,]
  train_data2 = train_data[,3:1713]
  test_data = Bulk_df[test_idx,]
  test_data2 = test_data[,3:1713]
  
  # make the pseudo counts of the data
  
  d = DGEList(t(train_data2), group= train_data$response)
  d = calcNormFactors(d)    
  d = estimateCommonDisp(d)
  d = estimateTagwiseDisp(d)
  
  #Find significant genes
  response <- exactTest(d, pair = c("no", "yes"))
  FDR <- p.adjust(response$table$PValue, method="BH")
  sum(FDR < 0.2)
  topTags(response)
  plotMD(response)
  abline(h=c(-1,1), col="blue")
  response$table$FDR = FDR
  
  dex_df = response$table
  dex_sel_genes = dex_df[which(dex_df$FDR < 0.2),] # take only genes that have FDR < 0.1
  dex_df$gene = row.names(dex_df)
  do_lasso = 1
  
  if (dim(dex_sel_genes)[1]<=1){
    next
  }
  
  train_data3 = train_data[,which(colnames(train_data) %in% rownames(dex_sel_genes))]
  train_data3$response = train_data$response
  
  if (do_lasso){
    nGene = dim(train_data3)[2] - 1
    dim2 = dim(train_data3)[2]
    var_counts = matrix(0,1,nGene+1)
    mean_coeffs = matrix(0,nGene+1,1)
    nSeeds = 100
    
    for(cSeed in 1:nSeeds){
      
      print(cSeed)
      set.seed(cSeed)
      
      # lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
      lambdas_to_try <- 10^seq(-3, 2, length.out = 100)
      ##
      lasso_cv <- cv.glmnet(as.matrix(train_data3[,1:nGene]), train_data3$response, alpha = 1, lambda = lambdas_to_try,
                            standardize = FALSE, nfolds = 10, family = "binomial")
      plot(lasso_cv)
      lambda_cv <- lasso_cv$lambda.min
      model_cv <- glmnet(as.matrix(train_data3[,1:nGene]), train_data3$response, alpha = 1, lambda = lambda_cv, standardize = FALSE, family = "binomial")
      coeffs <- predict(model_cv, as.matrix(train_data3[,1:nGene]), type = "coefficients")
      
      var_counts[which(coeffs!=0)] = var_counts[which(coeffs!=0)] + 1
      mean_coeffs = mean_coeffs + (coeffs/nSeeds)
      
    }
    
    mean_coeffs_weighted = mean_coeffs * (nSeeds/t(var_counts));
    
    var_props = var_counts / nSeeds
    colnames(var_props) = c("intercept",colnames(train_data3)[1:nGene])
    rownames(mean_coeffs) = c("intercept",colnames(train_data3)[1:nGene])
    rownames(mean_coeffs_weighted) = c("intercept",colnames(train_data3)[1:nGene])
    t(var_props)
    #mean_coeffs
    mean_coeffs_weighted #[2:nGene+1,]
    
    ## refit model
    var_props_nointerc = var_props[1,2:(nGene+1)]
    nonzero = which(var_props_nointerc>0)
    
    cSeed = 2001
    print(cSeed)
    set.seed(cSeed)
    
    model_cv <- glmnet(as.matrix(train_data3[,nonzero]), train_data3$response, alpha = 1, lambda = 0, standardize = FALSE, family = "binomial")
    coeffs <- predict(model_cv, as.matrix(train_data3[,nonzero]), type = "coefficients")
    
    coeffs
  }else{
    nGene = dim(train_data3)[2] - 1
    dim2 = dim(train_data3)[2]  
    nonzero = (1:nGene)
  }
  
  print(coeffs)
  
  # Testing set 
  test_data = Bulk_df[test_idx,]
  test_data3 = test_data[,which(colnames(test_data) %in% rownames(dex_sel_genes))]
  test_data3$response = test_data$response  
  test_data3$lasso_scores = matrix(0,dim(test_data3)[1],1)
  for(i in 1:dim(test_data3)[1]){
    test_data3$lasso_scores[i] = sum(as.matrix(test_data3[i,nonzero])*t(coeffs[2:dim(coeffs)[1],]))
  }  
  
  ## AUC analysis
  
  
  pROC_obj <- roc(test_data3$response,test_data3$lasso_scores,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE, direction = '<')
  
  pROC_obj
  
  #dev.off()
  
  aucs[cSplit] = pROC_obj$auc
  write.csv(aucs, "aucs_PBulk_300.csv")  
  
}

