library("survminer")
library("ggplot2")
require("survival")
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(edgeR)
library(pROC)

# get Sig DE genes for training set 

S100B_df = read.csv("S100B_av.csv", header = T)
S100B_df$SpotID <- S100B_df$X
Res_S100B = read.csv("S100B_response_table.csv", header = T)
Res_S100B$SpotID <- Res_S100B$ROI_ID_2
S100B_2 = merge(S100B_df, Res_S100B, by = "SpotID", all = FALSE)
S100B_2$group <- S100B_2$response # define the group you want to compare
data = S100B_2[,3:9286]
prop_test = 0.20
idx0 = which(S100B_2$group=="no")
idx1 = which(S100B_2$group=="yes")

aucs = read.csv("aucs_S100B_300.csv")

good_models = which((aucs$x)>=0.7)
nGoodModels = length(good_models)
goodModelGeneCounts = rep(0,dim(data)[2])
names(goodModelGeneCounts) = colnames(data)

for (cSplit in good_models){
  
  print("current split:")
  print(cSplit)
  print("")
  
  set.seed(cSplit)
  idx0_shuffled = sample(idx0, length(idx0), replace = FALSE)
  idx1_shuffled = sample(idx1, length(idx1), replace = FALSE)
  
  test_idx = c(idx0_shuffled[1:(length(idx0)*prop_test)],idx1_shuffled[1:(length(idx1)*prop_test)])
  train_idx = setdiff(1:dim(data)[1],test_idx)
  train_data = S100B_2[train_idx,]
  train_data2 = train_data[,3:9286]
  test_data = S100B_2[test_idx,]
  test_data2 = test_data[,3:9286]
  
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
  write.csv(response$table, "S100B_train_set_Sig_Genes.csv")
  
  dex_df = response$table
  dex_sel_genes = dex_df[which(dex_df$FDR < 0.2),] # take only genes that have FDR < 0.1
  dex_df$gene = row.names(dex_df)
  do_lasso = 1
  
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
      # plot(lasso_cv)
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
  
  test_data3 = test_data[,which(colnames(test_data) %in% rownames(dex_sel_genes))]
  test_data3$response = test_data$response
  test_data3$lasso_scores = matrix(0,dim(test_data3)[1],1)
  for(i in 1:dim(test_data3)[1]){
    test_data3$lasso_scores[i] = sum(as.matrix(test_data3[i,nonzero])*t(coeffs[2:dim(coeffs)[1],]))
  }
  
  pROC_obj <- roc(test_data3$response,test_data3$lasso_scores,
                  smoothed = TRUE,
                  # arguments for ci
                  ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                  # arguments for plot
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE)
  
  
  print(c(pROC_obj$auc,aucs$x[cSplit]))
  
  gene_idxs = which(names(goodModelGeneCounts) %in% rownames(coeffs)[2:length(rownames(coeffs))])
  goodModelGeneCounts[gene_idxs] = goodModelGeneCounts[gene_idxs] + 1
  
}

write.csv(goodModelGeneCounts,"goodModelCounts_S100B_from300.csv")
#goodModelGeneCounts = read.csv("goodModelCounts_S100B.csv",row.names = 1)
goodModelGeneCounts_ = as.matrix(goodModelGeneCounts)
names(goodModelGeneCounts_) = row.names(goodModelGeneCounts)
goodModelGeneCounts = goodModelGeneCounts_

goodModelGeneProps = goodModelGeneCounts / nGoodModels
names(goodModelGeneProps) = row.names(goodModelGeneCounts)
goodModelGeneProps = goodModelGeneProps[goodModelGeneProps>0]
sort(goodModelGeneProps)

final_genes = names(goodModelGeneProps)[which(goodModelGeneProps>0.04)]

train_data = S100B_2
train_data2 = train_data[,3:9286]
train_data3 = train_data[,which(colnames(train_data) %in% final_genes)]
train_data3$response = train_data$response

cSeed = 3001
print(cSeed)
set.seed(cSeed)

model_cv <- glmnet(as.matrix(train_data3[,1:length(final_genes)]), train_data3$response, alpha = 0.5, lambda = 0.00001, standardize = FALSE, family = "binomial")
final_coeffs <- predict(model_cv, as.matrix(train_data3), type = "coefficients")

print(final_coeffs)


train_data3$lasso_scores = matrix(0,dim(train_data3)[1],1)
for(i in 1:dim(train_data3)[1]){
  train_data3$lasso_scores[i] = sum(as.matrix(train_data3[i,1:length(final_genes)])*t(final_coeffs[2:dim(final_coeffs)[1],]))
}

## AUC analysis

#data3 %>%
#  roc(response, lasso_score)

pROC_obj <- roc(train_data3$response,train_data3$lasso_scores,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE, direction = '<')

pROC_obj





