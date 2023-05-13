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
S100B_2 = merge(S100B_df, Res_S100B, by = "SpotID")
S100B_2$group <- S100B_2$response # define the group you want to compare
data = S100B_2[,3:9286]
prop_test = 0.33
idx0 = which(S100B_2$group=="no")
idx1 = which(S100B_2$group=="yes")
test_idx = c(idx0[1:(length(idx0)*prop_test)],idx1[1:(length(idx1)*prop_test)])
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
sum(FDR < 0.1)
topTags(response)
plotMD(response)
abline(h=c(-1,1), col="blue")
response$table$FDR = FDR
write.csv(response$table, "S100B_train_set_Sig_Genes.csv")

# Fibroblast_genes = read.csv("Fibroblast.csv")# Fibroblast signature gene
# dex_df = read.csv("S100B_train_set_Sig_Genes.csv")
# names(dex_df)[1] = "Genes"
# dex_df = dex_df[which(dex_df$Genes %in% Fibroblast_genes$x),] # gene dex genes from the S100 counts file.
# dex_df = dex_df[order(dex_df$PValue),]
# dex_df = dex_df[which(dex_df$logFC<0),] # take only up regulated genes. 
# dex_sel_genes = dex_df$Genes[1:20]
# do_lasso = 0


dex_df = response$table
dex_df = dex_df[which(dex_df$FDR < 0.05),] # take only genes that have FDR < 0.05
dex_df$gene = row.names(dex_df)
dex_sel_genes = dex_df$gene
do_lasso = 0

dex_sel_genes = dex_df[which(dex_df$FDR < 0.1),] # take only genes that have FDR < 0.1
dex_df$gene = row.names(dex_df)
dex_sel_genes = dex_df$gene
do_lasso = 1

train_data3 = train_data[,which(colnames(train_data) %in% dex_sel_genes)]
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
    lambdas_to_try <- 10^seq(-3, 0, length.out = 100)
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

coeffs

# 20 x 1 sparse Matrix of class "dgCMatrix"
# s0
# (Intercept) 38.709388094
# IGFBP7      -0.028325549
# BCAN         0.017471517
# DCT          0.017975809
# EDNRB       -0.032674313
# ERBB3       -0.099629885
# PLAT         0.001009160
# YWHAE       -0.096888640
# S100A4       0.015941932
# CTSL         0.022249899
# IGFBP5      -0.015113786
# MYC         -0.112831852
# APOD         0.008093680
# H4C12        0.009889335
# B2M          0.009690327
# CD74        -0.003181190
# MT2A        -0.008360141
# HLA.C        0.019328614
# HSPA1B       0.061612741
# MYH9        -0.003054007

## heatmap

# Training set 
train_data4 = train_data3[,nonzero]
train_data4 <- log2(train_data4)
rownames(train_data4) = sapply(rownames(train_data4),function(x) 
  strsplit(as.character(x),split = "\\\\")[[1]][1])
train_data4 <- train_data4[order(train_data$response),]
train_data5 = data.frame("Response" = sort(train_data$response))
rownames(train_data5) = rownames(train_data4)
#rownames(train_data4) = train_data5$Response # name matching
pheatmap(train_data4,annotation_row = train_data5, main = "Signature genes S100B Comp",
         cluster_rows = T, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
         # cutree_cols = 2,
         cutree_rows = 2, 
         fontsize = 8)
pheatmap(train_data4,annotation_row = train_data5, main = "Signature genes S100B Comp",
         cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
         cutree_cols = 2,
         #cutree_rows = 3, 
         fontsize = 8)

## heatmap

# Testing set 
test_data = S100B_2[test_idx,]
test_data2 = test_data[,which(colnames(test_data) %in% dex_sel_genes)]
test_data2 = test_data2[,nonzero]
test_data2 <- log2(test_data2)
rownames(test_data2) = sapply(rownames(test_data2),function(x) 
  strsplit(as.character(x),split = "\\\\")[[1]][1])
test_data2 <- test_data2[order(test_data$response),]
test_data3 = data.frame("Response" = sort(test_data$response))
rownames(test_data3) = rownames(test_data2)
#rownames(test_data3) = train_data5$Response # name matching
pheatmap(test_data2,annotation_row = test_data3, main = "Signature genes S100B Comp",
         cluster_rows = T, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
         # cutree_cols = 2,
         cutree_rows = 2, 
         clustering_method = "complete",
         # clustering_distance_cols = "correlation" or "euclidean", ,
         # clustering_distance_rows = "correlation",
         fontsize = 8)

pheatmap(test_data2,annotation_row = test_data3, main = "Signature genes S100B Comp",
         cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
         cutree_cols = 2,
         #cutree_rows = 3, 
         fontsize = 8, clustering_method = "complete")

#data3$lasso_score <- predict(model_cv, as.matrix(data3[,nonzero]))

S100B_2$lasso_scores = matrix(0,dim(S100B_2)[1],1)
for(i in 1:dim(S100B_2)[1]){
  S100B_2$lasso_scores[i] = sum(as.matrix(S100B_2[i,nonzero])*t(coeffs[2:dim(coeffs)[1],]))
}

## AUC analysis

#data3 %>%
#  roc(response, lasso_score)

pROC_obj <- roc(S100B_2$response[test_idx],S100B_2$lasso_scores[test_idx],
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

pROC_obj
