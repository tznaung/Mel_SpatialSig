library(pROC)
library(glmnet)
library(RColorBrewer)
library(ggplot2)

# load signatures

S100B_Sig_sm = read.csv("final_coeffs_S100B_sm.csv")
CD45_Sig_sm = read.csv("final_coeffs_CD45_sm.csv")
CD68_Sig_sm = read.csv("final_coeffs_CD68_sm.csv")
pBulk_Sig_sm = read.csv("final_coeffs_pseudoBulk_sm.csv")

# load data to test

S100B_data = read.csv("S100B_df.csv", row.names = 2)
S100B_df = S100B_data[,3:9286]
CD45_data = read.csv("CD45_df.csv", row.names = 2)
CD45_df = CD45_data[,3:3738]
CD68_data = read.csv("CD68_df.csv", row.names = 2)
CD68_df = CD68_data[,3:1713]
pBulk_data = read.csv("PseudoBulk_df.csv", row.names = 1)
pBulk_df = pBulk_data[,1:7320]
pStroma_data = read.csv("PseudoStroma_df.csv", row.names = 2)
pStroma_df = pStroma_data[,2:2336]

###############################################################

# generate lasso scores for each signatures in each dataset
# S100B signature in S100B compartment

S100B_nGene = dim(S100B_Sig_sm)[1] - 1 # remove first row
S100B_final_genes = S100B_Sig_sm$X[2:(S100B_nGene+1)]
S100B_final_coeffs = S100B_Sig_sm$s0[2:(S100B_nGene+1)]

S100B_S100B_idx = rep(0,S100B_nGene)
for(i in 1:S100B_nGene){
  idx = which(colnames(S100B_df)==S100B_final_genes[i])
  if(length(idx)>0){
    S100B_S100B_idx[i] = idx
  }
}

S100B_S100B_final_genes = S100B_final_genes[which(S100B_S100B_idx>0)]
S100B_S100B_final_coeffs = S100B_final_coeffs[which(S100B_S100B_idx>0)]
S100B_S100B_nGene = length(S100B_S100B_final_genes)

## AUC analysis

S100B_S100B_AUC = S100B_df[,S100B_S100B_idx[S100B_S100B_idx>0]]
S100B_S100B_AUC$response = S100B_data$response

S100B_S100B_AUC$lasso_scores = matrix(0,dim(S100B_S100B_AUC)[1],1)
for(i in 1:dim(S100B_S100B_AUC)[1]){
  S100B_S100B_AUC$lasso_scores[i] = sum(as.matrix(S100B_S100B_AUC[i,1:S100B_S100B_nGene])*S100B_S100B_final_coeffs[1:S100B_S100B_nGene])
}

pROC_S100B_S100B <- roc(S100B_S100B_AUC$response,S100B_S100B_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_S100B_S100B

###############################################################

# generate lasso scores for each signatures in each dataset
# S100B signature in CD45 compartment

S100B_nGene = dim(S100B_Sig_sm)[1] - 1 # remove first row
S100B_final_genes = S100B_Sig_sm$X[2:(S100B_nGene+1)]
S100B_final_coeffs = S100B_Sig_sm$s0[2:(S100B_nGene+1)]

CD45_S100B_idx = rep(0,S100B_nGene)
for(i in 1:S100B_nGene){
  idx = which(colnames(CD45_df)==S100B_final_genes[i])
  if(length(idx)>0){
    CD45_S100B_idx[i] = idx
  }
}

CD45_S100B_final_genes = S100B_final_genes[which(CD45_S100B_idx>0)]
CD45_S100B_final_coeffs = S100B_final_coeffs[which(CD45_S100B_idx>0)]
CD45_S100B_nGene = length(CD45_S100B_final_genes)

## AUC analysis

CD45_S100B_AUC = CD45_df[,CD45_S100B_idx[CD45_S100B_idx>0]]
CD45_S100B_AUC$response = CD45_data$response

CD45_S100B_AUC$lasso_scores = matrix(0,dim(CD45_S100B_AUC)[1],1)
for(i in 1:dim(CD45_S100B_AUC)[1]){
  CD45_S100B_AUC$lasso_scores[i] = sum(as.matrix(CD45_S100B_AUC[i,1:CD45_S100B_nGene])*CD45_S100B_final_coeffs[1:CD45_S100B_nGene])
}

pROC_CD45_S100B <- roc(CD45_S100B_AUC$response,CD45_S100B_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_CD45_S100B

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# S100B in CD68 compartment

CD68_S100B_idx = rep(0,S100B_nGene)
for(i in 1:S100B_nGene){
  idx = which(colnames(CD68_df)==S100B_final_genes[i])
  if(length(idx)>0){
    CD68_S100B_idx[i] = idx
  }
}

CD68_S100B_final_genes = S100B_final_genes[which(CD68_S100B_idx>0)]
CD68_S100B_final_coeffs = S100B_final_coeffs[which(CD68_S100B_idx>0)]
CD68_S100B_nGene = length(CD68_S100B_final_genes)

## AUC analysis

CD68_S100B_AUC = CD68_df[,CD68_S100B_idx[CD68_S100B_idx>0]]
CD68_S100B_AUC$response = CD68_data$response

CD68_S100B_AUC$lasso_scores = matrix(0,dim(CD68_S100B_AUC)[1],1)
for(i in 1:dim(CD68_S100B_AUC)[1]){
  CD68_S100B_AUC$lasso_scores[i] = sum(as.matrix(CD68_S100B_AUC[i,1:CD68_S100B_nGene])*CD68_S100B_final_coeffs[1:CD68_S100B_nGene])
}

pROC_CD68_S100B <- roc(CD68_S100B_AUC$response,CD68_S100B_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_CD68_S100B

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# S100B in pseudoBulk 

pBulk_S100B_idx = rep(0,S100B_nGene)
for(i in 1:S100B_nGene){
  idx = which(colnames(pBulk_df)==S100B_final_genes[i])
  if(length(idx)>0){
    pBulk_S100B_idx[i] = idx
  }
}

pBulk_S100B_final_genes = S100B_final_genes[which(pBulk_S100B_idx>0)]
pBulk_S100B_final_coeffs = S100B_final_coeffs[which(pBulk_S100B_idx>0)]
pBulk_S100B_nGene = length(pBulk_S100B_final_genes)

## AUC analysis

pBulk_S100B_AUC = pBulk_df[,pBulk_S100B_idx[pBulk_S100B_idx>0]]
pBulk_S100B_AUC$response = pBulk_data$response

pBulk_S100B_AUC$lasso_scores = matrix(0,dim(pBulk_S100B_AUC)[1],1)
for(i in 1:dim(pBulk_S100B_AUC)[1]){
  pBulk_S100B_AUC$lasso_scores[i] = sum(as.matrix(pBulk_S100B_AUC[i,1:pBulk_S100B_nGene])*pBulk_S100B_final_coeffs[1:pBulk_S100B_nGene])
}

pROC_pBulk_S100B <- roc(pBulk_S100B_AUC$response,pBulk_S100B_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_pBulk_S100B

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# S100B in pseudoStroma 

pStroma_S100B_idx = rep(0,S100B_nGene)
for(i in 1:S100B_nGene){
  idx = which(colnames(pStroma_df)==S100B_final_genes[i])
  if(length(idx)>0){
    pStroma_S100B_idx[i] = idx
  }
}

pStroma_S100B_final_genes = S100B_final_genes[which(pStroma_S100B_idx>0)]
pStroma_S100B_final_coeffs = S100B_final_coeffs[which(pStroma_S100B_idx>0)]
pStroma_S100B_nGene = length(pStroma_S100B_final_genes)

## AUC analysis

pStroma_S100B_AUC = pStroma_df[,pStroma_S100B_idx[pStroma_S100B_idx>0]]
pStroma_S100B_AUC$response = pStroma_data$response

pStroma_S100B_AUC$lasso_scores = matrix(0,dim(pStroma_S100B_AUC)[1],1)
for(i in 1:dim(pStroma_S100B_AUC)[1]){
  pStroma_S100B_AUC$lasso_scores[i] = sum(as.matrix(pStroma_S100B_AUC[i,1:pStroma_S100B_nGene])*pStroma_S100B_final_coeffs[1:pStroma_S100B_nGene])
}

pROC_pStroma_S100B <- roc(pStroma_S100B_AUC$response,pStroma_S100B_AUC$lasso_scores,
                        smoothed = TRUE,
                        # arguments for ci
                        ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                        # arguments for plot
                        plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                        print.auc=TRUE, show.thres=TRUE)

pROC_pStroma_S100B

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD45 signature in CD45 compartment 

CD45_nGene = dim(CD45_Sig_sm)[1] - 1 # remove first row
CD45_final_genes = CD45_Sig_sm$X[2:(CD45_nGene+1)]
CD45_final_coeffs = CD45_Sig_sm$s0[2:(CD45_nGene+1)]

CD45_CD45_idx = rep(0,CD45_nGene)
for(i in 1:CD45_nGene){
  idx = which(colnames(CD45_df)==CD45_final_genes[i])
  if(length(idx)>0){
    CD45_CD45_idx[i] = idx
  }
}

CD45_CD45_final_genes = CD45_final_genes[which(CD45_CD45_idx>0)]
CD45_CD45_final_coeffs = CD45_final_coeffs[which(CD45_CD45_idx>0)]
CD45_CD45_nGene = length(CD45_CD45_final_genes)

## AUC analysis

CD45_CD45_AUC = CD45_df[,CD45_CD45_idx[CD45_CD45_idx>0]]
CD45_CD45_AUC$response = CD45_data$response

CD45_CD45_AUC$lasso_scores = matrix(0,dim(CD45_CD45_AUC)[1],1)
for(i in 1:dim(CD45_CD45_AUC)[1]){
  CD45_CD45_AUC$lasso_scores[i] = sum(as.matrix(CD45_CD45_AUC[i,1:CD45_CD45_nGene])*CD45_CD45_final_coeffs[1:CD45_CD45_nGene])
}

pROC_CD45_CD45 <- roc(CD45_CD45_AUC$response,CD45_CD45_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_CD45_CD45

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD45 signature in S100B compartment 

CD45_nGene = dim(CD45_Sig_sm)[1] - 1 # remove first row
CD45_final_genes = CD45_Sig_sm$X[2:(CD45_nGene+1)]
CD45_final_coeffs = CD45_Sig_sm$s0[2:(CD45_nGene+1)]

S100B_CD45_idx = rep(0,CD45_nGene)
for(i in 1:CD45_nGene){
  idx = which(colnames(S100B_df)==CD45_final_genes[i])
  if(length(idx)>0){
    S100B_CD45_idx[i] = idx
  }
}

S100B_CD45_final_genes = CD45_final_genes[which(S100B_CD45_idx>0)]
S100B_CD45_final_coeffs = CD45_final_coeffs[which(S100B_CD45_idx>0)]
S100B_CD45_nGene = length(S100B_CD45_final_genes)

## AUC analysis

S100B_CD45_AUC = S100B_df[,S100B_CD45_idx[S100B_CD45_idx>0]]
S100B_CD45_AUC$response = S100B_data$response

S100B_CD45_AUC$lasso_scores = matrix(0,dim(S100B_CD45_AUC)[1],1)
for(i in 1:dim(S100B_CD45_AUC)[1]){
  S100B_CD45_AUC$lasso_scores[i] = sum(as.matrix(S100B_CD45_AUC[i,1:S100B_CD45_nGene])*S100B_CD45_final_coeffs[1:S100B_CD45_nGene])
}

pROC_S100B_CD45 <- roc(S100B_CD45_AUC$response,S100B_CD45_AUC$lasso_scores,
                          smoothed = TRUE,
                          # arguments for ci
                          ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                          # arguments for plot
                          plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                          print.auc=TRUE, show.thres=TRUE)

pROC_S100B_CD45

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD45 signature in CD68 compartment 

CD68_CD45_idx = rep(0,CD45_nGene)
for(i in 1:CD45_nGene){
  idx = which(colnames(CD68_df)==CD45_final_genes[i])
  if(length(idx)>0){
    CD68_CD45_idx[i] = idx
  }
}

CD68_CD45_final_genes = CD45_final_genes[which(CD68_CD45_idx>0)]
CD68_CD45_final_coeffs = CD45_final_coeffs[which(CD68_CD45_idx>0)]
CD68_CD45_nGene = length(CD68_CD45_final_genes)

## AUC analysis

CD68_CD45_AUC = CD68_df[,CD68_CD45_idx[CD68_CD45_idx>0]]
CD68_CD45_AUC$response = CD68_data$response

CD68_CD45_AUC$lasso_scores = matrix(0,dim(CD68_CD45_AUC)[1],1)
for(i in 1:dim(CD68_CD45_AUC)[1]){
  CD68_CD45_AUC$lasso_scores[i] = sum(as.matrix(CD68_CD45_AUC[i,1:CD68_CD45_nGene])*CD68_CD45_final_coeffs[1:CD68_CD45_nGene])
}

pROC_CD68_CD45 <- roc(CD68_CD45_AUC$response,CD68_CD45_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_CD68_CD45

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD45 signature in pseudoBulk 

pBulk_CD45_idx = rep(0,CD45_nGene)
for(i in 1:CD45_nGene){
  idx = which(colnames(pBulk_df)==CD45_final_genes[i])
  if(length(idx)>0){
    pBulk_CD45_idx[i] = idx
  }
}

pBulk_CD45_final_genes = CD45_final_genes[which(pBulk_CD45_idx>0)]
pBulk_CD45_final_coeffs = CD45_final_coeffs[which(pBulk_CD45_idx>0)]
pBulk_CD45_nGene = length(pBulk_CD45_final_genes)

## AUC analysis

pBulk_CD45_AUC = pBulk_df[,pBulk_CD45_idx[pBulk_CD45_idx>0]]
pBulk_CD45_AUC$response = pBulk_data$response

pBulk_CD45_AUC$lasso_scores = matrix(0,dim(pBulk_CD45_AUC)[1],1)
for(i in 1:dim(pBulk_CD45_AUC)[1]){
  pBulk_CD45_AUC$lasso_scores[i] = sum(as.matrix(pBulk_CD45_AUC[i,1:pBulk_CD45_nGene])*pBulk_CD45_final_coeffs[1:pBulk_CD45_nGene])
}

pROC_pBulk_CD45 <- roc(pBulk_CD45_AUC$response,pBulk_CD45_AUC$lasso_scores,
                      smoothed = TRUE,
                      # arguments for ci
                      ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                      # arguments for plot
                      plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                      print.auc=TRUE, show.thres=TRUE)

pROC_pBulk_CD45

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD45 signature in pseudoBulk 

pStroma_CD45_idx = rep(0,CD45_nGene)
for(i in 1:CD45_nGene){
  idx = which(colnames(pStroma_df)==CD45_final_genes[i])
  if(length(idx)>0){
    pStroma_CD45_idx[i] = idx
  }
}

pStroma_CD45_final_genes = CD45_final_genes[which(pStroma_CD45_idx>0)]
pStroma_CD45_final_coeffs = CD45_final_coeffs[which(pStroma_CD45_idx>0)]
pStroma_CD45_nGene = length(pStroma_CD45_final_genes)

## AUC analysis

pStroma_CD45_AUC = pStroma_df[,pStroma_CD45_idx[pStroma_CD45_idx>0]]
pStroma_CD45_AUC$response = pStroma_data$response

pStroma_CD45_AUC$lasso_scores = matrix(0,dim(pStroma_CD45_AUC)[1],1)
for(i in 1:dim(pStroma_CD45_AUC)[1]){
  pStroma_CD45_AUC$lasso_scores[i] = sum(as.matrix(pStroma_CD45_AUC[i,1:pStroma_CD45_nGene])*pStroma_CD45_final_coeffs[1:pStroma_CD45_nGene])
}

pROC_pStroma_CD45 <- roc(pStroma_CD45_AUC$response,pStroma_CD45_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_pStroma_CD45

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD68 signature in CD68 compartment 

CD68_nGene = dim(CD68_Sig_sm)[1] - 1 # remove first row
CD68_final_genes = CD68_Sig_sm$X[2:(CD68_nGene+1)]
CD68_final_coeffs = CD68_Sig_sm$s0[2:(CD68_nGene+1)]

CD68_CD68_idx = rep(0,CD68_nGene)
for(i in 1:CD68_nGene){
  idx = which(colnames(CD68_df)==CD68_final_genes[i])
  if(length(idx)>0){
    CD68_CD68_idx[i] = idx
  }
}

CD68_CD68_final_genes = CD68_final_genes[which(CD68_CD68_idx>0)]
CD68_CD68_final_coeffs = CD68_final_coeffs[which(CD68_CD68_idx>0)]
CD68_CD68_nGene = length(CD68_CD68_final_genes)

## AUC analysis

CD68_CD68_AUC = CD68_df[,CD68_CD68_idx[CD68_CD68_idx>0]]
CD68_CD68_AUC$response = CD68_data$response

CD68_CD68_AUC$lasso_scores = matrix(0,dim(CD68_CD68_AUC)[1],1)
for(i in 1:dim(CD68_CD68_AUC)[1]){
  CD68_CD68_AUC$lasso_scores[i] = sum(as.matrix(CD68_CD68_AUC[i,1:CD68_CD68_nGene])*CD68_CD68_final_coeffs[1:CD68_CD68_nGene])
}

pROC_CD68_CD68 <- roc(CD68_CD68_AUC$response,CD68_CD68_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_CD68_CD68

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD68 signature in S100B compartment 

CD68_nGene = dim(CD68_Sig_sm)[1] - 1 # remove first row
CD68_final_genes = CD68_Sig_sm$X[2:(CD68_nGene+1)]
CD68_final_coeffs = CD68_Sig_sm$s0[2:(CD68_nGene+1)]

S100B_CD68_idx = rep(0,CD68_nGene)
for(i in 1:CD68_nGene){
  idx = which(colnames(S100B_df)==CD68_final_genes[i])
  if(length(idx)>0){
    S100B_CD68_idx[i] = idx
  }
}

S100B_CD68_final_genes = CD68_final_genes[which(S100B_CD68_idx>0)]
S100B_CD68_final_coeffs = CD68_final_coeffs[which(S100B_CD68_idx>0)]
S100B_CD68_nGene = length(S100B_CD68_final_genes)

## AUC analysis

S100B_CD68_AUC = S100B_df[,S100B_CD68_idx[S100B_CD68_idx>0]]
S100B_CD68_AUC$response = S100B_data$response

S100B_CD68_AUC$lasso_scores = matrix(0,dim(S100B_CD68_AUC)[1],1)
for(i in 1:dim(S100B_CD68_AUC)[1]){
  S100B_CD68_AUC$lasso_scores[i] = sum(as.matrix(S100B_CD68_AUC[i,1:S100B_CD68_nGene])*S100B_CD68_final_coeffs[1:S100B_CD68_nGene])
}

pROC_S100B_CD68 <- roc(S100B_CD68_AUC$response,S100B_CD68_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_S100B_CD68

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD68 signature in CD45 compartment 

CD45_CD68_idx = rep(0,CD68_nGene)
for(i in 1:CD68_nGene){
  idx = which(colnames(CD45_df)==CD68_final_genes[i])
  if(length(idx)>0){
    CD45_CD68_idx[i] = idx
  }
}

CD45_CD68_final_genes = CD68_final_genes[which(CD45_CD68_idx>0)]
CD45_CD68_final_coeffs = CD68_final_coeffs[which(CD45_CD68_idx>0)]
CD45_CD68_nGene = length(CD45_CD68_final_genes)

## AUC analysis

CD45_CD68_AUC = CD45_df[,CD45_CD68_idx[CD45_CD68_idx>0]]
CD45_CD68_AUC$response = CD45_data$response

CD45_CD68_AUC$lasso_scores = matrix(0,dim(CD45_CD68_AUC)[1],1)
for(i in 1:dim(CD45_CD68_AUC)[1]){
  CD45_CD68_AUC$lasso_scores[i] = sum(as.matrix(CD45_CD68_AUC[i,1:CD45_CD68_nGene])*CD45_CD68_final_coeffs[1:CD45_CD68_nGene])
}

pROC_CD45_CD68 <- roc(CD45_CD68_AUC$response,CD45_CD68_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_CD45_CD68

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD68 signature in pseudoBulk 

pBulk_CD68_idx = rep(0,CD68_nGene)
for(i in 1:CD68_nGene){
  idx = which(colnames(pBulk_df)==CD68_final_genes[i])
  if(length(idx)>0){
    pBulk_CD68_idx[i] = idx
  }
}

pBulk_CD68_final_genes = CD68_final_genes[which(pBulk_CD68_idx>0)]
pBulk_CD68_final_coeffs = CD68_final_coeffs[which(pBulk_CD68_idx>0)]
pBulk_CD68_nGene = length(pBulk_CD68_final_genes)

## AUC analysis

pBulk_CD68_AUC = pBulk_df[,pBulk_CD68_idx[pBulk_CD68_idx>0]]
pBulk_CD68_AUC$response = pBulk_data$response

pBulk_CD68_AUC$lasso_scores = matrix(0,dim(pBulk_CD68_AUC)[1],1)
for(i in 1:dim(pBulk_CD68_AUC)[1]){
  pBulk_CD68_AUC$lasso_scores[i] = sum(as.matrix(pBulk_CD68_AUC[i,1:pBulk_CD68_nGene])*pBulk_CD68_final_coeffs[1:pBulk_CD68_nGene])
}

pROC_pBulk_CD68 <- roc(pBulk_CD68_AUC$response,pBulk_CD68_AUC$lasso_scores,
                      smoothed = TRUE,
                      # arguments for ci
                      ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                      # arguments for plot
                      plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                      print.auc=TRUE, show.thres=TRUE)

pROC_pBulk_CD68

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# CD68 signature in pseudoBulk 

pStroma_CD68_idx = rep(0,CD68_nGene)
for(i in 1:CD68_nGene){
  idx = which(colnames(pStroma_df)==CD68_final_genes[i])
  if(length(idx)>0){
    pStroma_CD68_idx[i] = idx
  }
}

pStroma_CD68_final_genes = CD68_final_genes[which(pStroma_CD68_idx>0)]
pStroma_CD68_final_coeffs = CD68_final_coeffs[which(pStroma_CD68_idx>0)]
pStroma_CD68_nGene = length(pStroma_CD68_final_genes)

## AUC analysis

pStroma_CD68_AUC = pStroma_df[,pStroma_CD68_idx[pStroma_CD68_idx>0]]
pStroma_CD68_AUC$response = pStroma_data$response

pStroma_CD68_AUC$lasso_scores = matrix(0,dim(pStroma_CD68_AUC)[1],1)
for(i in 1:dim(pStroma_CD68_AUC)[1]){
  pStroma_CD68_AUC$lasso_scores[i] = sum(as.matrix(pStroma_CD68_AUC[i,1:pStroma_CD68_nGene])*pStroma_CD68_final_coeffs[1:pStroma_CD68_nGene])
}

pROC_pStroma_CD68 <- roc(pStroma_CD68_AUC$response,pStroma_CD68_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_pStroma_CD68

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoBulk signature in pseudoBulk compartment 

pBulk_nGene = dim(pBulk_Sig_sm)[1] - 1 # remove first row
pBulk_final_genes = pBulk_Sig_sm$X[2:(pBulk_nGene+1)]
pBulk_final_coeffs = pBulk_Sig_sm$s0[2:(pBulk_nGene+1)]

pBulk_pBulk_idx = rep(0,pBulk_nGene)
for(i in 1:pBulk_nGene){
  idx = which(colnames(pBulk_df)==pBulk_final_genes[i])
  if(length(idx)>0){
    pBulk_pBulk_idx[i] = idx
  }
}

pBulk_pBulk_final_genes = pBulk_final_genes[which(pBulk_pBulk_idx>0)]
pBulk_pBulk_final_coeffs = pBulk_final_coeffs[which(pBulk_pBulk_idx>0)]
pBulk_pBulk_nGene = length(pBulk_pBulk_final_genes)

## AUC analysis

pBulk_pBulk_AUC = pBulk_df[,pBulk_pBulk_idx[pBulk_pBulk_idx>0]]
pBulk_pBulk_AUC$response = pBulk_data$response

pBulk_pBulk_AUC$lasso_scores = matrix(0,dim(pBulk_pBulk_AUC)[1],1)
for(i in 1:dim(pBulk_pBulk_AUC)[1]){
  pBulk_pBulk_AUC$lasso_scores[i] = sum(as.matrix(pBulk_pBulk_AUC[i,1:pBulk_pBulk_nGene])*pBulk_pBulk_final_coeffs[1:pBulk_pBulk_nGene])
}

pROC_pBulk_pBulk <- roc(pBulk_pBulk_AUC$response,pBulk_pBulk_AUC$lasso_scores,
                        smoothed = TRUE,
                        # arguments for ci
                        ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                        # arguments for plot
                        plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                        print.auc=TRUE, show.thres=TRUE)

pROC_pBulk_pBulk

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoBulk signature in S100B compartment 

pBulk_nGene = dim(pBulk_Sig_sm)[1] - 1 # remove first row
pBulk_final_genes = pBulk_Sig_sm$X[2:(pBulk_nGene+1)]
pBulk_final_coeffs = pBulk_Sig_sm$s0[2:(pBulk_nGene+1)]

S100B_pBulk_idx = rep(0,pBulk_nGene)
for(i in 1:pBulk_nGene){
  idx = which(colnames(S100B_df)==pBulk_final_genes[i])
  if(length(idx)>0){
    S100B_pBulk_idx[i] = idx
  }
}

S100B_pBulk_final_genes = pBulk_final_genes[which(S100B_pBulk_idx>0)]
S100B_pBulk_final_coeffs = pBulk_final_coeffs[which(S100B_pBulk_idx>0)]
S100B_pBulk_nGene = length(S100B_pBulk_final_genes)

## AUC analysis

S100B_pBulk_AUC = S100B_df[,S100B_pBulk_idx[S100B_pBulk_idx>0]]
S100B_pBulk_AUC$response = S100B_data$response

S100B_pBulk_AUC$lasso_scores = matrix(0,dim(S100B_pBulk_AUC)[1],1)
for(i in 1:dim(S100B_pBulk_AUC)[1]){
  S100B_pBulk_AUC$lasso_scores[i] = sum(as.matrix(S100B_pBulk_AUC[i,1:S100B_pBulk_nGene])*S100B_pBulk_final_coeffs[1:S100B_pBulk_nGene])
}

pROC_S100B_pBulk <- roc(S100B_pBulk_AUC$response,S100B_pBulk_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_S100B_pBulk

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoBulk signature in CD45 compartment 

CD45_pBulk_idx = rep(0,pBulk_nGene)
for(i in 1:pBulk_nGene){
  idx = which(colnames(CD45_df)==pBulk_final_genes[i])
  if(length(idx)>0){
    CD45_pBulk_idx[i] = idx
  }
}

CD45_pBulk_final_genes = pBulk_final_genes[which(CD45_pBulk_idx>0)]
CD45_pBulk_final_coeffs = pBulk_final_coeffs[which(CD45_pBulk_idx>0)]
CD45_pBulk_nGene = length(CD45_pBulk_final_genes)

## AUC analysis

CD45_pBulk_AUC = CD45_df[,CD45_pBulk_idx[CD45_pBulk_idx>0]]
CD45_pBulk_AUC$response = CD45_data$response

CD45_pBulk_AUC$lasso_scores = matrix(0,dim(CD45_pBulk_AUC)[1],1)
for(i in 1:dim(CD45_pBulk_AUC)[1]){
  CD45_pBulk_AUC$lasso_scores[i] = sum(as.matrix(CD45_pBulk_AUC[i,1:CD45_pBulk_nGene])*CD45_pBulk_final_coeffs[1:CD45_pBulk_nGene])
}

pROC_CD45_pBulk <- roc(CD45_pBulk_AUC$response,CD45_pBulk_AUC$lasso_scores,
                        smoothed = TRUE,
                        # arguments for ci
                        ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                        # arguments for plot
                        plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                        print.auc=TRUE, show.thres=TRUE)

pROC_CD45_pBulk

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoBulk signature in CD65 compartment 

CD68_pBulk_idx = rep(0,pBulk_nGene)
for(i in 1:pBulk_nGene){
  idx = which(colnames(CD68_df)==pBulk_final_genes[i])
  if(length(idx)>0){
    CD68_pBulk_idx[i] = idx
  }
}

CD68_pBulk_final_genes = pBulk_final_genes[which(CD68_pBulk_idx>0)]
CD68_pBulk_final_coeffs = pBulk_final_coeffs[which(CD68_pBulk_idx>0)]
CD68_pBulk_nGene = length(CD68_pBulk_final_genes)

## AUC analysis

CD68_pBulk_AUC = CD68_df[,CD68_pBulk_idx[CD68_pBulk_idx>0]]
CD68_pBulk_AUC$response = CD68_data$response

CD68_pBulk_AUC$lasso_scores = matrix(0,dim(CD68_pBulk_AUC)[1],1)
for(i in 1:dim(CD68_pBulk_AUC)[1]){
  CD68_pBulk_AUC$lasso_scores[i] = sum(as.matrix(CD68_pBulk_AUC[i,1:CD68_pBulk_nGene])*CD68_pBulk_final_coeffs[1:CD68_pBulk_nGene])
}

pROC_CD68_pBulk <- roc(CD68_pBulk_AUC$response,CD68_pBulk_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_CD68_pBulk

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoBulk signature in pseudoStroma 

pStroma_pBulk_idx = rep(0,pBulk_nGene)
for(i in 1:pBulk_nGene){
  idx = which(colnames(pStroma_df)==pBulk_final_genes[i])
  if(length(idx)>0){
    pStroma_pBulk_idx[i] = idx
  }
}

pStroma_pBulk_final_genes = pBulk_final_genes[which(pStroma_pBulk_idx>0)]
pStroma_pBulk_final_coeffs = pBulk_final_coeffs[which(pStroma_pBulk_idx>0)]
pStroma_pBulk_nGene = length(pStroma_pBulk_final_genes)

## AUC analysis

pStroma_pBulk_AUC = pStroma_df[,pStroma_pBulk_idx[pStroma_pBulk_idx>0]]
pStroma_pBulk_AUC$response = pStroma_data$response

pStroma_pBulk_AUC$lasso_scores = matrix(0,dim(pStroma_pBulk_AUC)[1],1)
for(i in 1:dim(pStroma_pBulk_AUC)[1]){
  pStroma_pBulk_AUC$lasso_scores[i] = sum(as.matrix(pStroma_pBulk_AUC[i,1:pStroma_pBulk_nGene])*pStroma_pBulk_final_coeffs[1:pStroma_pBulk_nGene])
}

pROC_pStroma_pBulk <- roc(pStroma_pBulk_AUC$response,pStroma_pBulk_AUC$lasso_scores,
                       smoothed = TRUE,
                       # arguments for ci
                       ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                       # arguments for plot
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE)

pROC_pStroma_pBulk

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoStroma signature in pseudoStroma compartment 

pStroma_nGene = dim(pStroma_Sig_sm)[1] - 1 # remove first row
pStroma_final_genes = pStroma_Sig_sm$X[2:(pStroma_nGene+1)]
pStroma_final_coeffs = pStroma_Sig_sm$s0[2:(pStroma_nGene+1)]

pStroma_pStroma_idx = rep(0,pStroma_nGene)
for(i in 1:pStroma_nGene){
  idx = which(colnames(pStroma_df)==pStroma_final_genes[i])
  if(length(idx)>0){
    pStroma_pStroma_idx[i] = idx
  }
}

pStroma_pStroma_final_genes = pStroma_final_genes[which(pStroma_pStroma_idx>0)]
pStroma_pStroma_final_coeffs = pStroma_final_coeffs[which(pStroma_pStroma_idx>0)]
pStroma_pStroma_nGene = length(pStroma_pStroma_final_genes)

## AUC analysis

pStroma_pStroma_AUC = pStroma_df[,pStroma_pStroma_idx[pStroma_pStroma_idx>0]]
pStroma_pStroma_AUC$response = pStroma_data$response

pStroma_pStroma_AUC$lasso_scores = matrix(0,dim(pStroma_pStroma_AUC)[1],1)
for(i in 1:dim(pStroma_pStroma_AUC)[1]){
  pStroma_pStroma_AUC$lasso_scores[i] = sum(as.matrix(pStroma_pStroma_AUC[i,1:pStroma_pStroma_nGene])*pStroma_pStroma_final_coeffs[1:pStroma_pStroma_nGene])
}

pROC_pStroma_pStroma <- roc(pStroma_pStroma_AUC$response,pStroma_pStroma_AUC$lasso_scores,
                          smoothed = TRUE,
                          # arguments for ci
                          ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                          # arguments for plot
                          plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                          print.auc=TRUE, show.thres=TRUE)

pROC_pStroma_pStroma

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoStroma signature in S100B compartment 

pStroma_nGene = dim(pStroma_Sig_sm)[1] - 1 # remove first row
pStroma_final_genes = pStroma_Sig_sm$X[2:(pStroma_nGene+1)]
pStroma_final_coeffs = pStroma_Sig_sm$s0[2:(pStroma_nGene+1)]

S100B_pStroma_idx = rep(0,pStroma_nGene)
for(i in 1:pStroma_nGene){
  idx = which(colnames(S100B_df)==pStroma_final_genes[i])
  if(length(idx)>0){
    S100B_pStroma_idx[i] = idx
  }
}

S100B_pStroma_final_genes = pStroma_final_genes[which(S100B_pStroma_idx>0)]
S100B_pStroma_final_coeffs = pStroma_final_coeffs[which(S100B_pStroma_idx>0)]
S100B_pStroma_nGene = length(S100B_pStroma_final_genes)

## AUC analysis

S100B_pStroma_AUC = S100B_df[,S100B_pStroma_idx[S100B_pStroma_idx>0]]
S100B_pStroma_AUC$response = S100B_data$response

S100B_pStroma_AUC$lasso_scores = matrix(0,dim(S100B_pStroma_AUC)[1],1)
for(i in 1:dim(S100B_pStroma_AUC)[1]){
  S100B_pStroma_AUC$lasso_scores[i] = sum(as.matrix(S100B_pStroma_AUC[i,1:S100B_pStroma_nGene])*S100B_pStroma_final_coeffs[1:S100B_pStroma_nGene])
}

pROC_S100B_pStroma <- roc(S100B_pStroma_AUC$response,S100B_pStroma_AUC$lasso_scores,
                        smoothed = TRUE,
                        # arguments for ci
                        ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                        # arguments for plot
                        plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                        print.auc=TRUE, show.thres=TRUE)

pROC_S100B_pStroma

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoStroma signature in CD45 compartment 

CD45_pStroma_idx = rep(0,pStroma_nGene)
for(i in 1:pStroma_nGene){
  idx = which(colnames(CD45_df)==pStroma_final_genes[i])
  if(length(idx)>0){
    CD45_pStroma_idx[i] = idx
  }
}

CD45_pStroma_final_genes = pStroma_final_genes[which(CD45_pStroma_idx>0)]
CD45_pStroma_final_coeffs = pStroma_final_coeffs[which(CD45_pStroma_idx>0)]
CD45_pStroma_nGene = length(CD45_pStroma_final_genes)

## AUC analysis

CD45_pStroma_AUC = CD45_df[,CD45_pStroma_idx[CD45_pStroma_idx>0]]
CD45_pStroma_AUC$response = CD45_data$response

CD45_pStroma_AUC$lasso_scores = matrix(0,dim(CD45_pStroma_AUC)[1],1)
for(i in 1:dim(CD45_pStroma_AUC)[1]){
  CD45_pStroma_AUC$lasso_scores[i] = sum(as.matrix(CD45_pStroma_AUC[i,1:CD45_pStroma_nGene])*CD45_pStroma_final_coeffs[1:CD45_pStroma_nGene])
}

pROC_CD45_pStroma <- roc(CD45_pStroma_AUC$response,CD45_pStroma_AUC$lasso_scores,
                          smoothed = TRUE,
                          # arguments for ci
                          ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                          # arguments for plot
                          plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                          print.auc=TRUE, show.thres=TRUE)

pROC_CD45_pStroma

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoStroma signature in CD68 compartment 

CD68_pStroma_idx = rep(0,pStroma_nGene)
for(i in 1:pStroma_nGene){
  idx = which(colnames(CD68_df)==pStroma_final_genes[i])
  if(length(idx)>0){
    CD68_pStroma_idx[i] = idx
  }
}

CD68_pStroma_final_genes = pStroma_final_genes[which(CD68_pStroma_idx>0)]
CD68_pStroma_final_coeffs = pStroma_final_coeffs[which(CD68_pStroma_idx>0)]
CD68_pStroma_nGene = length(CD68_pStroma_final_genes)

## AUC analysis

CD68_pStroma_AUC = CD68_df[,CD68_pStroma_idx[CD68_pStroma_idx>0]]
CD68_pStroma_AUC$response = CD68_data$response

CD68_pStroma_AUC$lasso_scores = matrix(0,dim(CD68_pStroma_AUC)[1],1)
for(i in 1:dim(CD68_pStroma_AUC)[1]){
  CD68_pStroma_AUC$lasso_scores[i] = sum(as.matrix(CD68_pStroma_AUC[i,1:CD68_pStroma_nGene])*CD68_pStroma_final_coeffs[1:CD68_pStroma_nGene])
}

pROC_CD68_pStroma <- roc(CD68_pStroma_AUC$response,CD68_pStroma_AUC$lasso_scores,
                         smoothed = TRUE,
                         # arguments for ci
                         ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                         # arguments for plot
                         plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                         print.auc=TRUE, show.thres=TRUE)

pROC_CD68_pStroma

###############################################################
###############################################################

# generate lasso scores for each signatures in each dataset
# pseudoStroma signature in pseudoBulk 

pBulk_pStroma_idx = rep(0,pStroma_nGene)
for(i in 1:pStroma_nGene){
  idx = which(colnames(pBulk_df)==pStroma_final_genes[i])
  if(length(idx)>0){
    pBulk_pStroma_idx[i] = idx
  }
}

pBulk_pStroma_final_genes = pStroma_final_genes[which(pBulk_pStroma_idx>0)]
pBulk_pStroma_final_coeffs = pStroma_final_coeffs[which(pBulk_pStroma_idx>0)]
pBulk_pStroma_nGene = length(pBulk_pStroma_final_genes)

## AUC analysis

pBulk_pStroma_AUC = pBulk_df[,pBulk_pStroma_idx[pBulk_pStroma_idx>0]]
pBulk_pStroma_AUC$response = pBulk_data$response

pBulk_pStroma_AUC$lasso_scores = matrix(0,dim(pBulk_pStroma_AUC)[1],1)
for(i in 1:dim(pBulk_pStroma_AUC)[1]){
  pBulk_pStroma_AUC$lasso_scores[i] = sum(as.matrix(pBulk_pStroma_AUC[i,1:pBulk_pStroma_nGene])*pBulk_pStroma_final_coeffs[1:pBulk_pStroma_nGene])
}

pROC_pBulk_pStroma <- roc(pBulk_pStroma_AUC$response,pBulk_pStroma_AUC$lasso_scores,
                         smoothed = TRUE,
                         # arguments for ci
                         ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                         # arguments for plot
                         plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                         print.auc=TRUE, show.thres=TRUE)

pROC_pBulk_pStroma

###################################################################################

myList <- list(S100B_S100B_AUC = S100B_S100B_AUC, 
               CD45_S100B_AUC = CD45_S100B_AUC, 
               CD68_S100B_AUC = CD68_S100B_AUC,
               pBulk_S100B_AUC = pBulk_S100B_AUC,
               pStroma_S100B_AUC = pStroma_S100B_AUC,
               CD45_CD45_AUC = CD45_CD45_AUC,
               CD68_CD45_AUC = CD68_CD45_AUC,
               S100B_CD45_AUC = S100B_CD45_AUC,
               pBulk_CD45_AUC = pBulk_CD45_AUC,
               pStroma_CD45_AUC = pStroma_CD45_AUC,
               CD68_CD68_AUC = CD68_CD68_AUC,
               CD45_CD68_AUC = CD45_CD68_AUC,
               S100B_CD68_AUC = S100B_CD68_AUC,
               pBulk_CD68_AUC = pBulk_CD68_AUC,
               pStroma_CD68_AUC = pStroma_CD68_AUC,
               pBulk_pBulk_AUC = pBulk_pBulk_AUC,
               CD45_pBulk_AUC = CD45_pBulk_AUC,
               CD68_pBulk_AUC = CD68_pBulk_AUC,
               S100B_pBulk_AUC = S100B_pBulk_AUC,
               pStroma_pBulk_AUC = pStroma_pBulk_AUC,
               pStroma_pStroma_AUC = pStroma_pStroma_AUC,
               CD45_pStroma_AUC = CD45_pStroma_AUC,
               CD68_pStroma_AUC = CD68_pStroma_AUC,
               S100B_pStroma_AUC = S100B_pStroma_AUC,
               pBulk_pStroma_AUC = pBulk_pStroma_AUC)
for(i in names(myList)){
  write.csv(myList[[i]], paste0(i,".csv"))
}

###################################################################################

#####

CD45_sig = read.csv("Diff_signatures_in_CD45_comp_sm.csv")


Diff_Sig_CD45_comp <- roc(response ~ CD45_in_CD45 + CD68_in_CD45 + S100B_in_CD45 + pBulk_in_CD45 + pStroma_in_CD45,
                          data = CD45_sig)

Diff_Sig_CD45_comp2 <- ggroc(Diff_Sig_CD45_comp, size = 1.5)
Diff_Sig_CD45_comp2 + theme_classic()+
  scale_colour_manual(values = c("red", "blue",  "green", "darkorchid1", "turquoise"))
#scale_colour_brewer(palette="Set1")

#####

CD68_sig = read.csv("Diff_signatures_in_CD68_comp_sm.csv")


Diff_Sig_CD68_comp <- roc(response ~ CD68_in_CD68 + CD45_in_CD68 + S100B_in_CD68 + pBulk_in_CD68 + pStroma_in_CD68,
                          data = CD68_sig)

Diff_Sig_CD68_comp2 <- ggroc(Diff_Sig_CD68_comp, size = 1.5)
Diff_Sig_CD68_comp2 + theme_classic()+
  scale_colour_manual(values = c("red", "blue",  "green", "darkorchid1", "turquoise"))
#scale_colour_brewer(palette="Set1")

#####

S100B_sig = read.csv("Diff_signatures_in_S100B_comp_sm.csv")


Diff_Sig_S100B_comp <- roc(response ~ S100B_in_S100B + CD45_in_S100B + CD68_in_S100B + pBulk_in_S100B + pStroma_in_S100B,
                          data = S100B_sig)

Diff_Sig_S100B_comp2 <- ggroc(Diff_Sig_S100B_comp, size = 1.5)
Diff_Sig_S100B_comp2 + theme_classic()+
  scale_colour_manual(values = c("red", "blue",  "green", "darkorchid1", "turquoise"))
#scale_colour_brewer(palette="Set1")

#####

pBulk_sig = read.csv("Diff_signatures_in_pBulk_comp_sm.csv")


Diff_Sig_pBulk_comp <- roc(response ~ pBulk_in_pBulk + CD45_in_pBulk + CD68_in_pBulk + S100B_in_pBulk + pStroma_in_pBulk,
                           data = pBulk_sig)

Diff_Sig_pBulk_comp2 <- ggroc(Diff_Sig_pBulk_comp, size = 1.5)
Diff_Sig_pBulk_comp2 + theme_classic()+
  scale_colour_manual(values = c("red", "blue",  "green", "darkorchid1", "turquoise"))
#scale_colour_brewer(palette="Set1")

#####

pStroma_sig = read.csv("Diff_signatures_in_pStroma_comp_sm.csv")


Diff_Sig_pStroma_comp <- roc(response ~ pStroma_in_pStroma + CD45_in_pStroma + CD68_in_pStroma + S100B_in_pStroma + pBulk_in_pStroma,
                           data = pStroma_sig)

Diff_Sig_pStroma_comp2 <- ggroc(Diff_Sig_pStroma_comp, size = 1.5)
Diff_Sig_pStroma_comp2 + theme_classic()+
  scale_colour_manual(values = c("red", "blue",  "green", "darkorchid1", "turquoise"))
#scale_colour_brewer(palette="Set1")


