library("survminer")
library("ggplot2")
require("survival")
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(edgeR)
library(pROC)


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

final_coeffs = read.csv("final_coeffs_pseudoBulk.csv")
final_genes = final_coeffs$X[2:dim(final_coeffs)[1]]
nGene = length(final_genes)

train_data = Bulk_df
train_data2 = train_data[,3:1631]
train_data3 = train_data[,which(colnames(train_data) %in% final_genes)]
train_data3$response = train_data$response
train_data3$time = train_data$OS_FROM_START_OF_ITX
train_data3$VITAL = train_data$VITAL

train_data3$lasso_scores = matrix(0,dim(train_data3)[1],1)
for(i in 1:dim(train_data3)[1]){
  train_data3$lasso_scores[i] = sum(as.matrix(train_data3[i,1:nGene])*(final_coeffs$s0[2:dim(final_coeffs)[1]]))
} 

## heatmap

#pdf(file = paste("output_split_pseudoBulk_final.pdf", sep=""))

# Training set 
train_data4 = train_data3[,1:length(final_genes)]
train_data4 <- log2(train_data4)
rownames(train_data4) = train_data$X
#  sapply(rownames(train_data4),function(x) 
#  strsplit(as.character(x),split = "\\\\")[[1]][1])
train_data4 <- train_data4[order(train_data$response),]
train_data5 = data.frame("Response" = sort(train_data$response))
rownames(train_data5) = rownames(train_data4)
#rownames(train_data4) = train_data5$Response # name matching
pheatmap(train_data4,annotation_row = train_data5, main = "Signature genes pseudoBulk Comp",
         cluster_rows = T, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
         # cutree_cols = 2,
         cutree_rows = 4, 
         fontsize = 8)
# pheatmap(train_data4,annotation_row = train_data5, main = "Signature genes for pseudoBulk comp",
#          cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
#          cutree_cols = 2,
#          #cutree_rows = 3, 
#          fontsize = 8)

#dev.off()

#Finding the optimal cutpoin
library(pROC)
library(cutpointr)

Beta <- cutpointr(train_data3, lasso_scores, VITAL, 
                  method = maximize_metric, metric = sum_sens_spec)
summary(Beta)
plot(Beta)

#train_data3$binary_score = (train_data3$lasso_scores>=Beta$optimal_cutpoint)
#train_data3$binary_score = (train_data3$lasso_scores>=median(train_data3$lasso_scores))
train_data3$binary_score = (train_data3$lasso_scores>=quantile(train_data3$lasso_scores,probs=0.66))
train_data3$binary_resp = (train_data3$response == "yes")

# survival curve

fit <- survfit(Surv(time, VITAL) ~ binary_score, data = train_data3)
fit$n.risk
# Plot informative survival curves
A = ggsurvplot(fit, data = train_data3,
           title = "Pseudo-Bulk",
           pval = F, pval.method = TRUE,    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "Lasso",               # Change legend titles
           legend.labs = c("Low", "High"),  # Change legend labels
           palette = "jco",                    # Use JCO journal color palette
           risk.table = TRUE,                  # Add No at risk table
           cumevents = FALSE,                  # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE               # Hide tables y axis text
)
# find the hazard ratio with this function.
cox <- coxph(Surv(time, VITAL) ~ binary_score, data = train_data3)
cox
ggforest(cox)
A$plot <- A$plot+ 
  ggplot2::annotate("text", 
                    x = 1000, y = 0.6, # x and y coordinates of the text
                    label = "HR = 0.35 (0.12-1) \n p = 0.036*", size = 4)
A

