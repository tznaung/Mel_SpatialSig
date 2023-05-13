library(pheatmap)
library(RColorBrewer)

CD68 = read.csv("filter_5%_LOQ_CD68.csv", header = T)
CD45 = read.csv("filter_5%_LOQ_CD45.csv", header = T)
S100B = read.csv("filter_5%_LOQ_S100B.csv", header = T)

#get the average counts of both blocks and combine CD45 and CD68 to get stroma compartment
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

for(i in 1:dim(CD45_av)[1]){  
  CD45_av[i,] = 1000000 * CD45_av[i,] / sum(CD45_av[i,])
}
for(i in 1:dim(CD68_av)[1]){  
  CD68_av[i,] = 1000000 * CD68_av[i,] / sum(CD68_av[i,])
}
for(i in 1:dim(S100B_av)[1]){  
  S100B_av[i,] = 1000000 * S100B_av[i,] / sum(S100B_av[i,])
}

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
for(i in 1:dim(pseudo_bulk)[1]){  
  pseudo_bulk[i,] = 1000000 * pseudo_bulk[i,] / sum(pseudo_bulk[i,])
}
nSubj = length(subj_names)

pseudo_bulk_norm = log2(pseudo_bulk + 1)

pseudo_bulk_zscore = scale(pseudo_bulk_norm)

pseudo_bulk_zscore_df = as.data.frame(pseudo_bulk_zscore)
colnames(pseudo_bulk_zscore_df) = allnames
rownames(pseudo_bulk_zscore_df) = subj_names

## diff signatures test

published_sig = read.csv("Signatures.csv")
published_sig = published_sig[,1:7]
nSig = dim(published_sig)[2]
sig_scores = matrix(0,nSubj,nSig)
for (i in 1:nSig){
  idxs = which(allnames %in% published_sig[,i])
  for (j in 1:nSubj){
    sig_scores[,i] = rowSums(pseudo_bulk_zscore[,idxs])
  }
}

sig_scores_df = as.data.frame(sig_scores)
colnames(sig_scores_df) = colnames(published_sig)
rownames(sig_scores_df) = subj_names
write.csv(sig_scores_df, "sig_scores_df.csv")

Response = read.csv("Response.csv")

resp_vec = rep(0,nSubj)
for (i in 1:nSubj){
  idx = which(Response$ROI_ID_2 == subj_names[i])
  idx = idx[1]
  if (Response$response[idx]=="yes"){
    resp_vec[i] = 1
  }
}
resp = which(resp_vec==1)
nonresp = which(resp_vec==0)

# do t-tests
pvals = matrix(0,1,nSig)
mean_diffs = matrix(0,1,nSig)
for (i in 1:nSig){
  resp_vec = sig_scores[resp,i]
  nonresp_vec = sig_scores[nonresp,i]
  t = t.test(resp_vec,nonresp_vec, alternative = "greater")
  pvals[1,i] = t$p.value
  mean_diffs[1,i] = mean(resp_vec) - mean(nonresp_vec)  
}

colnames(pvals) = colnames(published_sig)
colnames(mean_diffs) = colnames(published_sig)
pvals
mean_diffs

Pvals_MeanDiff = rbind(pvals, mean_diffs)
Pvals_MeanDiff
write.csv(Pvals_MeanDiff, "Pvals_MeanDiff.csv")

# heatmap

df = read.csv("sig_scores_df_2.csv", row.names = 1)
#  sapply(rownames(train_data4),function(x) 
#  strsplit(as.character(x),split = "\\\\")[[1]][1])
df2 <- df[order(df$response),]
df3 = data.frame("Response" = sort(df$response))
rownames(df3) = rownames(df2)
#rownames(train_data4) = train_data5$Response # name matching
pheatmap(scale(df2[,1:7]),annotation_row = df3, main = "Signature genes bulk",
         cluster_rows = T, color = colorRampPalette(rev(brewer.pal(n = 11, name ="PiYG")))(15),
         # cutree_cols = 2,
         #cutree_rows = 2, 
         fontsize = 8, 
         clustering_method = "complete",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
pheatmap(scale(df2[,1:7]),annotation_row = df3, main = "Signature genes for bulk",
         cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(7),
         cutree_cols = 2,
         #cutree_rows = 3, 
         fontsize = 8, angle_col = "315")
?pheatmap
