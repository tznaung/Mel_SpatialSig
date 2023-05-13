S100B_av = read.csv("S100B_av.csv")
CD45_av = read.csv("CD45_av.csv")
CD68_av = read.csv("CD68_av.csv")
stroma_av = read.csv("pseudo_stroma_av.csv")
stroma_av_2 = stroma_av[, 1:1643]

mean_S100B = colMeans(S100B_av)
mean_CD45 = colMeans(CD45_av)
mean_CD68 = colMeans(CD68_av)
mean_stroma = colMeans(stroma_av_2)

nTopGenes = 640  

# get compartment genes
S100B_genes = names(sort(mean_S100B,decreasing = TRUE))[1:nTopGenes]
CD45_genes = names(sort(mean_CD45,decreasing = TRUE))[1:nTopGenes]
CD68_genes = names(sort(mean_CD68,decreasing = TRUE))[1:nTopGenes]

# get r/nr-sig genes
coeffs = read.csv("final_coeffs_pseudoBulk.csv")
resp_sig_genes = coeffs$X[2:dim(coeffs)[1]]
coeffs = read.csv("final_coeffs_pseudoStroma.csv")
resp_sig_genes = union(resp_sig_genes,coeffs$X[2:dim(coeffs)[1]])
coeffs = read.csv("final_coeffs_S100B.csv")
resp_sig_genes = union(resp_sig_genes,coeffs$X[2:dim(coeffs)[1]])
coeffs = read.csv("final_coeffs_CD45.csv")
resp_sig_genes = union(resp_sig_genes,coeffs$X[2:dim(coeffs)[1]])
coeffs = read.csv("final_coeffs_CD68.csv")
resp_sig_genes = union(resp_sig_genes,coeffs$X[2:dim(coeffs)[1]])
coeffs = read.csv("final_coeffs_pseudoStroma.csv")
resp_sig_genes = union(resp_sig_genes,coeffs$X[2:dim(coeffs)[1]])

# allgenes = union(resp_sig_genes,union(S100B_genes,union(CD45_genes,CD68_genes)))
nGene = length(allgenes)

# make compartment sig matrix
sig_matrix = matrix(0,3,nGene)
colnames(sig_matrix) = allgenes
rownames(sig_matrix) = c("S100B","CD45","CD68")
for (i in 1:nTopGenes){
  idx = which(allgenes==S100B_genes[i])
  idx2 = which(names(mean_S100B)==S100B_genes[i])
  sig_matrix[1,idx] = mean_S100B[idx2]
  idx = which(allgenes==CD45_genes[i])
  idx2 = which(names(mean_CD45)==CD45_genes[i])
  sig_matrix[2,idx] = mean_CD45[idx2]
  idx = which(allgenes==CD68_genes[i])
  idx2 = which(names(mean_CD68)==CD68_genes[i])
  sig_matrix[3,idx] = mean_CD68[idx2]  
}

for (i in 1:length(resp_sig_genes)){
  idx = which(allgenes==resp_sig_genes[i])
  idx2 = which(names(mean_S100B)==resp_sig_genes[i])
  if (length(idx2)>0){
    sig_matrix[1,idx] = mean_S100B[idx2]
  }
  idx2 = which(names(mean_CD45)==resp_sig_genes[i])
  if (length(idx2)>0){
    sig_matrix[2,idx] = mean_CD45[idx2]
  }
  idx2 = which(names(mean_CD68)==resp_sig_genes[i])
  if (length(idx2)>0){
    sig_matrix[3,idx] = mean_CD68[idx2]
  }  
}

# save it (PS remember to convert to .txt, and add GeneNames to first col)
write.csv(t(sig_matrix),"sig_matrix.csv")

# load external data
load("~/Desktop/DSP_Melanoma/0_NC_2NPC/12_validation_2/deconvolution_data.RData")
nPat = dim(Gide_data$TPM)[2]

#write out data, change - to . in names, and read back in
write.csv(Gide_data$TPM,"Gide_TPM.csv")
Gide_data_TPM = read.csv("Gide_TPM.csv",row.names = 1)

#make mixture matrix
gida_mat= as.data.frame(matrix(0,nGene,nPat))
for(i in 1:nGene){
  idx = which(allgenes[i] == rownames(Gide_data_TPM))
  if (length(idx)>0){
    gida_mat[i,] = Gide_data_TPM[idx,]
  }
}
rownames(gida_mat) = allgenes
colnames(gida_mat) = colnames(Gide_data_TPM)

#save it (PS remember to convert to .txt, and add GeneNames to first col)
write.csv(gida_mat,paste("gide_",as.numeric(nGene),".csv",sep=""))
dim(gida_mat)

#save gene list (PS reformat if nec, so it's only one column
write.table(allgenes,"gene_list.txt",sep = "\t")


