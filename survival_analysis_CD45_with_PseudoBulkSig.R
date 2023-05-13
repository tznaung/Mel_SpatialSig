library("survminer")
library("ggplot2")
require("survival")
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(edgeR)
library(pROC)

CD45_df = read.csv("CD45_av.csv", header = T)
CD45_df$SpotID <- CD45_df$X
Res_CD45 = read.csv("CD45_response_table.csv", header = T)
Res_CD45$SpotID <- Res_CD45$ROI_ID_2
CD45_2 = merge(CD45_df, Res_CD45, by = "SpotID", all = FALSE)
CD45_2$group <- CD45_2$response # define the group you want to compare
data = CD45_2[,3:3738]
prop_test = 0.20
idx0 = which(CD45_2$group=="no")
idx1 = which(CD45_2$group=="yes")

final_coeffs = read.csv("final_coeffs_pseudoBulk.csv")
# final_coeffs_2 = final_coeffs[2:14,]
final_genes = final_coeffs$X[2:dim(final_coeffs)[1]]
nGene = length(final_genes)

train_data =  CD45_2
train_data2 = train_data[,3:3738]
train_data3 = train_data[,which(colnames(train_data) %in% final_genes)]
final_genes_2 = colnames(train_data3)
nGene_2 = length(final_genes_2)
# final_coeffs_2 = final_coeffs[c(3:10,12:14),]
train_data3$response = train_data$response
train_data3$time = train_data$OS_FROM_START_OF_ITX
train_data3$VITAL = train_data$VITAL

train_data3$lasso_scores = matrix(0,dim(train_data3)[1],1)
for(i in 1:dim(train_data3)[1]){
  train_data3$lasso_scores[i] = sum(as.matrix(train_data3[i,1:nGene_2])*(final_coeffs$s0[2:dim(final_coeffs)[1]]))
} 

## heatmap

#pdf(file = paste("output_split_CD45_final.pdf", sep=""))

# Training set 
train_data4 = train_data3[,1:length(final_genes_2)]
train_data5 <- log2(train_data4)
rownames(train_data5) = train_data$X
#  sapply(rownames(train_data4),function(x) 
#  strsplit(as.character(x),split = "\\\\")[[1]][1])
train_data5 <- train_data5[order(train_data$response),]
train_data6 = data.frame("Response" = sort(train_data$response))
rownames(train_data6) = rownames(train_data5)
#rownames(train_data4) = train_data5$Response # name matching
pheatmap(train_data5,annotation_row = train_data6, main = "Signature genes CD45 Comp",
         cluster_rows = T, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
         # cutree_cols = 2,
         cutree_rows = 4, 
         fontsize = 8)
# pheatmap(train_data4,annotation_row = train_data5, main = "Signature genes for CD45 comp",
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
           title = "CD45 compartment w/ PseudoBulk signatures",
           pval = F, pval.method = TRUE,    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "Lasso",               # Change legend titles
           legend.labs = c("Low", "High"),  # Change legend labels
           palette = "lancet",                    # Use JCO journal color palette
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
                    x = 1000, y = 0.25, # x and y coordinates of the text
                    label = "HR = 0.54 (0.15-2) \n p = 0.359", size = 4)
A

# CD45 Compartment with PseudoBulk Signatures
