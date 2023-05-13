library("survminer")
library("ggplot2")
require("survival")
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(edgeR)
library(pROC)

CD68_df = read.csv("CD68_av.csv", header = T)
CD68_df$SpotID <- CD68_df$X
Res_CD68 = read.csv("CD68_response_table.csv", header = T)
Res_CD68$SpotID <- Res_CD68$ROI_ID_2
CD68_2 = merge(CD68_df, Res_CD68, by = "SpotID", all = FALSE)
CD68_2$group <- CD68_2$response # define the group you want to compare
data = CD68_2[,3:1713]
prop_test = 0.20
idx0 = which(CD68_2$group=="no")
idx1 = which(CD68_2$group=="yes")

final_coeffs = read.csv("final_coeffs_CD68.csv")
final_genes = final_coeffs$X[2:dim(final_coeffs)[1]]
nGene = length(final_genes)

train_data =  CD68_2
train_data2 = train_data[,3:1713]
train_data3 = train_data[,which(colnames(train_data) %in% final_genes)]
train_data3$response = train_data$response
train_data3$time = train_data$OS_FROM_START_OF_ITX
train_data3$VITAL = train_data$VITAL

train_data3$lasso_scores = matrix(0,dim(train_data3)[1],1)
for(i in 1:dim(train_data3)[1]){
  train_data3$lasso_scores[i] = sum(as.matrix(train_data3[i,1:nGene])*(final_coeffs$s0[2:dim(final_coeffs)[1]]))
} 

## heatmap

#pdf(file = paste("output_split_CD68_final.pdf", sep=""))

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
pheatmap(train_data4,annotation_row = train_data5, main = "Signature genes CD68 Comp",
         cluster_rows = T, color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn")))(15),
         # cutree_cols = 2,
         cutree_rows = 4, 
         fontsize = 8)
# pheatmap(train_data4,annotation_row = train_data5, main = "Signature genes for CD68 comp",
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
           title = "CD68 compartment",
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
                    x = 1000, y = 0.6, # x and y coordinates of the text
                    label = "HR = 0.29 (0.084-1) \n p = 0.02*", size = 4)
A
