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

final_coeffs = read.csv("final_coeffs_CD68_8_genes.csv")
final_genes = final_coeffs$X[2:dim(final_coeffs)[1]]
nGene = length(final_genes)

train_data =  CD68_2
train_data2 = train_data[,3:1713]
idxs = rep(0,nGene)                 #### copy next 5 lines for reordering genes to match coefficients
for (i in 1:nGene){
  idxs[i] = which(colnames(train_data)==final_genes[i])
}
train_data3 = train_data[,idxs]
#train_data3 = train_data[,which(colnames(train_data) %in% final_genes)]
train_data3$response = train_data$response
train_data3$OS_FROM_START_OF_ITX = train_data$OS_FROM_START_OF_ITX
train_data3$VITAL = train_data$VITAL

train_data3$lasso_scores = matrix(0,dim(train_data3)[1],1)
for(i in 1:dim(train_data3)[1]){
  train_data3$lasso_scores[i] = sum(as.matrix(train_data3[i,1:nGene])*(final_coeffs$s0[2:dim(final_coeffs)[1]]))
} 



#Finding the optimal cutpoin
library(pROC)
library(cutpointr)

Beta <- cutpointr(train_data3, lasso_scores, VITAL, 
                  method = maximize_metric, metric = sum_sens_spec)
summary(Beta)
plot(Beta)

#train_data3$binary_score = (train_data3$lasso_scores>=Beta$optimal_cutpoint)
#train_data3$binary_score = (train_data3$lasso_scores>=median(train_data3$lasso_scores))
train_data3$binary_score_tertile = (train_data3$lasso_scores>=quantile(train_data3$lasso_scores,probs=0.66))
#train_data3$binary_score = train_data3$binary_score * 1
train_data3$binary_resp = (train_data3$response == "yes")
write.csv(train_data3, "CD68_survival_data.csv")
train_data_CD68 = read.csv("CD68_survival_1st_last_tertile.csv")

# train_data3$bin_response = (train_data3$response=="yes")*1
# pROC_obj <- roc(train_data3$response,train_data3$lasso_scores,
#                 smoothed = TRUE,
#                 # arguments for ci
#                 ci=TRUE, ci.alpha=0.9, stratified=FALSE,
#                 # arguments for plot
#                 plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
#                 print.auc=TRUE, show.thres=TRUE, direction = '<')
# pROC_obj

#tertile
fit <- survfit(Surv(OS_FROM_START_OF_ITX, VITAL) ~ binary_score_tertile, data = train_data_CD68)
fit$n.risk
# Plot informative survival curves
A = ggsurvplot(fit, data = train_data_CD68,
               title = "CD68 compartment",
               pval = FALSE, pval.method = FALSE,    # Add p-value &  method name
               surv.median.line = "hv",            # Add median survival lines
               legend.title = "Lasso",               # Change legend titles
               legend.labs = c("Low", "High"),  # Change legend labels
               palette = "lancet",                    # Use JCO journal color palette
               risk.table = FALSE,                  # Add No at risk table
               cumevents = FALSE,                  # Add cumulative No of events table
               tables.height = 0.15,               # Specify tables height
               tables.theme = theme_cleantable(),  # Clean theme for tables
               tables.y.text = FALSE               # Hide tables y axis text
)
# find the hazard ratio with this function.
cox <- coxph(Surv(OS_FROM_START_OF_ITX, VITAL) ~ binary_score_tertile, data = train_data_CD68)
cox
ggforest(cox)
A
A$plot <- A$plot+ 
  ggplot2::annotate("text", 
                    x = 1300, y = 0.6, # x and y coordinates of the text
                    label = "HR = 0.5 (0.056-0.75) \n p = 0.008** (Log-Rank)", size = 4)
A

