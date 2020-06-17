rm(list=ls())
setwd("../7_EvolRate_Models/")

outputMatrix <- matrix(0, nrow = 12, ncol = 9)
colnames(outputMatrix) <- c('Variable', 'Pearson r', 'Pval_Pearson','Spearman rho', 'Pval_spearman', "Student's t", 'Pval student', "Wilcoxon's Z", 'Pval Wilcox')

# First run script "1_Bayesian_analysis.R" to obtain 'rates_VARIABLE.csv' files (unless existing files in root directory)

rates_PC1 = log10(read.csv("rates_PC1.csv", header = TRUE, row.names = 1, col.names = c('','rates_PC1')))
rates_PC2 = log10(read.csv("rates_PC2.csv", header = TRUE, row.names = 1, col.names = c('','rates_PC2')))
rates_PC3 = log10(read.csv("rates_PC3.csv", header = TRUE, row.names = 1, col.names = c('','rates_PC3')))
rates_PC4 = log10(read.csv("rates_PC4.csv", header = TRUE, row.names = 1, col.names = c('','rates_PC4')))
rates_PC5 = log10(read.csv("rates_PC5.csv", header = TRUE, row.names = 1, col.names = c('','rates_PC5')))
rates_PC6 = log10(read.csv("rates_PC6.csv", header = TRUE, row.names = 1, col.names = c('','rates_PC6')))
rates_LD1 = log10(read.csv("rates_LD1.csv", header = TRUE, row.names = 1, col.names = c('','rates_LD1')))
rates_LD2 = log10(read.csv("rates_LD2.csv", header = TRUE, row.names = 1, col.names = c('','rates_LD2')))
rates_LD3 = log10(read.csv("rates_LD3.csv", header = TRUE, row.names = 1, col.names = c('','rates_LD3')))
rates_LD4 = log10(read.csv("rates_LD4.csv", header = TRUE, row.names = 1, col.names = c('','rates_LD4')))
rates_LD5 = log10(read.csv("rates_LD5.csv", header = TRUE, row.names = 1, col.names = c('','rates_LD5')))
rates_LD6 = log10(read.csv("rates_LD6.csv", header = TRUE, row.names = 1, col.names = c('','rates_LD6')))
rates_LI = log10(read.csv("rates_localInfo.csv", header = TRUE, row.names = 1, col.names = c('','rates_LI')))

all_rates_PC <- as.data.frame(cbind(rates_PC1,rates_PC2,rates_PC3,rates_PC4,rates_PC5,rates_PC6, rates_LI))
all_rates_LD <- as.data.frame(cbind(rates_LD1,rates_LD2,rates_LD3,rates_LD4,rates_LD5,rates_LD6, rates_LI))



outputMatrix[1:6,1] <- colnames(all_rates_PC)[1:6]
outputMatrix[7:12,1] <- colnames(all_rates_LD)[1:6]



Shap_rates_PC <- lapply(all_rates_PC, shapiro.test) 
Shap_rates_LD <- lapply(all_rates_LD, shapiro.test) 


correlations_PCrates <- cor(all_rates_PC, use="complete.obs", method="pearson")
correlations_LDrates <- cor(all_rates_LD, use="complete.obs", method="pearson")

pairs(all_rates_PC)
pairs(all_rates_LD)


#pdf('rates_PC_LI.pdf')
par(mar=c(5,5,5,5))
counter_rowMatrix <- 0
counter_colName <- 0

for (variable in all_rates_PC[1:6]){
  counter_colName <- counter_colName +1
  plot(variable, all_rates_PC$rates_LI, xlab=colnames(all_rates_PC)[counter_colName], ylab= 'LI rates', cex.axis= 1.4, cex.lab=2, font.lab=2, main = NULL)
  abline(lm(rates_LI~variable, data <- all_rates_PC), col = 'red', lw=2)
  cor_pearson <- cor.test(variable, all_rates_PC$rates_LI, method = 'pearson')
  cor_spearman <- cor.test(variable, all_rates_PC$rates_LI, method = 'spearman')
  mtext(paste0("Spearman's Rho = ",round(cor_spearman$estimate, digits = 2),'\n'), 3, adj=1, font = 4, cex = 1.5)
  counter_rowMatrix <- counter_rowMatrix+1
  outputMatrix[counter_rowMatrix, 2:5] <- c(cor_pearson$estimate,
                            cor_pearson$p.value,
                            cor_spearman$estimate,
                            cor_spearman$p.value)
  studentTest <- t.test(variable, all_rates_PC$rates_LI)
  wilcoxTest <- wilcox.test(variable, all_rates_PC$rates_LI,paired = T)
  outputMatrix[counter_rowMatrix,6:9] <- c(studentTest$statistic, studentTest$p.value, qnorm(wilcoxTest$p.value/2), wilcoxTest$p.value)  
  
  }
#dev.off()

#pdf('rates_LD_LI.pdf')
par(mar=c(5,5,5,5))
#counter_rowMatrix <- 0 #this counter needs not be reset, as we keep feeding data into main matrix across loops
counter_colName <- 0

for (variable in all_rates_LD[1:6]){
  counter_colName <- counter_colName +1
  plot(variable, all_rates_LD$rates_LI, xlab=colnames(all_rates_LD)[counter_colName], ylab= 'LI rates', cex.axis= 1.4, cex.lab=2, font.lab=2, main = NULL)
  abline(lm(rates_LI~variable, data <- all_rates_LD), col = 'red', lw=2)
  cor_pearson <- cor.test(variable, all_rates_LD$rates_LI, method = 'pearson')
  cor_spearman <- cor.test(variable, all_rates_LD$rates_LI, method = 'spearman')
  mtext(paste0("Spearman's Rho = ",round(cor_spearman$estimate, digits = 2),'\n'), 3, adj=1, font = 4, cex = 1.5)
  counter_rowMatrix <- counter_rowMatrix+1
  outputMatrix[counter_rowMatrix, 2:5] <- c(cor_pearson$estimate,
                                            cor_pearson$p.value,
                                            cor_spearman$estimate,
                                            cor_spearman$p.value)
  studentTest <- t.test(variable, all_rates_LD$rates_LI)
  wilcoxTest <- wilcox.test(variable, all_rates_LD$rates_LI,paired = T)
  outputMatrix[counter_rowMatrix,6:9] <- c(studentTest$statistic, studentTest$p.value, qnorm(wilcoxTest$p.value/2), wilcoxTest$p.value)  
  
}
#dev.off()


write.csv(outputMatrix, 'OutputMatrix.csv')
p.adjust(as.vector(outputMatrix[,7]),method = 'bonferroni', n=length(outputMatrix[,7])) 

p.adjust(as.vector(outputMatrix[,9]),method = 'bonferroni', n=length(outputMatrix[,9]))
