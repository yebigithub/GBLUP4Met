boxplot(a1~Batch, data = met_control)



library(tidyverse)

treatment = "Control"
treatment = 'Stress' 

met = met_control
met = met_stress
y = met
pdf(paste0("../temp/", treatment, "_met_boxplot_Batch.pdf"), height = 10, width = 10)

par(mfrow=c(3, 3))
for (j in 5:ncol(y)){
  
  cat("Now running met = ", colnames(y)[j], "\n")
  formu = as.formula(paste0(colnames(y)[j], "~Batch"))
  boxplot(formu, data = y)
  
}

dev.off()