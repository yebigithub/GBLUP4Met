library(tidyverse)

# treatment = "Control"
treatment = 'Stress'

# met = met_rr_Control
met = met_rr_Stress

y0 <- dplyr::select(met, -c("NSFTV_ID", "Treatment"))
y <- y0

# y = y0 = met_rr_control[,-c(1:2)]
# y = y0 = met_rr_stress[,-c(1:2)]

pdf(paste0("../temp/", treatment, "_met_density_plot_lmer.pdf"), height = 10, width = 10)

  par(mfrow=c(3, 3))
  for (j in 1:ncol(y)){
    
    cat("Now running met = ", colnames(y0)[j], "\n")
    ySingle <- y[, j]
    hist(ySingle, breaks = 10, prob = T, main = paste0(colnames(y0)[j], "_", j))
    lines(density(ySingle), lwd = 2, col = 'red')
  }
  
dev.off()

