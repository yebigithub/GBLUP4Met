# ARGS <- commandArgs(trailingOnly = TRUE)
# ia = ARGS[1]
# ib = ARGS[2]



#################################################################
path.output = "~/Library/CloudStorage/OneDrive-VirginiaTech/Research/Codes/research/RiceUNLMetabolites/GBLUP4Met/outputs/GBLUP_subpopulation/"
path.geno = "../../Geno" 
library(BGLR)
library(tidyverse)

#################################################################
pred_func <- function(treatment, method){

  # G matrix
  load("../../Geno/GL.RData")
  G <- GL[[treatment]]
  EVD_G <- eigen(G)
  
  # Set up kernel for GBLUP
  ETA1 <- list(
    G = list(V=EVD_G$vectors, d=EVD_G$values, model='RKHS')
  )
  subss = c('aus', 'indica', 'temperate-japonica', 'tropical-japonica')
  
  nCV = as.numeric(length(subss))
  
  # predictive correlation
  corR_GBLUP <- matrix(0, ncol = 66, nrow = nCV) #GBLUP
  
  met0<-list()
  # CV
  for (i in 1:nCV) {

    # random-sampling to decide testing & reference accessions
    met0 <- read.csv(paste0("../../Met/CrossValidation_subpop/cv_",i,"/met_cv_",i,".csv"))
    met <- met0 %>% filter(Treatment == treatment)
    
    name <- colnames(met)[!(colnames(met) %in% c('NSFTV_ID', 'Treatment', 'set', 'Subpopu'))]
    colnames(corR_GBLUP) = name
    
    index <- which(met$set=="test")# random sampling
    y0 <- dplyr::select(met, -c("NSFTV_ID", "Treatment", "set", "Subpopu"))
    y <- y0
    y[index, ] <- NA
    
    for (j in 1:ncol(y0)){

      cat("Now running nCV = ", subss[i],"trait = ", colnames(y0)[j], "\n")
      ySingle <- y[, j]
      
      # #GBLUP
      if(method == 'bglr'){
      fit1 <- BGLR(y=ySingle, ETA=ETA1, nIter = 10000,
                   burnIn = 3000, thin = 5, verbose = F)
      pred1 = fit1$yHat
      corR_GBLUP[i,j] <- cor(pred1[index], y0[index,j])
      }
      
      }
    #save results for every CV iteration.
    saveRDS(corR_GBLUP[i, ], file=file.path(path.output, paste0(treatment, "_", method, "_cv_",i,".RDS")))
  }
  # save(corR_GBLUP, file=file.path(path.output, paste0(treatment, "_", method,"_r2_", R20, ".rda")))
}

#################################################################
for (treatment0 in c('Control', 'Stress')){
  pred_func(treatment = treatment0, method = "bglr")
}

