# ARGS <- commandArgs(trailingOnly = TRUE)
# ia = ARGS[1]
# ib = ARGS[2]
# treatment0 = ARGS[3]


#################################################################
path.output = "~/Library/CloudStorage/OneDrive-VirginiaTech/Research/Codes/research/RiceUNLMetabolites/GBLUP4Met/outputs/GBLUP_delta/"
path.geno = "../../Geno" 

library(BGLR)
library(tidyverse)

#################################################################
pred_func <- function(method){

  # G matrix
  load("../../Geno/GL.RData")
  G <- GL[[1]]
  EVD_G <- eigen(G)
  
  # Set up kernel for GBLUP
  ETA1 <- list(
    G = list(V=EVD_G$vectors, d=EVD_G$values, model='RKHS')
  )
  nCV = 100
  
  # predictive correlation
  corR_GBLUP <- matrix(0, ncol = 66, nrow = nCV) #GBLUP
  colnames(corR_GBLUP) = name

  # CV
  for (i in 1:nCV) {

    # random-sampling to decide testing & reference accessions
    met0 = read.csv(paste0("../../Met/CrossValidation/cv_",i,"/met_cv_",i,".csv"))
    met_control0 <- met0 %>% filter(Treatment == 'Control')
    met_control = met_control0[, !(colnames(met_control0) %in% c('NSFTV_ID', 'Treatment', 'set'))]
    met_stress0 <- met0 %>% filter(Treatment == 'Stress')
    met_stress = met_stress0[, !(colnames(met_stress0) %in% c('NSFTV_ID', 'Treatment', 'set'))]
    
    met_delta = data.frame(NSFTV_ID = met_control0$NSFTV_ID, 
                            set = met_control0$set,
                            met_control - met_stress)
    met = met_delta

    name <- colnames(met)[!(colnames(met) %in% c('NSFTV_ID', 'Treatment', 'set', 'Subpopu'))]
    colnames(corR_GBLUP) = name
    
    index <- which(met$set=="test")# random sampling
    y0 <- dplyr::select(met, -c("NSFTV_ID", "set"))
    y <- y0
    y[index, ] <- NA
    
    for (j in 1:ncol(y0)){

      cat("Now running nCV = ", i,"trait = ", colnames(y0)[j], "\n")
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
    saveRDS(corR_GBLUP[i, ], file=file.path(path.output, paste0("delta_cv_",i,".RDS")))
  }
  # save(corR_GBLUP, file=file.path(path.output, paste0(treatment, "_", method,"_r2_", R20, ".rda")))
}

#################################################################
pred_func(method = 'bglr')


