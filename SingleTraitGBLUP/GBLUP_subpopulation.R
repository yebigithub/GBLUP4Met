# ARGS <- commandArgs(trailingOnly = TRUE)
ia = ARGS[1]
ib = ARGS[2]
treatment0 = 'Control'


#################################################################
path.output = "~/Library/CloudStorage/OneDrive-VirginiaTech/Research/Codes/research/RiceUNLMetabolites/GBLUP4Met/outputs/GBLUP_subpopulation/"
path.geno = "../../Geno" 
subpopulation <- read_csv("../../Met/raw_data/subpopulation.csv")


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
  nCV = 100
  name <- c(paste0("a",1:10),paste0("b",1:10),paste0("c",1:10),paste0("d",1:10),paste0("e",1:10),paste0("f",1:10),paste0("g",1:10),paste0("h",1:3))
  # met2remove=c(13,28,30,38,62,67)
  # name = name[-met2remove]
  # predictive correlation
  corR_GBLUP <- matrix(0, ncol = length(name), nrow = nCV) #GBLUP
  colnames(corR_GBLUP) = name
  subpopulations = unique(subpopulation$Subpopu)
  met0<-list()
  # CV
  for (i in subpopulations) {

    # random-sampling to decide testing & reference accessions
    met0[[i]] = read.csv(paste0("../../Met/CrossValidation/cv_",i,"/met_cv_",i,".csv"))
    met <- met0[[i]] %>% filter(Treatment == treatment)
    index <- which(met$set=="test")# random sampling
    y0 <- dplyr::select(met, -c("NSFTV_ID", "Treatment", "set"))
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
    saveRDS(corR_GBLUP[i, ], file=file.path(path.output, paste0(treatment, "_", method, "_cv_",i,".RDS")))
  }
  # save(corR_GBLUP, file=file.path(path.output, paste0(treatment, "_", method,"_r2_", R20, ".rda")))
}

#################################################################

pred_func(treatment = treatment0, method = "bglr")


