ARGS <- commandArgs(trailingOnly = TRUE)
ia = ARGS[1]
ib = ARGS[2]
treatment0 = ARGS[3]

#################################################################
path.output = "./outputs_GK"
path.geno = "../Geno" 

library(BGLR)
library(tidyverse)

#################################################################
pred_func <- function(treatment){

  # G matrix
  load("../Geno/GKL.RData")
  G0.25 <- GKL[[treatment]]$GK0.25
  G1.5 <- GKL[[treatment]]$GK1.5
  
  EVD_G0.25 <- eigen(G0.25)
  EVD_G1.5 <- eigen(G1.5)
  
  # Set up kernel for GK
  ETA1 <- list(
    G0.25 = list(V=EVD_G0.25$vectors, d=EVD_G0.25$values, model='RKHS'),
    G1.5 = list(V=EVD_G1.5$vectors, d=EVD_G1.5$values, model='RKHS')
    
  )
  
  nCV = 100
  name <- c(paste0("a",1:10),paste0("b",1:10),paste0("c",1:10),paste0("d",1:10),paste0("e",1:10),paste0("f",1:10),paste0("g",1:10),paste0("h",1:3))
  met2remove=c(13,28,30,38,62,67)
  name = name[-met2remove]
  
  # predictive correlation
  corR_GK <- matrix(0, ncol = length(name), nrow = nCV) #GK
  colnames(corR_GK) = name
  
  met0<-list()
  # CV
  for (i in ia:ib) {

    # random-sampling to decide testing & reference accessions
    met0[[i]] = read.csv(paste0("../Met/CrossValidation/cv_",i,"/met_cv_",i,".csv"))
    met <- met0[[i]] %>% filter(Treatment == treatment)
    index <- which(met$set=="test")# random sampling
    y0 <- select(met, -c("NSFTV_ID", "Treatment", "set"))
    y <- y0
    y[index, ] <- NA
    for (j in 1:ncol(y0)){

      cat("Now running nCV = ", i,"trait = ", colnames(y0)[j], "\n")
      ySingle <- y[, j]
      
      #GK
      fit1 <- BGLR(y=ySingle, ETA=ETA1, nIter = 10000,
                   burnIn = 3000, thin = 5, verbose = F)
      pred1 = fit1$yHat
      corR_GK[i,j] <- cor(pred1[index], y0[index,j])
      
    }
    #save results for every CV iteration.
    saveRDS(corR_GK[i, ], file=file.path(path.output, paste0(treatment, "_cv_",i,".RDS")))
  }
  
  # save(corR_GBLUP, file=file.path(path.output, paste0(treatment, "_GBLUP.rda")))
}

#################################################################

pred_func(treatment = treatment0)
