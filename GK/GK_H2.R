############################################################
path.output = "./outputs_H2_GK"
path.geno = "../Geno" 

library(BGLR)
library(tidyverse)

############################################################
H2_func <- function(treatment){
  
  # G matrix
  load("../Geno/GKL.RData")
  G0.25 <- GKL[[treatment]]$GK0.25
  G1.5 <- GKL[[treatment]]$GK1.5
  
  EVD_G0.25 <- eigen(G0.25)
  EVD_G1.5 <- eigen(G1.5)
  
  # Set up kernel for GBLUP
  ETA1 <- list(
    G0.25 = list(V=EVD_G0.25$vectors, d=EVD_G0.25$values, model='RKHS'),
    G1.5 = list(V=EVD_G1.5$vectors, d=EVD_G1.5$values, model='RKHS')
  )
  
  met0<-list()
  H2L<-list()
 # CV
  i = 1
    
    # random-sampling to decide testing & reference accessions
    met0[[i]] = read.csv(paste0("../Met/CrossValidation/cv_",i,"/met_cv_",i,".csv"))
    met <- met0[[i]] %>% filter(Treatment == treatment)
    index <- which(met$set=="test")# random sampling
    y0 <- select(met, -c("NSFTV_ID", "Treatment", "set"))
    y <- y0
    
    for (j in 1:ncol(y0)){
      
      cat("Now running met = ", colnames(y0)[j], "\n")
      ySingle <- y[, j]
      
      #GBLUP
      fit1 <- BGLR(y=ySingle, ETA=ETA1, nIter = 30000,
                   burnIn = 10000, thin = 5, verbose = F)
      
      varE = fit1$varE
      varG0.25 = fit1$ETA$G0.25$varU
      varG1.5 = fit1$ETA$G1.5$varU
      H2 = (varG0.25+varG1.5) / (varG0.25 + varG1.5 + varE)
      H2vec = c(H2, varG0.25, varG1.5, varE); names(H2vec) = c('H2', 'varG0.25', 'varG1.5', 'varE')
      H2L[[j]] = H2vec
      
    }
    save(H2L, file=file.path(path.output, paste0(treatment, "_H2.rda")))
}

############################################################
H2_func(treatment = 'Control')
H2_func(treatment = 'Stress')
