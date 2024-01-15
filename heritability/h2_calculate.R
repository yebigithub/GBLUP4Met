############################################################
## Sommer to calculate heritability
############################################################

path.output = "../temp"
path.geno = "../Geno" 

library(tidyverse)
library(sommer)


############################################################
h2_func <- function(treatment, method){
  
  # G matrix
  load("../Geno/GL.RData")
  G <- GL[[treatment]]
  EVD_G <- eigen(G)
  
  met0<-list()
  h2L<-list()
  # # CV
  i = 1
  
  # random-sampling to decide testing & reference accessions
  met0[[i]] = read.csv(paste0("../Met/CrossValidation/cv_",i,"/met_cv_",i,".csv"))
  met <- met0[[i]] %>% filter(Treatment == treatment)
  # index <- which(met$set=="test")# random sampling
  input = met
  y0 <- dplyr::select(met, -c("NSFTV_ID", "Treatment", "set"))
  y <- y0
  
  for (j in 1:ncol(y)){
    
    cat("Now running met = ", colnames(y0)[j], "\n")
    ySingle <- y[, j]
    
    if(method == "sommer"){
      
      fit_ans <- mmer(as.formula(paste0(colnames(y0)[j], "~1")),
                      random=~vsr(NSFTV_ID, Gu=G),
                      rcov=~units,nIters=20,
                      data=input, verbose = FALSE) 
      
      varE = fit_ans$sigma$units
      varG = fit_ans$sigma$`u:NSFTV_ID`
      h2 = varG / (varG + varE)
      h2vec = c(h2, varG, varE); names(h2vec) = c('h2', 'varG', 'varE')
      h2L[[j]] = h2vec       
    }
    
    }
  save(h2L, file=file.path(path.output, paste0(treatment,"_", method, "_h2.rda")))
}

############################################################
h2_func(treatment = 'Control', method = 'sommer')
h2_func(treatment = 'Stress', method = 'sommer')
############################################################

