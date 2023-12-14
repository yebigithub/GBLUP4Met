############################################################
path.output = "../temp"
path.geno = "../Geno" 

library(BGLR)
library(rrBLUP)
library(tidyverse)
library(sommer)
# library(gaston)
# library(EMMREML)

############################################################
h2_func <- function(treatment, method){
  
  # G matrix
  load("../Geno/GL.RData")
  G <- GL[[treatment]]
  EVD_G <- eigen(G)
  
  # Set up kernel for GBLUP
  # ETA1 <- list(
  #   # G = list(V=EVD_G$vectors, d=EVD_G$values, model='RKHS')
  #   G = list(K = G, model = 'RKHS')
  # )
  
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
      
      #GBLUP
      if(method == "bglr"){
        
        fit1 <- BGLR(y=ySingle, ETA=ETA1, nIter = 10000,
                     burnIn = 3000, thin = 1, verbose = F, R2=0.05)
        
        varE = fit1$varE
        varG = fit1$ETA$G$varU
        h2 = varG / (varG + varE)
        h2vec = c(h2, varG, varE); names(h2vec) = c('h2', 'varG', 'varE')
        h2L[[j]] = h2vec
      }
      
      #rrBLUP
      if(method == "rrblup"){
          
        fit2 <- mixed.solve(y=ySingle, K=G, method = 'REML')
        
        varE = fit2$Ve
        varG = fit2$Vu
        h2 = varG / (varG + varE)
        h2vec = c(h2, varG, varE); names(h2vec) = c('h2', 'varG', 'varE')
        h2L[[j]] = h2vec
      }
      
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
      
      # if(method == 'gaston'){
      #   
      #   fit3 <- lmm.aireml(Y=ySingle, K=G, EMsteps=5)
      #   
      #   varE = fit3$sigma2
      #   varG = fit3$tau
      #   h2 = varG / (varG + varE)
      #   h2vec = c(h2, varG, varE); names(h2vec) = c('h2', 'varG', 'varE')
      #   h2L[[j]] = h2vec
      # }
      #   
      # if(method == 'emmreml'){
      #   
      #   fit4 = emmreml(y=ySingle, X=matrix(rep(1, length(ySingle)), ncol=1), Z=diag(length(ySingle)), K=G)
      #   
      #   varE = fit4$Ve
      #   varG = fit4$Vu
      #   h2 = varG / (varG + varE)
      #   h2vec = c(h2, varG, varE); names(h2vec) = c('h2', 'varG', 'varE')
      #   h2L[[j]] = h2vec
      #   
      # }
        
    }
    save(h2L, file=file.path(path.output, paste0(treatment,"_", method, "_h2.rda")))
}

############################################################
h2_func(treatment = 'Control', method = 'rrblup')
h2_func(reatment = 'Stress', method = 'rrblup')

h2_func(treatment = 'Control', method = 'bglr')
h2_func(treatment = 'Stress', method = 'bglr')

h2_func(treatment = 'Control', method = 'sommer')
h2_func(treatment = 'Stress', method = 'sommer')



#this loading is for testing preprocessing method
# input_load = function(method, treatment){
#   
#   met_rr = read.csv(paste0("../temp/", method, "/met_rr_", treatment, ".csv"))
#   return(met_rr)
#   
# }
# met_rr_Control = input_load("log_lmer", "Control")
# met_rr_Stress = input_load("log_lmer", "Stress")
