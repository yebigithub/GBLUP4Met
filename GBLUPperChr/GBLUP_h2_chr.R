path.output = "../temp"
path.geno = "../Geno" 

library(BGLR)
# library(sommer)
library(tidyverse)


h2_chr_func <- function(treatment, method, R2){
  
  # G matrix
  load("../Geno/G_chrLL.RData")
  
  G_chrL = G_chrLL[[treatment]]
  h2LL <- list()
  for (nchr in 1:12){
  
    
    ### For Sommer
    G <- G_chrL[[nchr]]$G_chr
    G_others <- G_chrL[[nchr]]$G_chr_others
    
    ## For BGLR
    EVD_G <- eigen(G_chrL[[nchr]]$G_chr)
    EVD_G_others <- eigen(G_chrL[[nchr]]$G_chr_others)

    ETA2 <- list(
      G_chr = list(V=EVD_G$vectors, d=EVD_G$values, model='RKHS'),
      G_others = list(V=EVD_G_others$vectors, d=EVD_G_others$values, model='RKHS')
    )
    
    
    ## Loading metabolite files.
    met0 <- list()
    i = 1
    # random-sampling to decide testing & reference accessions
    met0[[i]] = read.csv(paste0("../Met/CrossValidation/cv_",i,"/met_cv_",i,".csv"))
    met <- met0[[i]] %>% filter(Treatment == treatment)
    met$NSFTV_IDD <- met$NSFTV_ID
    y0 <- dplyr::select(met, -c("NSFTV_ID", "Treatment", "set", 'NSFTV_IDD'))
    y <- y0
    
    h2L <- list()
    for (j in 1:ncol(y0)){
      
      cat("Now running met =", colnames(y0)[j],"chr =", nchr, "\n")
      
      ySingle <- y[, j]
      
      if (method == 'bglr'){
      #GBLUP per chr
      fit2 <- BGLR(y=ySingle, ETA=ETA2, nIter = 10000,
                   burnIn = 3000, thin = 5, verbose = F, R2=R2)
      
      varE = fit2$varE
      varG = fit2$ETA$G_chr$varU
      varO = fit2$ETA$G_others$varU
      h2 = varG / (varG + varE + varO)
      }
      
      if(method == 'sommer'){
        
        tryCatch({
        ans.ADE <- mmer(as.formula(paste0(colnames(y0)[j],"~1")),
                        random=~vsr(NSFTV_ID,Gu=G) + vsr(NSFTV_IDD,Gu=G_others),
                        rcov=~units, nIters=20,
                        data=met,verbose = FALSE)
        }, error = function(e){
          ans.ADE <<- list()
          cat("ERROR:", conditionMessage(e), '\n')
        })
        
        if(length(ans.ADE) == 0){
          varG = NA
          varE = NA
          varO = NA
          h2 = NA
        }else{
        varG = summary(ans.ADE)$varcomp[1,1]
        varO = summary(ans.ADE)$varcomp[2,1]
        varE = summary(ans.ADE)$varcomp[3,1]
        h2 = varG / (varG + varE + varO)
        
        }
      }
      

      h2vec = c(h2, varG, varO, varE); names(h2vec) = c('h2', 'varG', 'varO', 'varE')
      h2L[[j]] = h2vec
      cat("The h2vec is:", h2vec, "\n")
      
    }
    
    h2LL[[nchr]] = h2L
  }
  save(h2LL, file = file.path(path.output, paste0(treatment, "_h2_chr_R2_", R2, ".rda")))
}


h2_chr_func(treatment = 'Control', method='bglr', R2 = 0.05)
h2_chr_func(treatment = 'Stress', method='bglr', R2 = 0.05)
