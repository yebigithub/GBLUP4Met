ARGS <- commandArgs(trailingOnly = TRUE)
ia = ARGS[1]
ib = ARGS[2]
treatment0 = ARGS[3]


############################################################
path.output = "./outputs"
path.geno = "../Geno" 

library(BGLR)
# library(sommer)
library(tidyverse)

############################################################
multitrait_func <- function(treatment, method){
  
  # G matrix
  load("../Geno/GL.RData")
  G <- GL[[treatment]]
  # EVD_G <- eigen(G)
  
  # Set up kernel for GBLUP
  ETA1 <- list(
    G = list(K = G, model='RKHS')
  )
  
  met0<-list()

  # No CV
  i = 1
    # random-sampling to decide testing & reference accessions
    met0[[i]] = read.csv(paste0("../Met/CrossValidation/cv_",i,"/met_cv_",i,".csv"))
    met <- met0[[i]] %>% filter(Treatment == treatment)
    y0 <- select(met, -c("NSFTV_ID", "Treatment", "set"))
    y <- y0
    
    multi_mets = combn(colnames(y0), 2)
    genetic_corrL = list()

    for (ll in ia:ib){
      
      met1 = multi_mets[1, ll]
      met2 = multi_mets[2, ll]
      cat("Now running", met1, met2, "\n")

      ymulti <- cbind(y[, met1], y[, met2])
      #
      # #BGLR-multi-trait
      if(method == "bglr-multi"){

        fit1 <- Multitrait(y = ymulti,
                           ETA = ETA1,
                           verbose = F,
                           nIter = 10000,
                           burnIn = 3000,
                           thin = 5)

        omega = fit1$ETA$G$Cov$Omega
        corr = omega[1,2]/(sqrt(omega[1,1]*omega[2,2]))
      }
      
      # #sommer-multi-trait
      # if(method == "sommer-multi"){
      #   
      #   tryCatch({
      #   ans.m <- mmer(as.formula(paste0('cbind(',met1,',',met2,')~1')),
      #                 random=~ vsr(NSFTV_ID, Gu=G, Gtc=unsm(2)),
      #                 rcov=~ vsr(units, Gtc=unsm(2)), nIters=20,
      #                 data=ymet, verbose = FALSE, tolParInv = 0)
      #   }, error=function(e){
      #     ans.m <<- list()
      #     cat("ERROR :",conditionMessage(e), "\n")
      #   })
      #   
      #   if(length(ans.m) == 0){
      #     corr = NA
      #   }else{corr = cov2cor(ans.m$sigma$`u:NSFTV_ID`)[1,2]}
      # 
      #   }
      
        corr_vec = c(met1, met2, corr); names(corr_vec) = c('met1', 'met2', 'corr')
        genetic_corrL[[ll]] = corr_vec
        save(corr_vec, file=file.path(path.output, paste0(treatment,"_", method, "_", met1, "_", met2, "_genetic_corr.rda")))
        cat(corr_vec)
        cat("\n ########################################################## \n")
    }
    # save(genetic_corrL, file=file.path(path.output, paste0(treatment,"_", method, "_genetic_corr.rda")))
}

############################################################

multitrait_func(treatment = treatment0, method = 'bglr-multi')


