# ARGS <- commandArgs(trailingOnly = TRUE)
# ia = ARGS[1]
# ib = ARGS[2]
# treatment0 = ARGS[3]


############################################################
path.output = "../temp"
path.geno = "../Geno" 

library(BGLR)
# library(sommer)
library(tidyverse)

############################################################
multitrait_func <- function(method, R2){
  
  # G matrix
  load("../Geno/Gop.RData")
  G <- G_op
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
  met = met0[[i]]
  met_con = met %>% filter(Treatment == "Control") %>% droplevels()
  met_trt = met %>% filter(Treatment == "Stress") %>% droplevels()
  met_con_op <- met_con[met_con$NSFTV_ID %in% rownames(G), ]
  met_trt_op <- met_trt[met_trt$NSFTV_ID %in% rownames(G), ]
  

  y1 <- select(met_con_op, -c("NSFTV_ID", "Treatment", "set"))
  y2 <- select(met_trt_op, -c("NSFTV_ID", "Treatment", "set"))


  genetic_corrL = list()
  
  for (ll in 1:ncol(y1)){
    
    met1 = colnames(y1)[ll]
    cat("Now running", met1, "\n")
    
    ymulti <- cbind(y1[, met1], y2[, met1])
    #
    # #BGLR-multi-trait
    if(method == "bglrmulti"){
      
      fit1 <- Multitrait(y = ymulti,
                         ETA = ETA1,
                         verbose = F,
                         nIter = 50000,
                         burnIn = 10000,
                         thin = 5,
                         R2 = R2)
      
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
    
    corr_vec = c(met1, met1, corr); names(corr_vec) = c('met1', 'met2', 'corr')
    genetic_corrL[[ll]] = corr_vec
    # save(corr_vec, file=file.path(path.output, paste0(treatment,"_", method, "_", met1, "_", met2, "_genetic_corr.rda")))
    cat(corr_vec)
    cat("\n ########################################################## \n")
  }
  save(genetic_corrL, file=file.path(path.output, paste0(method, "_genetic_corr_perTreatment", R2, ".rda")))
}

############################################################

multitrait_func(method="bglrmulti", R2 = 0.05)


