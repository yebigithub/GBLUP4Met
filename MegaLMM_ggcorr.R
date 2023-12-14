# ARGS <- commandArgs(trailingOnly = TRUE)
# ia = ARGS[1]
# ib = ARGS[2]
# treatment0 = ARGS[3]
#Select specific R version
# .libPaths()
# .libPaths(.libPaths()[3:1])
# 
# #install packages from github
# install.packages("devtools", repos = "http://cran.us.r-project.org")
# devtools::install_github('deruncie/MegaLMM')

#################################################################
path.output = "./outputs"
path.geno = "../Geno" 

library(MegaLMM)
library(tidyverse)

#################################################################
MegaLMM_pred_func <- function(treatment = "Control"){
  
  # G matrix
  load("../Geno/GL.RData")
  G <- GL[[treatment]]
  
  # load("../Geno/GKL.RData")
  # G0.25 <- GKL[[treatment]]$GK0.25
  # G1.5 <- GKL[[treatment]]$GK1.5
  
  # Set up kernel for GBLUP
  # ETA1 <- list(
  #   G = list(V=EVD_G$vectors, d=EVD_G$values, model='RKHS')
  # )
  name <- c(paste0("a",1:10),paste0("b",1:10),paste0("c",1:10),paste0("d",1:10),paste0("e",1:10),paste0("f",1:10),paste0("g",1:10),paste0("h",1:3))
  met2remove=c(13,28,30,33,38,62,67)
  name = name[-met2remove]
  # predictive correlation
  # corRL <- list()
  met0<-list()
  # CV
  for (i in 1:1) {
    # i=1
    # random-sampling to decide testing & reference accessions
    met0[[i]] = read.csv(paste0("../Met/CrossValidation/cv_",i,"/met_cv_",i,".csv"))
    met.wide0 <- met0[[i]] %>% filter(Treatment == treatment)
    met.wide <- dplyr::select(met.wide0, all_of(c(name, "NSFTV_ID", "set")))
    met.long <- tidyr::gather(met.wide, 
                              key="met",
                              value="value",
                              -c(NSFTV_ID,set),
                              factor_key = F)
    # met.long$NSFTV_IDD = met.long$NSFTV_ID
    # 
    # data_matrices = create_data_matrices(
    #   tall_data = met.long, 
    #   id_cols = c('NSFTV_ID','NSFTV_IDD', "set"), 
    #   names_from = 'met', 
    #   values_from = 'value' 
    # )
    
    data_matrices = create_data_matrices(
      tall_data = met.long, 
      id_cols = c('NSFTV_ID', "set"), 
      names_from = 'met', 
      values_from = 'value' 
    )
    
    
    Y = data_matrices$Y
    sample_data = data_matrices$data
    id_train = sample_data$NSFTV_ID[sample_data$set=="train"]
    id_test = sample_data$NSFTV_ID[sample_data$set=="test"]
    
    fold_ID = i
    Y_train = Y_testing = Y
    # Y_train[id_test, ] = NA
    # Y_testing[id_train, ] = NA

    run_parameters = MegaLMM_control(
      h2_divisions = 20, 
      burn = 0,  
      thin = 1,
      K = nnFA
      )
    
    # MegaLMM_state = setup_model_MegaLMM(
    #   Y = Y_train,
    #   formula = ~ (1|NSFTV_ID) + (1|NSFTV_IDD),
    #   data = sample_data,
    #   relmat = list(NSFTV_ID = G0.25,
    #                 NSFTV_ID = G1.5),
    #   run_parameters=run_parameters,
    #   run_ID = sprintf('MegaLMM_fold_%02d',fold_ID)
    # )
    
    MegaLMM_state = setup_model_MegaLMM(
      Y = Y_train,
      formula = ~ (1|NSFTV_ID),
      data = sample_data,
      relmat = list(NSFTV_ID = G),
      run_parameters=run_parameters,
      run_ID = sprintf('MegaLMM_fold_%02d',fold_ID)
    )
    
    Lambda_prior = list(
      sampler = sample_Lambda_prec_horseshoe, 
      prop_0 = 0.1,    
      delta = list(shape = 3, scale = 1),    
      delta_iterations_factor = 100   
    )
    
    priors = MegaLMM_priors(
      tot_Y_var = list(V = 0.5,   nu = 5),      
      tot_F_var = list(V = 18/20, nu = 20),     
      h2_priors_resids_fun = function(h2s,n)  1,  
      h2_priors_factors_fun = function(h2s,n) 1, 
      Lambda_prior = Lambda_prior
    )
    
    MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
    MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
    MegaLMM_state = initialize_MegaLMM(MegaLMM_state,verbose = T)
    MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2','tot_Eta_prec')
    MegaLMM_state$Posterior$posteriorMean_params = 'Eta_mean'
    MegaLMM_state$Posterior$posteriorFunctions = list(
      U = 'U_F %*% Lambda + U_R',
      G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
      R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])',
      h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])'
    )
    MegaLMM_state = clear_Posterior(MegaLMM_state) 
    
    n_iter = 250
    for(jj in 1:40) {
      print(sprintf('Sampling run %d',jj))
      MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter) 
      MegaLMM_state = save_posterior_chunk(MegaLMM_state)
      print(MegaLMM_state)
    }
    

    G_samples = load_posterior_param(MegaLMM_state,'G')
    
    # saveRDS(corR_GBLUP[i, ], file=file.path(path.output, paste0(treatment, "_", method, "_cv_",i,".RDS")))
  }

  return(G_samples)
}

#################################################################
nnFA = 4
varcov_control <- MegaLMM_pred_func(treatment = "Control")
save(varcov_control, file="../temp/varcov_control.RDS")

nnFA = 5
varcov_stress <- MegaLMM_pred_func(treatment = "Stress")
save(varcov_stress, file = "../temp/varcov_stress.RDS")

