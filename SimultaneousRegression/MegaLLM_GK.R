
#################################################################
path.geno = "../Geno" 

# devtools::install_github('deruncie/MegaLMM')

library(MegaLMM)
library(tidyverse)

#################################################################
pred_GK_func <- function(treatment, nFA, bw_mean){
  
  
  path.output = paste0("../temp/mega/GK/", treatment, "_", bw_mean)
  
  if(!dir.exists(path.output)){
    dir.create(path.output, recursive = TRUE)
    cat("Directory created successfully:", path.output, "\n")
    }
  # G matrix
  
  load("../Geno/GKL.RData")
  GKL <- GKL[[treatment]]
  G = GKL[[bw_mean]]
  
   nCV = 100
  name <- c(paste0("a",1:10),paste0("b",1:10),paste0("c",1:10),paste0("d",1:10),paste0("e",1:10),paste0("f",1:10),paste0("g",1:10),paste0("h",1:3))
  met2remove=c(13,28,30,33,38,62,67)
  name = name[-met2remove]
  # predictive correlation
  corRL <- list()
  met0<-list()
  # CV
  for (i in 18:18) {
    # i=1
    
    RDS_path = file.path(path.output, paste0(treatment, "_mega_cv_",i,".RDS"))
    
    if(file.exists(RDS_path)){
      cat(RDS_path," already exists \n")
      next
    }
    cat("Now is running cv", i, "******************************************************************************************************************************************** \n")
    # random-sampling to decide testing & reference accessions
    met0[[i]] = read.csv(paste0("../Met/CrossValidation/cv_",i,"/met_cv_",i,".csv"))
    met.wide0 <- met0[[i]] %>% filter(Treatment == treatment)
    met.wide <- dplyr::select(met.wide0, all_of(c(name, "NSFTV_ID", "set")))
    met.long <- tidyr::gather(met.wide, 
                              key="met",
                              value="value",
                              -c(NSFTV_ID,set),
                              factor_key = F)
     
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
    Y_train[id_test, ] = NA
    Y_testing[id_train, ] = NA
    
    run_parameters = MegaLMM_control(
      h2_divisions = 20, 
      burn = 0,  
      thin = 1,
      K = nFA 
    )

    MegaLMM_state = setup_model_MegaLMM(
      Y = Y_train,
      formula = ~ (1|NSFTV_ID),
      data = sample_data,
      relmat = list(NSFTV_ID = G),
      run_parameters=run_parameters,
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
    for(jj in 1:28) {
      print(sprintf('Sampling run %d',jj))
      MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter) 
      MegaLMM_state = save_posterior_chunk(MegaLMM_state)
      print(MegaLMM_state)
    }
    
    # Lambda_samples = load_posterior_param(MegaLMM_state,'Lambda')
    U_samples = load_posterior_param(MegaLMM_state,'U')
    
    U_hat = apply(U_samples[4000:7000,,], c(2,3), mean)
    # Eta_mean = load_posterior_param(MegaLMM_state,'Eta_mean')
    MegaLMM_Uhat_accuracy = diag(cor(Y_testing,U_hat,use='p'))
    # MegaLMM_Eta_mean_accuracy = diag(cor(Y_testing,Eta_mean,use='p'))
     
    corRL[[i]] = MegaLMM_Uhat_accuracy
    
    saveRDS(MegaLMM_Uhat_accuracy, file=file.path(path.output, paste0(treatment, "_mega_cv_",i,".RDS")))
  }
  save(corRL, file=file.path(path.output, paste0(treatment, "_MegaLMM.rda")))
}

#################################################################
nFA=5
pred_GK_func(treatment = 'Control', nFA = 5, bw_mean = "GK_0.8")
pred_GK_func(treatment = 'Control', nFA = 5, bw_mean = "GK_0.6")
pred_GK_func(treatment = 'Control', nFA = 5, bw_mean = "GK_0.4")
pred_GK_func(treatment = 'Control', nFA = 5, bw_mean = "GK_0.2")



nFA=5
pred_GK_func(treatment = 'Stress', nFA = 5, bw_mean = "GK_0.8")
pred_GK_func(treatment = 'Stress', nFA = 5, bw_mean = "GK_0.6")
pred_GK_func(treatment = 'Stress', nFA = 5, bw_mean = "GK_0.4")
pred_GK_func(treatment = 'Stress', nFA = 5, bw_mean = "GK_0.2")



