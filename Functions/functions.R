
############################################################
met_name_func <- function(mode){
met_names0 <- read.delim("/Users/yebi/Library/CloudStorage/OneDrive-VirginiaTech/Research/Codes/research/RiceUNLMetabolites/GBLUP4Met/Met/raw_data/met_names.txt", header=FALSE)[,1]
if(mode == "alphabet"){
name = c(paste0("a",1:10),paste0("b",1:10),paste0("c",1:10),paste0("d",1:10),
         paste0("e",1:10),paste0("f",1:10),paste0("g",1:10),paste0("h",1:3))
met_names0 = name
}
met2remove = c(13,28,30,33,38,62,67)
met_names = met_names0[-met2remove]
return(met_names)
}

############################################################
load_pred_corr_func <- function(path, namm, kernel){
  PredCorr_conL <- list()
  PredCorr_trtL <- list()
  for (i in 1:100){
    PredCorr_conL[[i]] <- readRDS(file=file.path(path,paste0("Control_",namm,"cv_",i,".RDS")))
    PredCorr_trtL[[i]] <- readRDS(file=file.path(path,paste0("Stress_",namm,"cv_",i,".RDS")))
  }
  
  corr.c.df <- do.call(rbind, PredCorr_conL)
  corr.c.df = corr.c.df[,1:66]
  colnames(corr.c.df)= met_names
  
  corr.t.df <- do.call(rbind, PredCorr_trtL)
  corr.t.df = corr.t.df[,1:66]
  colnames(corr.t.df) = met_names
  
  temp1 <- reshape2::melt(corr.c.df)
  temp2 <- reshape2::melt(corr.t.df)
  temp3 <- rbind.data.frame(temp1,temp2)
  temp3$group <- rep(c("Control", "Stress"), each = 6600)
  colnames(temp3) <- c("CV", "Met", "Corr", "Treatment")
  temp3$Kernel = kernel
  temp3$Met = factor(temp3$Met, levels = met_names)
  temp3 = temp3 %>% arrange(CV, Treatment)
  
  return(temp3)
}

############################################################
density_plot_generator = function(df=df, 
                                 my_colors=RColorBrewer::brewer.pal(11, 'Spectral')[c(10,3)],
                                 limits_range = limits_range,
                                 breaks_range = breaks_range, 
                                 x_name = "Heritability",
                                 lines = T,
                                 ylim){

   
  mu <- ddply(df, "Treatment", summarise, grp.mean=mean(value), grp.me = median(value))
  # head(mu)
  # my_colors <- RColorBrewer::brewer.pal(11, 'Spectral')[c(10,3)]
  
  ppp = ggplot(df, aes(x=value, color=Treatment, fill = Treatment)) + 
        geom_density(alpha=0.6)+
        scale_x_continuous(limits=limits_range, breaks=breaks_range) +
        labs(x=x_name, y = "Density") +
        scale_fill_manual(values = my_colors)+
        scale_color_manual(values = my_colors)+
        theme_bw()

  if(lines){
  ppp = ppp+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Treatment), linetype='solid')+
    geom_vline(data=mu, aes(xintercept=grp.me, color=Treatment), linetype='dashed')
  }
  if(length(ylim) != 0){
    ppp = ppp+ 
      coord_cartesian(ylim=ylim)
  }
  return(ppp)
}



############################################################
FAdf_prep_func = function(df, treatment, FmetL){
  dfcon = df %>% filter(Treatment == treatment) %>% droplevels()
  dfcon$Factor = NA
  for(i in 1:ncol(FmetL)){
    Factor = colnames(FmetL)[i]
    ix = FmetL[,i]
    dfcon$Factor[ix] = Factor 
  }
  dfcon = na.omit(dfcon)
  dfcon$Factor = factor(dfcon$Factor, levels=paste0("F",1:10))
  return(dfcon)
}

############################################################

FA_density_plot_generator = function(df,
                                     my_colors = RColorBrewer::brewer.pal(11, 'Spectral')[c(11,10,9,8,7,5,4,3,2,1)],
                                     limits_range = limits_range,
                                     breaks_range = breaks_range,
                                     x_name = "Heritability",
                                     ylim){
  df = df
  FmetL_con = readRDS("/Users/yebi/Library/CloudStorage/OneDrive-VirginiaTech/Research/Codes/research/RiceUNLMetabolites/GBLUP4Met/Met/FmetL_con.RDS")
  FmetL_trt = readRDS("/Users/yebi/Library/CloudStorage/OneDrive-VirginiaTech/Research/Codes/research/RiceUNLMetabolites/GBLUP4Met/Met/FmetL_trt.RDS")
  dfcon = FAdf_prep_func(df, 'Control', FmetL_con) ####functions.R
  dftrt = FAdf_prep_func(df, "Stress", FmetL_trt) ####functions.R
  
  dff = rbind.data.frame(dfcon, dftrt)
  # summary(dff$value)
  
  fa_ppp = ggplot(dff, aes(x=value, color=Factor, fill = Factor)) + 
    geom_density(alpha=0.6)+
    scale_x_continuous(limits=limits_range, breaks=breaks_range) +
    labs(x=x_name, y = "Density") +
    facet_grid(rows = vars(Treatment))+
    scale_fill_manual(values = my_colors)+
    scale_color_manual(values = my_colors)+
    theme_bw()
  
  if(length(ylim) != 0){
    fa_ppp = fa_ppp +
             coord_cartesian(ylim=ylim)
  }
  return(fa_ppp)
}

############################################################
pred_mean_generator = function(PreCorr){
  mmm_df = plyr::ddply(PreCorr, c("Treatment", "Met"), summarise, corr.mean=mean(Corr, na.rm=T))
  return(mmm_df)
}

############################################################
pred_dff_generator = function(df1, df2, Method){
          dff_df = data.frame(Met=df1$Met,
                                value = (df2$corr.mean-df1$corr.mean)/df1$corr.mean*100,
                              # value = df2$corr.mean,
                                Treatment = df1$Treatment,
                                Method = Method)
          
          return(dff_df)
}



####################################
FA_summary_func = function(df){
  df = df
  FmetL_con = readRDS("/Users/yebi/Library/CloudStorage/OneDrive-VirginiaTech/Research/Codes/research/RiceUNLMetabolites/GBLUP4Met/Met/FmetL_con.RDS")
  FmetL_trt = readRDS("/Users/yebi/Library/CloudStorage/OneDrive-VirginiaTech/Research/Codes/research/RiceUNLMetabolites/GBLUP4Met/Met/FmetL_trt.RDS")
  dfcon = FAdf_prep_func(df, 'Control', FmetL_con) ####functions.R
  dftrt = FAdf_prep_func(df, "Stress", FmetL_trt)
  dff = rbind.data.frame(dfcon, dftrt)
  
  dff_FA = ddply(dff, c("Factor"), summarise, 
                 Mean = round(mean(value),2),
                 Median = round(median(value),2),
                 Min = round(min(value),2),
                 Max = round(max(value),2))
  return(dff_FA)
}



##############################################
compare_plot_generator = function(df, low_thr, up_thr, limits_range, breaks_range){
  df_wide = pivot_wider(df, names_from = Treatment, values_from = value)
  df_wide$color = "grey40"
  aaa = abs(df_wide$Control-df_wide$Stress)
  df_wide$color[aaa>=low_thr] = "brown"
  df_wide$color[aaa<=up_thr] = "darkgreen"
  
  ggplot(df_wide, aes(x=Control, y=Stress)) + 
    geom_abline(slope=1, intercept=0, colour = "grey", linewidth=0.7)+
    # geom_point(colour = "midnightblue", alpha = 0.8) +
    geom_point(color = df_wide$color, alpha = 0.7, size = 1)+
    geom_point(color = df_wide$color, alpha = 0.7, size = 1)+
    theme_bw()+
    scale_y_continuous(limits=limits_range, breaks=breaks_range) +
    scale_x_continuous(limits=limits_range, breaks=breaks_range) +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8))+
    geom_text_repel(data=subset(df_wide, abs(Control-Stress)>=low_thr),
                    aes(Control,Stress,label=Met), col="brown", size = 2,
                    box.padding = 0.3,
                    min.segment.length = 0.1,
                    arrow = arrow(length = unit(0.008, "npc"), 
                                  type = "closed",
                                  ends = "first")) +
    geom_text_repel(data=subset(df_wide, abs(Control-Stress)<=up_thr),
                    aes(Control,Stress,label=Met), col="darkgreen", size = 2,
                    box.padding = 0.3,
                    min.segment.length = 0.1,
                    # point.padding = 2,
                    # segment.linetype = 2,
                    # segment.curvature = -1e-1,
                    arrow = arrow(length = unit(0.008, "npc"), 
                                  type = "open",
                                  ends = "first"))
}

#########################################
DL_extractor_func <- function(df, met_names, kernel){
  
  df_control = df[,1:66];colnames(df_control)=met_names;df_control$CV=1:100
  df_stress = df[,67:132];colnames(df_stress)=met_names;df_stress$CV=1:100
  df_control = pivot_longer(df_control, names_to = "Met", values_to = 'Corr', cols = 1:66)
  df_control$Treatment = "Control"
  df_stress = pivot_longer(df_stress, names_to = "Met", values_to = 'Corr', cols = 1:66)
  df_stress$Treatment = "Stress"
  df_long = rbind.data.frame(df_control, df_stress)
  df_long$Kernel = kernel
  df_long$Met = factor(df_long$Met, levels = met_names)
  df_long = df_long %>% arrange(CV, Treatment)
  
  return(df_long)
  }


############################################
DLs_RR_reader <- function(model, alpha){
  
  model_alphaL = list()
  for (i in 1:100){
    path = paste0("/Users/yebi/Library/CloudStorage/OneDrive-VirginiaTech/Research/Codes/research/RiceUNLMetabolites/GBLUP4Met/ARC_outputs/Met_same_genotypes_per_treatment/DLs_RR/", model, "_cv", i,"_alpha", alpha, ".csv")
    # print(path)
    model_rr = read.csv(path)
    model_alphaL[[i]] = model_rr
  }
  model_alpha = do.call(rbind, model_alphaL)
  
  return(model_alpha)
}

############################################
DLs_RR_mu_generator = function(model){
  DL_RR_1 = DLs_RR_reader(model, 1)
  DL_RR_100 = DLs_RR_reader(model, 100)
  DL_RR_1000 = DLs_RR_reader(model, 1000)
  DL_RR_10000 = DLs_RR_reader(model, 10000)
  
  PreCorr_DL_RR_1 = DL_extractor_func(df = DL_RR_1, met_names, paste0(model,"_RR_1"))
  PreCorr_DL_RR_100 = DL_extractor_func(df = DL_RR_100, met_names, paste0(model,"_RR_100"))
  PreCorr_DL_RR_1000 = DL_extractor_func(df = DL_RR_1000, met_names, paste0(model,"_RR_1000"))
  PreCorr_DL_RR_10000 = DL_extractor_func(df = DL_RR_10000, met_names, paste0(model,"_RR_10000"))
  
  DL_RR_1_mu <- pred_mean_generator(PreCorr_DL_RR_1)
  DL_RR_100_mu <- pred_mean_generator(PreCorr_DL_RR_100)
  DL_RR_1000_mu <- pred_mean_generator(PreCorr_DL_RR_1000)
  DL_RR_10000_mu <- pred_mean_generator(PreCorr_DL_RR_10000)
  DL_RR_comp = cbind(DL_RR_1_mu$corr.mean,
                        DL_RR_100_mu$corr.mean,
                        DL_RR_1000_mu$corr.mean,
                        DL_RR_10000_mu$corr.mean)
  
  DL_RR_mu = apply(DL_RR_comp, 1, max)
  DL_RR_mu = cbind.data.frame(DL_RR_1_mu[,1:2],
                                 corr.mean = DL_RR_mu)
  
  return(DL_RR_mu)
}


#############################################
DLs_RR_PreCorr_generator = function(model){
  DL_RR_1 = DLs_RR_reader(model, 1)
  DL_RR_100 = DLs_RR_reader(model, 100)
  DL_RR_1000 = DLs_RR_reader(model, 1000)
  DL_RR_10000 = DLs_RR_reader(model, 10000)
  
  PreCorr_DL_RR_1 = DL_extractor_func(df = DL_RR_1, met_names, paste0(model,"_RR_1"))
  PreCorr_DL_RR_100 = DL_extractor_func(df = DL_RR_100, met_names, paste0(model,"_RR_100"))
  PreCorr_DL_RR_1000 = DL_extractor_func(df = DL_RR_1000, met_names, paste0(model,"_RR_1000"))
  PreCorr_DL_RR_10000 = DL_extractor_func(df = DL_RR_10000, met_names, paste0(model,"_RR_10000"))
  
  comp = cbind(PreCorr_DL_RR_1$Corr,
               PreCorr_DL_RR_100$Corr,
               PreCorr_DL_RR_1000$Corr,
               PreCorr_DL_RR_10000$Corr)
  maxx = apply(comp, 1, max)
  PreCorr_DL = data.frame(CV = PreCorr_DL_RR_1$CV,
                          Met = PreCorr_DL_RR_1$Met,
                          Corr = comp,
                          Treatment = PreCorr_DL_RR_1$Treatment,
                          Kernel = model
                          )
  return(PreCorr_DL)
  
  }
