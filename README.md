# Genomic prediction of metabolic content in rice grain in response to warmer night conditions

Bi, Y., Walia, H., Obata, T., & Morota, G. (2024). Genomic prediction of metabolic content in rice grain in response to warmer night conditions. Crop Science, e21435. https://doi.org/10.1002/csc2.21435

## Abstract
It has been argued that metabolic content can be used as a selection marker to accelerate crop improvement because metabolic profiles in crops are often under genetic control. Evaluating the role of genetics in metabolic variation is a long-standing challenge. Rice, one of the world's most important staple crops, is known to be sensitive to recent increases in nighttime temperatures. Quantification of metabolic levels can help measure rice responses to high night temperature (HNT) stress. However, the extent of metabolic variation that can be explained by regression on whole-genome molecular markers remains to be evaluated. In the current study, we generated metabolic profiles for mature grains from a subset of rice diversity panel accessions grown under optimal and HNT conditions. Metabolite accumulation was low to moderately heritable, and genomic prediction accuracies of metabolite accumulation were within the expected upper limit set by their genomic heritability estimates. Genomic heritability estimates were slightly higher in the control group than in the HNT group. Genomic correlation estimates for the same metabolite accumulation between the control and HNT conditions indicated the presence of genotype-by-environment interactions. Reproducing kernel Hilbert spaces regression and image-based deep learning improved prediction accuracy, suggesting that some metabolite levels are under non-additive genetic control. Joint analysis of multiple metabolite accumulation simultaneously was effective in improving prediction accuracy by exploiting correlations among metabolites. The current study serves as an important first step in evaluating the cumulative effect of markers in influencing metabolic variation under control and HNT conditions. 


## 0. Data Preprocessing
- [.Rmd file](https://github.com/yebigithub/GBLUP4Met/blob/main/PreProcessing/DataPreprocessing.Rmd) Including metabolite and genotype data cleaning
## 1. Genomic heritability of metabolites
- [.R file](https://github.com/yebigithub/GBLUP4Met/blob/main/heritability/h2_calculate.R) Using sommer package to calculate heritability for metabolites.
- [.Rmd file](https://github.com/yebigithub/GBLUP4Met/blob/main/heritability/h2_plot_drawing.Rmd) Drawing heritability plots.

<p align="center">
  <img src='./heritability/Figure2.png' width='70%' height='70%' alt="">
</p>
Figure 2: Genomic heritability estimates of metabolite accumulation in control and high night temperature stress conditions. A) Scatter plot. B) Density plot. Solid and dashed lines indicate mean and median, respectively. C) Agreement of heritability estimates between control and high night temperature stress conditions. Metabolites in green and red colors indicate that the heritability difference between control and high night temperature stress conditions was small (< 0.05) and large (> 0.1)

## 2. Single-trait genomic prediction of metabolites
- [.R file](https://github.com/yebigithub/GBLUP4Met/blob/main/SingleTraitGBLUP/GBLUP_whole.R) Running Single trait GBLUP in cluster.
- [.Rmd file](https://github.com/yebigithub/GBLUP4Met/blob/main/SingleTraitGBLUP/SingleTraitGBLUP.Rmd) Drawing Single trait GBLUP plots.
- [.Rmd file](https://github.com/yebigithub/GBLUP4Met/blob/main/SingleTraitGBLUP/GK.Rmd) Selecting suitable bandwidth for RKHS.
- [.R file](https://github.com/yebigithub/GBLUP4Met/blob/main/SingleTraitGBLUP/GK.R) Runing Single trait RKHS in cluster.

<p align="center">
  <img src='./SingleTraitGBLUP/Figure4.png' width='70%' height='70%' alt="">
</p>
Figure 4: Genomic prediction accuracy of metabolite accumulation in control and high night
temperature stress conditions. A) Box plot. The horizontal line indicates the mean value.
B) Density plot. The solid and dashed lines indicate the mean and median, respectively.
C) Agreement of genomic prediction accuracy between control and high night temperature
stress conditions. Metabolite accumulations in green and red colors indicate that the genomic
prediction difference between control and high night temperature stress conditions was small
(< 0.05) and large (> 0.1).

## 3. Genomic correlation between the same metabolite in different treatments
- [.R file](./GenoCorrSameMet/GenoCorrSameMet.R) Running multi-trait genomic correlation.
- [.Rmd file](./GenoCorrSameMet/GenoCorrSameMet.Rmd) Drawing multi-trait genomic correlation plots.

<p align="center">
  <img src='./GenoCorrSameMet/Figure3.png' width='70%' height='70%' alt="">
</p>  
Figure 3: Genomic correlation estimates between the same metabolite accumulation measured under control and high night temperature stress conditions.  A) Scatter plot. B) Bar chart. Solid and dashed lines indicate mean and median, respectively.

## 4. Exporatory factor analysis
- [.Rmd file](./FactorAnalysis/FA4Met.Rmd) Factorial analysis to identify underlying latent factors controlling metabolites.

## 5. Simultaneous regression modeling of metabolites
- [.R file](./SimultaneousRegression/MegaLLM.R) Running MegaLMM for genomic prediction.
- [.R file](./SimultaneousRegression/MegaLLM_GK.R) Running MegaLMM for RKHS.
- [.Rmd file](./SimultaneousRegression/MegaLMM_drawing.Rmd) Drawing barplot, density plots for MegaLMM genomic prediction model.
- [.Rmd file](./SimultaneousRegression/GenomicCorr.Rmd) Drawing genomic correlation density plot.

<p align="center">
  <img src='./SimultaneousRegression/ggcorr_MegaLMM_density_-1.png' width='70%' height='70%' alt="">
</p>
Figure 7: Genomic correlation estimates between different metabolite accumulation in control and high night temperature stress conditions. The solid and dashed lines indicate mean and median, respectively.
<br/><br/>

<p align="center">
  <img src='./SimultaneousRegression/G_GK_Mega-1.png' width='70%' height='70%' alt="">
</p> 
Figure 8: Percentage difference of gain in prediction accuracy for multi-trait genomic best linear unbiased prediction (MegaLMM-G) and multi-trait reproducing kernel Hilbert spaces regression (MegaLMM-GK) relative to single-trait genomic best linear unbiased prediction (A). Density plots of percentage difference are shown for MegaLMM-G (B) and MegaLMM-GK (C).


## 6. Deep learning models
- [.ipynb](./DL/ConvertSNP2Image.ipynb) Shows examples about how to convert SNP tabular data into SNP images.
- [.py file](./DL/ConvertSNP2Image.py) Loop converting for SNPs in all chromosomes.
- [.py file](./DL/FeatureExtractorMulti.py) Convolutional neural network with multiple branches.
- [.Rmd file](./DL/DL_drawing.Rmd) Drawing barplot to compare performance of all deep learning models and RKHS.
  
<p align="center">
  <img src='./DL/workflow4DL-1.png' width='80%' height='80%' alt="">
</p> 
Figure1: Flowchart of converting single nucleotide polymorphisms to image data

<br/><br/>
<p align="center">
  <img src='./DL/SNPimagesExample-1.png' width='70%' height='60%' alt="">
</p> 
Figure 5: Example of a set of single nucleotide polymorphisms transformed into image data for a randomly selected genotype. Images of 12 chromosomes were processed in the multi-channel convolutional neural networks

<br/><br/>
<p align="center">
  <img src='./DL/DL_dff_barplot_0.1_all.jpg' width='90%' height='90%' alt="">
</p> 
Figure 6: Percentage difference of gain in prediction accuracy for single-trait reproducing kernel Hilbert spaces regression (RKHS), VGG16, ResNet50 EfficientNetB7, InceptionV3, MobileNetV2, and DenseNet201 relative to single-trait genomic best linear unbiased prediction.

## 7. Supplementary
- [.Rmd file](./Supplementary/PhenoCorr.Rmd) Calculating phenotypical correaliton between metabolites in control and stress conditions.
- [.Rmd file](./SimultaneousRegression/MegaLMM_ggcorr.R) Drawing MegaLMM genomic correlation heatmaps.
- [.Rmd file](./FactorAnalysis/FA4Met.Rmd) Drawing factorial analysis heatmaps.
- [.Rmd file](./Supplementary/FA_drawing.Rmd) Drawing factorial analysis density plots.
