# Role of genomics on regulating rice grain metabolic variability under warmer nights: A statistical and image-based deep learning approach

## Abstract
It has been argued that metabolites can be used to accelerate crop improvement because metabolic profiles in crops are generally under genetic control. Evaluating the role of genetics in metabolic variation is a longstanding challenge. Rice, one of the world's most important staple crops, is known to be sensitive to recent increases in nighttime temperatures. Quantification of metabolic levels can help measure rice responses to high nighttime temperature (HNT) stress. However, the extent of metabolic variation that can be explained by regression on whole-genome molecular markers remains to be answered. In the current study, primary metabolites of a rice diversity panel generated from grains using gas chromatography-mass spectrometry were used. The metabolites obtained were low to moderately heritable, and the genomic prediction accuracies of the metabolites were within the expected upper limit set by their genomic heritability estimates. Genomic heritability estimates were slightly higher in the control group than in the HNT group. Genomic correlation estimates for the same metabolites between the control and HNT conditions indicated the presence of genotype by environment interactions. Reproducing kernel Hilbert spaces regression and deep learning, which treat markers as images, improved prediction accuracy, suggesting that some metabolites are under non-additive genetic control. Joint analysis of multiple metabolites simultaneously was effective in improving prediction accuracy by exploiting correlations among metabolites. The current study serves as an important first step in evaluating the cumulative effects of the genome in regulating metabolic variation under control and HNT conditions. 

Preprint: link

## 0. Data Preprocessing
- [.Rmd file](https://github.com/yebigithub/GBLUP4Met/blob/main/PreProcessing/DataPreprocessing.Rmd) Including metabolite and genotype data cleaning
## 1. Genomic heritability of metabolites
- [.R file](https://github.com/yebigithub/GBLUP4Met/blob/main/heritability/h2_calculate.R) Using sommer package to calculate heritability for metabolites.
- [.Rmd file](https://github.com/yebigithub/GBLUP4Met/blob/main/heritability/h2_plot_drawing.Rmd) Drawing heritability plots.

<p align="center">
  <img src='https://github.com/yebigithub/GBLUP4Met/blob/main/heritability/Figure2.png' width='70%' height='70%' alt="Figure 2: Genomic heritability estimates of metabolite accumulation in control and high night temperature stress conditions. A) Scatter plot. B) Density plot. Solid and dashed lines indicate mean and median, respectively. C) Agreement of heritability estimates between control and high night temperature stress conditions. Metabolites in green and red colors indicate that the heritability difference between control and high night temperature stress conditions was small (< 0.05) and large (> 0.1)">
</p>

## 2. Single-trait genomic prediction of metabolites
- [.R file](https://github.com/yebigithub/GBLUP4Met/blob/main/SingleTraitGBLUP/GBLUP_whole.R) Running Single trait GBLUP in cluster.
- [.Rmd file](https://github.com/yebigithub/GBLUP4Met/blob/main/SingleTraitGBLUP/SingleTraitGBLUP.Rmd) Drawing Single trait GBLUP plots.
- [.Rmd file](https://github.com/yebigithub/GBLUP4Met/blob/main/SingleTraitGBLUP/GK.Rmd) Selecting suitable bandwidth for RKHS.
- [.R file](https://github.com/yebigithub/GBLUP4Met/blob/main/SingleTraitGBLUP/GK.R) Runing Single trait RKHS in cluster.
## 3. Genomic correlation between the same metabolite in different treatments
## 4. Exporatory factor analysis
## 5. Simultaneous regression modeling of metabolites
## 6. Deep learning models