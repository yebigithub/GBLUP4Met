#!/bin/bash

### multi.sh   ### <---------- Change me
###########################################################################
## environment & variable setup
####### job customization
#SBATCH --mail-type=end
#SBATCH --mail-user=yebi@vt.edu
#SBATCH --job-name="muti_trt"   ### <---------- Change me
#SBATCH -N 2
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -p normal_q
#SBATCH --mem-per-cpu=32G
#SBATCH -A multiomicquantgen    #### <------- change me
####### end of job customization
# end of environment & variable setup
###########################################################################
#### add modules on TC/Infer
module load containers/singularity/apptainer-wrapper
### from DT/CA, use module load singularity
module list
#end of add modules
###########################################################################
###print script to keep a record of what is done
cat GBLUP_con1.sh  ### <---------- Change me
###########################################################################
echo start running R
## note, on DT/CA, you should replace projects with groups

singularity exec --bind=/work,/projects \
    /projects/arcsingularity/ood-rstudio141717-bio_4.1.0.sif Rscript MultiTrait.R  ### <---------- Change me

exit;
