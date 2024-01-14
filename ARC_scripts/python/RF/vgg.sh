#!/bin/bash

#SBATCH --mail-type=end
#SBATCH --mail-user=yebi@vt.edu
#SBATCH --job-name="vgg"
#SBATCH -N 2
#SBATCH -n 8
#SBATCH --gres=gpu:1
#SBATCH -t 10:00:00
#SBATCH -p a100_normal_q
#SBATCH --mem-per-cpu=32G
#SBATCH -A multiomicquantgen    #### <------- Change me

module load apps  site/tinkercliffs/easybuild/setup
module load Anaconda3/2020.11
source activate tf_gpu3
python FeatureExtractor.py

exit;