## $1 = job name
## $2 = cv
## $3 = modell
## $4 = time

echo "#!/bin/bash

#SBATCH --mail-type=end
#SBATCH --mail-user=yebi@vt.edu
#SBATCH --job-name=\"$1\"
#SBATCH -N 2
#SBATCH -n 8
#SBATCH --gres=gpu:1
#SBATCH -t $4
#SBATCH -p a100_normal_q
#SBATCH --mem-per-cpu=32G
#SBATCH -A multiomicquantgen    #### <------- Change me

module load apps  site/tinkercliffs/easybuild/setup
module load Anaconda3/2020.11
source activate tf_gpu3
python FeatureExtractorMulti.py --cv $2 --modell $3

exit;

"
