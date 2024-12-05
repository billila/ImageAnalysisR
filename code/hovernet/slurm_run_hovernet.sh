#!/bin/sh 
# FILENAME: ila_hovernet

#SBATCH -A cis240854-gpu       # allocation name
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=4   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --gpus-per-node=4     # Number of GPUs per node
#SBATCH -J ila_hovernet          # Job name
#SBATCH -o myjob.o%j          # Name of stdout output file
#SBATCH -e myjob.e%j          # Name of stderr error file
#SBATCH -p gpu                # Queue (partition) name
#SBATCH --time=48:00:00
#SBATCH --mail-user=ilaria.billato@phd.unipd.it
#SBATCH --mail-type=all       # Send email to above address at begin and end of job

# Manage processing environment, load compilers, and applications.
module purge
module load modtree/gpu
module list

module load anaconda 
module list

conda activate
conda env list
conda activate /home/x-ibillato/.conda/envs/2024.02-py311/hovernet_new

conda env list

python /anvil/scratch/x-ibillato/hover_net/run_infer.py \
--gpu='0,1,2,3' \
--nr_types=6 \
--type_info_path=/anvil/scratch/x-ibillato/hover_net/type_info.json \
--batch_size=64 \
--model_mode=fast \
--model_path=/anvil/scratch/x-ibillato/hover_net/hovernet_fast_pannuke_type_tf2pytorch.tar \
--nr_inference_workers=16 \
--nr_post_proc_workers=16 \
wsi \
--input_dir=/anvil/scratch/x-ibillato/tcga_images/ \
--output_dir=/anvil/scratch/x-ibillato/output/ \
--save_thumb \
--save_mask

