#!/bin/bash
#SBATCH --job-name=step1
#SBATCH --nodes=4
#SBATCH --time=6:00:00
#SBATCH --constraint=haswell
#SBATCH --qos=regular
#SBATCH -e step1.err
#SBATCH -o step1.out

source ./setup.sh
cd /global/cfs/cdirs/desc-fsi/ssi_vs_resim/submit/reprocess/
bps submit /global/homes/j/jemeyers/src/ssi_vs_resim/reprocess/01_bps_drp_step1.yaml
