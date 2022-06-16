#!/bin/bash
#SBATCH --job-name=step3
#SBATCH --nodes=4
#SBATCH --time=6:00:00
#SBATCH --constraint=haswell
#SBATCH --qos=regular
#SBATCH -e step3.err
#SBATCH -o step3.out

source ./setup.sh
cd /global/cfs/cdirs/desc-fsi/ssi_vs_resim/submit/reprocess/
bps submit /global/homes/j/jemeyers/src/ssi_vs_resim/reprocess/03_bps_drp_step3.yaml
