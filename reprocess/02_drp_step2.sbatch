#!/bin/bash
#SBATCH --job-name=step2
#SBATCH --nodes=4
#SBATCH --time=6:00:00
#SBATCH --constraint=haswell
#SBATCH --qos=regular
#SBATCH -e step2.err
#SBATCH -o step2.out

source ./setup.sh
cd /global/cfs/cdirs/desc-fsi/ssi_vs_resim/submit/reprocess/
bps submit /global/homes/j/jemeyers/src/ssi_vs_resim/reprocess/02_bps_drp_step2.yaml
