#!/bin/bash
#SBATCH --job-name=injected-step1
#SBATCH --nodes=6
#SBATCH --time=4:00:00
#SBATCH --constraint=haswell
#SBATCH --qos=regular
#SBATCH -e injected-step1.err
#SBATCH -o injected-step1.out

source ./setup.sh
cd /global/cfs/cdirs/desc-fsi/ssi_vs_resim/submit/reprocess/
bps submit /global/homes/j/jemeyers/src/ssi_vs_resim/reprocess/06_bps_injected_visits_step1.yaml
