#!/bin/bash
#SBATCH --job-name=coadd-inject
#SBATCH --nodes=6
#SBATCH --time=6:00:00
#SBATCH --constraint=haswell
#SBATCH --qos=regular
#SBATCH -e coadd-inject.err
#SBATCH -o coadd-inject.out

source ./setup.sh
cd /global/cfs/cdirs/desc-fsi/ssi_vs_resim/submit/reprocess/
bps submit /global/homes/j/jemeyers/src/ssi_vs_resim/reprocess/05_bps_coadd-inject.yaml
