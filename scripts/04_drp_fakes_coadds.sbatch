#!/bin/bash
#SBATCH --job-name=ssi_coadd_4430_00
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --constraint=haswell
#SBATCH --qos=regular
#SBATCH -e ssi_coadd_4430_00.err
#SBATCH -o ssi_coadd_4430_00.out

source ./04_setup.sh
bps submit /global/homes/j/jemeyers/src/ssi_vs_resim/scripts/04_bps_drp_fakes_coadds.yaml
