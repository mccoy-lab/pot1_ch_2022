#!/bin/bash

#======================================================#
# Submit 10,000 simulation runs as an sbatch array job #
#======================================================#

#SBATCH -N 1
#SBATCH --time=15:00:00
#SBATCH --partition=defq
#SBATCH --mem=10G
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=1-10000%100

SLIM=~/code/SLiM-3.7.1/bin/slim

# simulation script
simfile="hpc_simulation.slim"

# make directory for output files
outdir=${simtype}_${telo_mean}
outfile=$outdir/${simtype}_${telo_mean}_${SLURM_ARRAY_TASK_ID}.out
mkdir -p $outdir

# run simulation in SLiM; telo_mean and shorten_mean
# are specified in run_simulation_wrapper.sh
$SLIM \
	-d gen_time_days=280 \
	-d telo_mean=$telo_mean \
	-d telo_sd=1250 \
	-d shorten_mean=$shorten_mean \
	-d "shorten_sd=70/6" \
	-d death_thresh=4000 \
	-d pAsym=0.9 \
	-d pAsym_mut=0.8 \
	-d pSym=0.525 \
	-d pSym_mut=0.7 \
	-d mut_mult=1.1 \
	$simfile \
	> $outfile
