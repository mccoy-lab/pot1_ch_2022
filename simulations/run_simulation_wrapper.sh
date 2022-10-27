#!/bin/bash

#======================================================#
# Run simulations for each telomere length background, #
# passing in background-specific parameters to sbatch  #
#======================================================#

# 1st percentile, short telomere
sbatch --job-name=wt_8kb \
       --export=ALL,simtype="wt",telo_mean=8600,shorten_mean=100 \
       run_simulation.sh

# 50th percentile, normal telomere
sbatch --job-name=wt_11kb \
       --export=ALL,simtype="wt",telo_mean=11000,shorten_mean=100 \
       run_simulation.sh

# 99th percentile, long telomere
sbatch --job-name=wt_13kb \
       --export=ALL,simtype="wt",telo_mean=13400,shorten_mean=100 \
       run_simulation.sh

# 99th percentile, POT1 mutation carrier
sbatch --job-name=pot1_13kb \
        --export=ALL,simtype="pot1",telo_mean=13400,shorten_mean=-50 \
        run_simulation.sh