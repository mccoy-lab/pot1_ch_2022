# Simulations of hematopoietic progenitor cells

## Scripts in this repo

* `hpc_simulation.slim`: SLiM code for simulation
* `run_simulation.sh`: bash script for running one simulation
* `run_simulation_wrapper.sh`: slurm wrapper script for parallelizing multiple simulations on a computing cluster

## Description of simulation

### Overview

We used simulations, constructed with the software SLiM ([Haller and Messer, 2019](https://academic.oup.com/mbe/article/36/3/632/5229931)), to test the role of telomere length (TL) on the longevity of a single heterozygous gain-of-function somatic mutation arising in one of 100,000 **hematopoietic progenitor cells (HPCs)** at birth.

We ran 10,000 simulations for each of four telomere length backgrounds - the mean clinically validated flow-FISH measured average TLs ([Alder et al., 2018](https://www.pnas.org/doi/10.1073/pnas.1720427115)) at:

1. 1st percentile (8.6 kb)
2. 50th percentile (11.0 kb)
3. 99th percentile (13.4 kb)
4. POT1 mutant background (13.4 kb)

### Setup

We begin each simulation with a population of 100,000 cells. We select one cell at birth to carry a CH driver mutation, which confers a higher probability of entering the cell cycle. We did not account for secondary mutations that may have been gained over an individual's lifetime.

We draw every cell’s telomere length from a normal distribution with mean μ (i.e., 8.6 kb, 11.0 kb, etc.) and standard deviation 1.25 kb. This standard deviation was determined from the distribution of telomere lengths within an individual at birth, as estimated in [Mitchell et al.](https://www.nature.com/articles/s41586-022-04786-y)

### Events

Within each day of the simulation, the same set of events occur:

#### 1. Cell fate decision

##### Division vs. non division

We set the HPC division rate at once every 280 days (daily division probability = 1/280), as described in [Mitchell et al.](https://www.nature.com/articles/s41586-022-04786-y)

Progenitor cells with the driver clonal mutation had a higher probability of entering the cell cycle (1.1x daily division probability), based on literature estimating the advantage of _JAK2_ and _DNMT3A_ hotspot mutations.

##### Fate decision

Probabilities of fate decisions were derived from [Lahouel et al.](https://www.pnas.org/doi/10.1073/pnas.1914589117) as follows:

* 0.9 asymmetric division
* 0.0525 symmetric division
* 0.0475 differentiation

In cells with driver mutations, the probabilities were:

* 0.8 asymmetric division
* 0.14 symmetric division
* 0.06 differentiation

#### 2. Telomere length change

All dividing cells experience a 100 bp/division telomere shortening rate, estimated from [Mitchell et al.](https://www.nature.com/articles/s41586-022-04786-y). We used a 50 bp/division gain in the _POT1_ mutant simulation, estimated from the data in Figure 4 of this paper.

#### 3. Cell death

Every cell with TL <4 kb experienced cell death with a probability dependent entirely on telomere length. Likelihood of death increased linearly from 0 (certain survival) to 1 (certain death) between 4 kb and 2 kb ([Alder et al., 2018](https://www.pnas.org/doi/10.1073/pnas.1720427115)).
