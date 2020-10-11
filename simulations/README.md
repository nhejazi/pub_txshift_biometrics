To reproduce the results of the simulation experiments reported in the main
manuscript and the supporting information, examine the scripts in the `R`
directory. Note that the simulation experiments are both time-intensive and
computationally intensive; however, the simulation experiments may be re-run
using the scripts in the `R` subdirectory. There were two general sets of
experiments: (1) those assessing estimator performance across differences in
strategies for estimation of the sampling mechanism, and (2) those assessing
performance of our proposed estimation strategies in the context of simulated
data based on the HVTN 505 trial.

## Performance Experiments

For the simulation experiments described in the main manuscript (Section 4) and
Section 1b of the Supporting Information, the directory `R/main_performance`
contains scripts to generate data for and extract results from different
estimator variants, while the script `R/03a_main_performance_sims.R` can be used
to generate plots from the results. All other scripts contain auxiliary
functions for setting up the simulation study. As the simulation experiments
reported were run circa July 2019 on an HPC cluster, scripts may contain some
hardcoded paths that will need to be altered to re-run the experiments; the data
produced and analyzed is made available in the directory `data/main`.

## HVTN 505 Experiments

For the simulation experiments described in Sections 2a and 2b of the Supporting
Information, the directory `R/supp_hvtn505` contains two scripts to generate the
results for the scenario with a moderate effect of exposure and no effect of
exposure, respectively. Analogous to the above, scripts to summarize estimator
performance in plots are `R/03b_hvtn505_null_effect_sims.R` and
`R/03c_hvtn505_moderate_effect_sims.R`. Scripts providing auxiliary support for
setting up the simulation study are shared with the scenario described above.
These simulation experiments were most recently run circa October 2019, and the
corresponding scripts may contain hardcoded paths specific to the HPC cluster on
which they were run; the data produced and analyzed is made available in the
directories `data/supp_2a` and `data/supp_2b`.

---

At that time, R version 3.5.1 was used
on [UC Berkeley's Savio high-performance computing
cluster](https://research-it.berkeley.edu/services/high-performance-computing/system-overview).
