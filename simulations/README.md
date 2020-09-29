To reproduce the results of the simulation experiments reported in the main
manuscript and the supporting information, examine the scripts in the `R`
directory. Note that the simulation experiments are both time-intensive and
computationally intensive; however, the simulation experiments may be re-run as
per the following.

prior to posting of
the pre-print and submission to JRSS-B.

## Performance Experiments

The simulation experiments reported run circa July 2019,

For the simulation experiments described in the main manuscript (section 4) and
1b of the supporting information:

```r
install.packages("here")
library(here)
source(here("R", "00_install_pkgs.R"))
source(here("R", "03_run_simulation_param.R"))
```

```r
install.packages("here")
library(here)
source(here("R", "00_install_pkgs.R"))
source(here("R", "03_run_simulation_param.R"))
```

## HVTN 505 Experiments

The simulation experiments in the suplement October 2019

For the simulation experiments described in section 2a and 2b of the supporting
information:

```r
install.packages("here")
library(here)
source(here("R", "00_install_pkgs.R"))
source(here("R", "03_run_simulation_param.R"))
```

```r
install.packages("here")
library(here)
source(here("R", "00_install_pkgs.R"))
source(here("R", "03_run_simulation_param.R"))
```

---

At that time, [R version 3.5.1]() was used
on [UC Berkeley's Savio high-performance computing
cluster](https://research-it.berkeley.edu/services/high-performance-computing/system-overview).
