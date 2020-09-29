# lie to pkgbuild, as per Jeremy
pkgbuild:::cache_set("has_compiler", TRUE)

# set user-specific package library
if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths("/global/scratch/nhejazi/R")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# from CRAN
install.packages(c("here", "remotes", "future", "future.apply", "doFuture",
                   "foreach", "tidyverse", "data.table", "Rsolnp", "nnls",
                   "Rcpp", "doRNG", "origami"),
                 lib = "/global/scratch/nhejazi/R")

# use remotes to install from GitHub
remotes::install_github(c("osofr/simcausal@master",
                          "osofr/condensier@master",
                          "tlverse/hal9001@master",
                          "nhejazi/haldensify@master",
                          "tlverse/sl3@v1.2.0",
                          "nhejazi/txshift@master"
                         ),
                        lib = "/global/scratch/nhejazi/R")

# update all packages
update.packages(lib.loc = "/global/scratch/nhejazi/R")
