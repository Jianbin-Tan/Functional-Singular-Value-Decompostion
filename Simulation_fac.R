# Packages and functions
library(snowfall)
source("Function.R")

# Simulation
sfInit(parallel = T, cpus = 40)
sfExport("basis_num", "n", "obs_point", "rat")
sfSource("Function.R")
Result <- sfLapply(seq(1, 10000, length.out = 100), sim_func_fac,
                   basis_num = basis_num, n = n, 
                   obs_point = obs_point,
                   rat = rat)

sfStop()

save(Result, file = paste0("Result/Result_fac_", n, "_", obs_point, ".rda"), version = 2)

