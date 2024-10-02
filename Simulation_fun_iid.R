# Packages and functions
library(snowfall)
source("Function.R")

# Norm of functions
norm <- rowMeans(sapply(1:200, function(k){
  data_gen_fun_iid(basis_num, n, obs_point, rat = NA, norm = NA)$norm
}))

# Simulation
sfInit(parallel = T, cpus = 50)
sfExport("basis_num", "n", "obs_point", "rat", "norm")
sfSource("Function.R")
Result <- sfLapply(seq(1, 10000, length.out = 100), sim_func_iid,
                   basis_num = basis_num, n = n,
                   obs_point = obs_point, norm = norm,
                   rat = rat)

sfStop()

save(Result, file = paste0("Result/Result_iid_", n, "_", obs_point, ".rda"), version = 2)
