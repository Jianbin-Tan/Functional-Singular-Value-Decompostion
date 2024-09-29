# Packages and functions
library(snowfall)
source("Function.R")
source("Jame_FC.R")

# Norm of functions
norm <- rowMeans(sapply(1:200, function(k){
  data_gen_clu(basis_num, n, obs_point, rat = NA, norm = NA)$norm
}))

# Simulation
cl <- parallel:: makeCluster(parallel::detectCores() - 1)
invisible(gc())
doParallel:: registerDoParallel(cl)

Result <- foreach(seed = seq(1, 10000, length.out = 2)) %dopar% {
  sim_func_clu(basis_num = basis_num, n = n,
               obs_point = obs_point, norm = norm,
               rat = rat, seed = seed)
}


save(Result, file = paste0("Result/Result_", n, "_", obs_point, "_", "clu", ".rda"), version = 2)
