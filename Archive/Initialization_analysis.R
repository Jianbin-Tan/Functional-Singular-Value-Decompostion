# Set address
rm(list=ls())
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Load functions and packages
source("Function.R")

## Data generation
data_gen_fun_init <- function(basis_num, n, obs_poi, rat, norm){
  
  sigma <- 2 * exp((basis_num):1 / 2)
  time_grid <- seq(0, 1, length.out = 101)
  basis <- fourier(time_grid, nbasis = basis_num + 3)[,2:(basis_num + 1)]
  
  a <- t(sapply(1:n, function(i){
    sapply(1:basis_num, function(k){
      rnorm(1)
    })
  }))
  
  a[,1] <- a[,1] / sqrt(sum(a[,1] ^ 2))
  for(k in 2:basis_num){
    a[,k] <- a[,k] - a[,1:(k-1)] %*% (t(a[,1:(k-1)]) %*% a[,k]) 
    a[,k] <- a[,k] / sqrt(sum(a[,k] ^ 2))
  }
  
  a <- t(t(a) * sigma) 
  
  fda_full <- basis %*% t(a)  
  
  if(is.numeric(rat) == T){
    obs_sig <- sqrt(norm * rat)
  }else{
    norm <- (colSums(fda_full ^ 2)) * 0.01
    obs_sig <- rep(1, n)
  }
  
  obs_dat <- lapply(1:n, function(i){
    time_num <- sample((obs_poi - 2):(obs_poi + 2), 1)
    time <- sort(sample(time_grid, time_num))
    time_mark <- sapply(1:time_num, function(k) which.min(abs(time_grid - time[k])))
    obs_point <- fda_full[time_mark,i] + rnorm(time_num, 0, obs_sig[i])
    return(list(time = time, obs = obs_point))
  })
  
  return(list(obs_dat = obs_dat,
              fda_full = fda_full,
              basis = basis,
              a = a,
              basis_num = basis_num,
              n = n,
              obs_point = obs_poi,
              rat = rat,
              time_grid = time_grid,
              obs_sig = obs_sig,
              norm = norm
  ))
}

## FSVD with random initialization
FSVD_rk_random <- function(Ly, dat_t, lambda, time_grid, init_num, abs){
  
  n <- length(Ly)
  Y <- unlist(Ly)
  
  fit <- lapply(1:30, function(g){
    
    a <- rnorm(n)
    a <- a / sqrt(sum(a ^ 2))
    
    a_coef <- sapply(1:length(dat_t$sam_mark), function(k) a[dat_t$sam_mark[k]])
    
    fit <- smooth.spline(x = dat_t$time_sam_tol, y = Y / a_coef, w = a_coef ^ 2 / dat_t$sam_num, 
                         all.knots = T, lambda = lambda * sum(a ^ 2), cv = F)
    
    phi <- predict(fit, x = time_grid)$y
    phi_t <- phi
    phi <- phi_t + 100
    iter <- 1
    
    pen <- 0
    
    while((sum(abs(phi - phi_t) ^ 2) / sum(abs(phi) ^ 2) > abs) & (iter < init_num)){
      iter <- 1 + iter
      phi <- phi_t
      
      a <- sapply(1:n, function(i){
        mean(Ly[[i]] * phi[dat_t$time_mark[[i]]]) / (mean(phi[dat_t$time_mark[[i]]] ^ 2) + pen * lambda)
      })
      a <- a / sqrt(sum(a ^ 2))
      
      a_coef <- sapply(1:length(dat_t$sam_mark), function(k) a[dat_t$sam_mark[k]])
      
      fit <- smooth.spline(x = dat_t$time_sam_tol, y = Y / a_coef, w = a_coef ^ 2 / dat_t$sam_num, 
                           all.knots = T, lambda = lambda * sum(a ^ 2), cv = F)
      
      phi_t <- predict(fit, x = time_grid)$y
      phi_t_2 <- predict(fit, x = time_grid, deriv = 2)$y
      
      pen <- sum((phi_t_2) ^ 2 * 0.01)
      
      na_mark <- sum(is.na(phi_t)) != 0
      if(na_mark){
        phi_t <- phi
      }
    } 
    
    phi <- phi_t
    
    a <- sapply(1:n, function(i){
      mean(Ly[[i]] * phi[dat_t$time_mark[[i]]]) / (mean(phi[dat_t$time_mark[[i]]] ^ 2) + pen * lambda)
    })
    a <- a / sqrt(sum(a ^ 2))
    
    rho <- sqrt(sum(phi ^ 2) * 0.01)
    phi <- phi / rho
    
    res <- lapply(1:n, function(i){
      Ly[[i]] - rho * a[i] * phi[dat_t$time_mark[[i]]]
    })
    
    loss <- sum(sapply(1:n, function(i){
      mean((Ly[[i]] - rho * a[i] * phi[dat_t$time_mark[[i]]]) ^ 2)
    })) + pen * lambda
    
    return(list(a = a, rho = rho, phi = phi, res = res, loss = loss))
  })
  
  mark <- which.min(sapply(1:length(fit), function(k) fit[[k]]$loss))
  
  
  return(fit[[mark]])
}

FSVD_CV_random <- function(Ly, dat_t, lambda, tran_datset, time_grid){
  
  n <- length(Ly)
  fit_p <- FSVD_rk_random(Ly, dat_t, lambda, time_grid, init_num = 500, abs = 10 ^ (-5))
  
  res <- sapply(1:length(tran_datset), function(k){
    fit <- FSVD_rk(Ly = tran_datset[[k]]$Ly, dat_t = tran_datset[[k]]$dat_t, 
                   phi = fit_p$phi, lambda, time_grid, init_num = 50, abs = 10 ^ (-3))
    pre_res <- sapply(1:n, function(i){
      if(length(tran_datset[[k]]$mark[[i]]) < length(Ly[[i]])){
        mark <- setdiff(1:length(Ly[[i]]), tran_datset[[k]]$mark[[i]])
        fit_dat <- fit$phi[dat_t$time_mark[[i]][mark]] * fit$rho * fit$a[i]
        mean((fit_dat - Ly[[i]][mark]) ^ 2)
      }else{
        0
      } 
    })
    return(mean(pre_res))
  })
  
  return(sum(res))
}

FSVD_tune_random <- function(Ly, time_grid, dat_t, tran_datset){
  
  m <- mean(dat_t$sam_num)
  n <- length(Ly)
  
  tran_datset <- lapply(1:length(tran_datset), function(k){
    Ly_i <- lapply(1:n, function(i) Ly[[i]][tran_datset[[k]]$mark[[i]]])
    return(list(Lt = tran_datset[[k]]$Lt, Ly = Ly_i, dat_t = tran_datset[[k]]$dat_t, mark = tran_datset[[k]]$mark))
  })
  
  time_tun <- exp(seq(log(10^(-6)), log(10 ^ (1)), length.out = 20)) 
  
  loss <- lapply(time_tun, function(lambda){
    FSVD_CV_random(Ly, dat_t, lambda, tran_datset, time_grid)
  })
  
  mark <- which.min(sapply(1:length(time_tun), function(i) loss[[i]]))
  lambda <- time_tun[mark]
  
  return(lambda)
}

FSVD_random <- function(Ly, Lt, R_max, R_pre, num_sel,
                 time_grid = seq(0, 1, length.out = 101),
                 Large_data = F
){
  
  n <- length(Ly)
  time_grid_mat <- unique(unlist(Lt))
  
  dat_t <- tran_dat(Lt, time_grid)
  m <- mean(sapply(1:n, function(i) length(Ly[[i]])))
  
  # Condition check
  if (length(time_grid_mat) < 5) {
    warning("Warning: Total number of the unique time points is less than 5.")
  }else if (sum(sapply(1:n, function(i) length(Ly[[i]])) < 2) > 0){
    warning("Warning: The number of non-zero observations in Ly[[i]] is less than 2 for some subjects.")
  }else{
    
    tran_datset <- lapply(1:5, function(k){
      mark <- lapply(1:n, function(i){
        if(k <= length(Lt[[i]])){
          setdiff(1:length(Lt[[i]]), seq(k, length(Lt[[i]]), 5))
        }else{
          1:length(Lt[[i]])
        }
      })
      Lt_i <- lapply(1:n, function(i) Lt[[i]][mark[[i]]])
      dat_t_i <- tran_dat(Lt_i, time_grid)
      return(list(Lt = Lt_i, dat_t = dat_t_i, mark = mark))
    })
    
    time_grid_mat_mark <- lapply(1:n, function(i){
      sapply(1:length(Lt[[i]]), function(k) which.min(abs(Lt[[i]][k] - time_grid_mat)))
    })
    time_grid_mat_tol_mark <- sapply(1:length(time_grid_mat), function(k) which.min(abs(time_grid_mat[k] - time_grid)))
    dat_raw <- sapply(1:n, function(i){
      A <- rep(NA, length(time_grid_mat))
      A[time_grid_mat_mark[[i]]] <- Ly[[i]]
      return(A)
    })
    
    ## R = 1
    lambda <- FSVD_tune_random(Ly, time_grid, dat_t, tran_datset)
    fit_FSVD <- FSVD_rk_random(Ly, dat_t, lambda = lambda, time_grid, init_num = 500, abs = 10 ^ (-5))
    
    R <- 1
    
    Rho <- fit_FSVD$rho
    Phi <- matrix(fit_FSVD$phi, length(time_grid))
    A <- matrix(fit_FSVD$a, n)
    
    ## R = 2
    lambda <- FSVD_tune_random(fit_FSVD$res, time_grid, dat_t, tran_datset)
    fit_FSVD <- FSVD_rk_random(fit_FSVD$res, dat_t, lambda, time_grid, init_num = 500, abs = 10 ^ (-5))
    
    while(R < R_max){
      R <- R + 1
      
      Rho <- c(Rho, fit_FSVD$rho)
      Phi <- cbind(Phi, fit_FSVD$phi)
      A <- cbind(A, fit_FSVD$a)
      
      lambda <- FSVD_tune_random(fit_FSVD$res, time_grid, dat_t, tran_datset)
      fit_FSVD <- FSVD_rk_random(fit_FSVD$res, dat_t, lambda, time_grid, init_num = 500, abs = 10 ^ (-5))
      # print(R)
    }
    
    R_max <- max(min(R_max, sum(cumsum(Rho[-R_max] <= Rho[-1] * 0.95) == 0) + 1), 2)
    
    if(num_sel == "FM"){
      
      IC <- sapply(1:R_max, function(r){
        log(1 / length(unlist(Ly)) * sum(sapply(1:n, function(i) sum((Ly[[i]] - matrix(Phi[dat_t$time_mark[[i]],1:r], nrow = length(dat_t$time_mark[[i]])) %*% diag(Rho[1:r], nrow = r) %*% c(A[i,1:r])) ^ 2)))) + r * log(min(n, m)) / min(n, m)
      })
      
      if(is.numeric(R_pre) == T){
        R <- R_pre
      }else{
        R <- which.min(IC)
        R <- max(R, 2)
      }
      
      W <- Phi[,1:R]
      fit <- qr(W)
      Fac_serial <- qr.Q(fit) * 10
      
      Loading <- A[,1:R] %*% diag(Rho[1:R], nrow = R) %*% t(qr.R(fit) / 10)
      
      fit <- qr(Loading)
      Loading <- qr.Q(fit)
      
      Fac_serial <- Fac_serial %*% t(qr.R(fit))
      
      return(list(
        Loading = Loading,
        Fac_serial = Fac_serial,
        Rho = Rho,
        A = A,
        Phi = Phi, 
        R = R,
        dat_t  = dat_t
      ))
      
    }else if(num_sel == "FD"){
      
      IC <- sapply(1:R_max, function(r){
        sum(sapply(1:n, function(i) length(Ly[[i]]) * log(mean((Ly[[i]] - matrix(Phi[dat_t$time_mark[[i]],1:r], nrow = length(dat_t$time_mark[[i]])) %*% diag(Rho[1:r], nrow = r) %*% c(A[i,1:r])) ^ 2)))) + 2 * n * r
      })
      
      if(is.numeric(R_pre) == T){
        R <- R_pre
      }else{
        R <- which.min(IC)
        R <- max(R, 2)
      }
      
      W <- Phi[,1:R]
      fit <- qr(W)
      Intric_basis <- qr.Q(fit) * 10
      
      Score <- A[,1:R] %*% diag(Rho[1:R], nrow = R) %*% t(qr.R(fit) / 10)
      
      W <- Phi
      fit <- qr(W)
      Intric_basis <- qr.Q(fit) * 10
      
      return(list(
        Intric_basis = Intric_basis,
        Score = Score,
        Rho = Rho,
        A = A,
        Phi = Phi, 
        R = R,
        dat_t  = dat_t
      ))
      
    }else{
      
      IC <- Rho[-R_max] / Rho[-1] 
      
      if(is.numeric(R_pre) == T){
        R <- R_pre
      }else{
        R <- which.max(IC)
      }
      
      return(list(
        Rho = Rho,
        A = A,
        Phi = Phi, 
        R = R,
        dat_t  = dat_t
      ))
    }
  }
}

# Implementation
basis_num <- 3 # Number of basis functions used
rat <- 0.05 # Noise level

sample_mark <- c(50, 100, 150) # Sample size
point_mark <- c(6, 8, 10) # Mean number of observed time points

imp_fun <- function(seed, basis_num, n, obs_point, rat, norm){
  
  set.seed(seed)
  dat_col <- data_gen_fun_init(basis_num, n, obs_point, rat, norm)
  dat <- dat_col$obs_dat
  
  Ly <- lapply(1:n, function(i) dat[[i]]$obs)
  Lt <- lapply(1:n, function(i) dat[[i]]$time)
  time_grid <- seq(0, 1, 0.01)
  
  fit_1 <- FSVD(Ly, Lt, R_max = 3, R_pre = 3, num_sel = "SVD")
  fit_2 <- FSVD_random(Ly, Lt, R_max = 3, R_pre = 3, num_sel = "SVD")
  
  error_1 <- sapply(1:basis_num, function(k) sqrt(1 - sum(fit_1$Phi[,k] * dat_col$basis[,k]) ^ 2 / (sum(fit_1$Phi[,k] ^ 2)) / (sum(dat_col$basis[,k] ^ 2))))
  error_2 <- sapply(1:basis_num, function(k) sqrt(1 - sum(fit_2$Phi[,k] * dat_col$basis[,k]) ^ 2 / (sum(fit_2$Phi[,k] ^ 2)) / (sum(dat_col$basis[,k] ^ 2))))
  
  return(data.frame(error_1, error_2))
}

# Implementation
library(snowfall)

# Simulation
for(n in sample_mark){
  for(obs_point in point_mark){
    norm <- rowMeans(sapply(1:200, function(k){
      data_gen_fun_init(basis_num, n, obs_point, rat = NA, norm = NA)$norm
    }))
    
    sfInit(parallel = T, cpus = 40)
    sfExport("basis_num", "n", "obs_point", "rat", "norm")
    sfExport("imp_fun", "FSVD_rk_random", "FSVD_CV_random", "FSVD_tune_random", "FSVD_random", "data_gen_fun_init")
    sfSource("Function.R")
    Result <- sfLapply(seq(1, 10000, length.out = 100), imp_fun,
                       basis_num = basis_num, n = n,
                       obs_point = obs_point, norm = norm,
                       rat = rat)
    
    save(Result, file = paste0("Result/Result_init_analysis_", n, "_", obs_point, ".rda"), version = 2)
    sfStop()
  }
}

# Plot
data_plot <- NULL
for(n in sample_mark){
  for(obs_point in point_mark){
    load(paste0("Result/Result_init_analysis_", n, "_", obs_point, ".rda"))
    
    for(k in 1:3){
      result <- sapply(1:length(Result), function(i) Result[[i]][k,1])
      
      data_plot <- rbind(data_plot, data.frame(
        value = c(sapply(1:length(Result), function(i) Result[[i]][k,1]), sapply(1:length(Result), function(i) Result[[i]][k,2])),
        method = c(rep("Completion", 100), rep("Random", 100)),
        comp = rep(paste0("Comp ", k), 200),
        n = rep(n, 200),
        obs = rep(paste0("{", obs_point - 2, ",...,", obs_point + 2, "}"), 200)
      ))
    }
  }
}

p_1 <- ggplot(data_plot[data_plot$comp == "Comp 1",]) +
  geom_boxplot(
    aes(x = method, y = value),
    width = 0.6, size = 0.4, outliers = F, alpha = 0.85
  ) +
  facet_wrap(~n + obs, ncol = 3) +
  labs(x = "", y = "Dist", title = "(A) Component 1") +
  theme_bw() +
  theme(text=element_text(size=15),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 0, size = 9),
        panel.border = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5)) 

p_2 <- ggplot(data_plot[data_plot$comp == "Comp 2",]) +
  geom_boxplot(
    aes(x = method, y = value),
    width = 0.6, size = 0.4, outliers = F, alpha = 0.85
  ) +
  facet_wrap(~n + obs, ncol = 3) +
  labs(x = "", y = "Dist", title = "(B) Component 2") +
  theme_bw() +
  theme(text=element_text(size=15),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 0, size = 9),
        panel.border = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5)) 

p_3 <- ggplot(data_plot[data_plot$comp == "Comp 3",]) +
  geom_boxplot(
    aes(x = method, y = value),
    width = 0.6, size = 0.4, outliers = F, alpha = 0.85
  ) +
  facet_wrap(~n + obs, ncol = 3) +
  labs(x = "", y = "Dist", title = "(C) Component 3") +
  theme_bw() +
  theme(text=element_text(size=15),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 0, size = 9),
        panel.border = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5)) 

p <- gridExtra::arrangeGrob(p_1, p_2, p_3, nrow = 1)

ggsave(paste0("Figure/", "init_analysis", ".pdf"), p, width = 16, height = 6, dpi = 300)
