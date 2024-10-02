# Package
library(mclust, quietly=TRUE)
library(mvtnorm)
library(ggplot2)
library(mgcv)
library(fda)
library(fdapace)
library(xtable)
library(MASS)
library(Matrix)
library(HDTSA)
library(aricode)
library(fdaMocca)
library(glmnet)
library(fields)
library(foreach)
library(deSolve)
library(filling)

################################################################################
# Data generation
## Intrinsic basis functions' case
### Heterogeneous cases
data_gen_fun <- function(basis_num, n, obs_poi, rat, norm){
  
  sigma <- 2 * exp((basis_num):1 / 2)
  time_grid <- seq(0, 1, length.out = 101)
  basis <- fourier(time_grid, nbasis = basis_num + 3)[,2:(basis_num + 1)]
  
  a <- t(sapply(1:n, function(i){
    sapply(1:basis_num, function(k){
      sin(k * pi * (i + (n/ 4)) / (2 * n))
    })
  }))
  
  a[,1] <- a[,1] / sqrt(sum(a[,1] ^ 2))
  for(k in 2:basis_num){
    a[,k] <- a[,k] - a[,1:(k-1)] %*% (t(a[,1:(k-1)]) %*% a[,k]) 
    a[,k] <- a[,k] / sqrt(sum(a[,k] ^ 2))
  }
  
  a <- t(sapply(1:n, function(i){
    sapply(1:basis_num, function(k){
      rnorm(1, a[i,k],  abs(a[i,k]))
    })
  }))
  
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

### i.i.d. cases
data_gen_fun_iid <- function(basis_num, n, obs_poi, rat, norm){
  
  sigma <- 2 * exp((basis_num):1 / 2)
  time_grid <- seq(0, 1, length.out = 101)
  basis <- fourier(time_grid, nbasis = basis_num + 3)[,2:(basis_num + 1)]
  
  a <- t(sapply(1:n, function(i){
    sapply(1:basis_num, function(k){
      0
    })
  }))
  
  # a[,1] <- a[,1] / sqrt(sum(a[,1] ^ 2))
  # for(k in 2:basis_num){
  #   a[,k] <- a[,k] - a[,1:(k-1)] %*% (t(a[,1:(k-1)]) %*% a[,k]) 
  #   a[,k] <- a[,k] / sqrt(sum(a[,k] ^ 2))
  # }
  
  a <- t(sapply(1:n, function(i){
    sapply(1:basis_num, function(k){
      rnorm(1, a[i,k],  sqrt(1/n))
    })
  }))
  
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
    obs_point <- fda_full[time_mark,i] + rnorm(time_num, 0, mean(obs_sig))
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

## Intrinsic basis vectors' case
data_gen_fac <- function(basis_num, n, obs_poi, rat){
  
  sigma <- 2 * exp((basis_num):1 / 2) 
  time_grid <- seq(0, 1, length.out = 101)
  basis <- fourier(time_grid, nbasis = basis_num + 3)
  
  a <- t(sapply(1:n, function(i){
    sapply(1:basis_num, function(k){
      sin(k * pi * (i + (n/ 4)) / (2 * n))
    })
  }))
  
  a[,1] <- a[,1] / sqrt(sum(a[,1] ^ 2))
  for(k in 2:basis_num){
    a[,k] <- a[,k] - a[,1:(k-1)] %*% (t(a[,1:(k-1)]) %*% a[,k]) 
    a[,k] <- a[,k] / sqrt(sum(a[,k] ^ 2))
  }
  
  xi <- matrix(runif(ncol(basis) * basis_num, -1, 1), ncol(basis))
  xi[,1] <- xi[,1] / sqrt(sum(xi[,1] ^ 2))
  for(k in 2:basis_num){
    xi[,k] <- xi[,k] - xi[,1:(k-1)] %*% (t(xi[,1:(k-1)]) %*% xi[,k])
    xi[,k] <- xi[,k] / sqrt(sum(xi[,k] ^ 2))
  }
  
  basis <- basis %*% xi
  
  fda_full <- t(t(basis) * sigma)
  fda_full <- t(a %*% t(fda_full))
  
  obs_sig <- sqrt((colSums(fda_full ^ 2)) * 0.01 * rat)
  
  obs_dat <- lapply(1:n, function(i){
    time_num <- sample((obs_poi - 2):(obs_poi + 2), 1)
    time <- sort(sample(time_grid, time_num))
    time_mark <- sapply(1:time_num, function(k) which.min(abs(time_grid - time[k])))
    obs_point <- fda_full[time_mark,i] + rnorm(time_num, 0, obs_sig[i])
    return(list(time = time, obs = obs_point))
  })
  
  tur_dat <- sapply(1:n, function(i){
    obs_dat[[i]]$tur
  })
  
  time_mark <- lapply(1:n, function(i){
    obs_dat[[i]]$time_mark
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
              tur_dat = tur_dat,
              time_mark = time_mark,
              norm = norm
  ))
}

## Functional clustering case
data_gen_clu <- function(basis_num, n, obs_poi, rat, norm){
  
  sigma <- 2 * exp((basis_num):1 / 2)
  time_grid <- seq(0, 1, length.out = 51)
  basis <- fourier(time_grid, nbasis = basis_num + 3)[,2:(basis_num + 1)]
  
  cluster_num <- basis_num
  cluter_label <- sample(1:cluster_num, n, replace = T)
  
  a <- matrix(runif(cluster_num * basis_num, -1, 1), cluster_num)
  a <- t(sapply(1:n, function(i){
    a[cluter_label[i],]
  }))
  
  a[,1] <- a[,1] / sqrt(sum(a[,1] ^ 2))
  for(k in 2:basis_num){
    a[,k] <- a[,k] - a[,1:(k-1)] %*% (t(a[,1:(k-1)]) %*% a[,k]) 
    a[,k] <- a[,k] / sqrt(sum(a[,k] ^ 2))
  }
  
  a <- t(t(a) * sigma) 
  
  xi <- sapply(1:basis_num, function(k){
    Sigma <- diag(rep(sigma[k] ^ 2 / 5 / n, n))
    xi <- rmvnorm(1, mean = rep(0, n), sigma = Sigma)
    return(xi)
  })
  
  xi <- xi + a
  
  fda_full <- basis %*% t(xi)  
  
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
    obs_point <- fda_full[time_mark,i] + rnorm(time_num, 0, sqrt(mean(obs_sig ^ 2)))
    return(list(time = time, obs = obs_point, time_mark = time_mark))
  })
  
  return(list(obs_dat = obs_dat,
              fda_full = fda_full,
              basis = basis,
              xi = xi,
              basis_num = basis_num,
              n = n,
              cluter_label = cluter_label,
              obs_point = obs_poi,
              rat = rat,
              time_grid = time_grid,
              obs_sig = obs_sig,
              norm = norm
  ))
}

################################################################################
# FSVD
## Auxiliary code
tran_dat <- function(Lt, time_grid){
  
  n <- length(Lt)
  time_sam_tol <- unlist(Lt)
  
  sam_mark <- unlist(lapply(1:n, function(i){
    rep(i, length(Lt[[i]]))
  }))
  
  time_mark <- lapply(1:n, function(i){
    sapply(1:length(Lt[[i]]), function(k) which.min(abs(time_grid - Lt[[i]][k])))
  })
  
  sam_num <- unlist(lapply(1:n, function(i){
    rep(length(Lt[[i]]), length(Lt[[i]]))
  }))
  
  return(list(sam_mark = sam_mark,
              time_sam_tol = time_sam_tol,
              time_mark = time_mark,
              sam_num = sam_num,
              time_grid = time_grid
  ))
}

initize_FSVD <- function(Ly, Lt, time_grid, dat_t, Y){
  
  a <- svd(Y)[[3]][,1]
  
  a_coef <- sapply(1:length(dat_t$sam_mark), function(k) a[dat_t$sam_mark[k]])
  Y <- unlist(Ly)
  fit <- smooth.spline(x = dat_t$time_sam_tol, y = Y / a_coef, w = a_coef ^ 2 / dat_t$sam_num, 
                       all.knots = T, cv = F)
  phi <- predict(fit, x = time_grid)$y
  
  return(phi)
}

FSVD_rk <- function(Ly, dat_t, phi, lambda, time_grid, init_num, abs){
  
  n <- length(Ly)
  Y <- unlist(Ly)
  
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
    
    # rho_t <- sqrt(sum(phi_t ^ 2) * 0.01)
    # phi_t <- phi_t
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
  
  return(list(a = a, rho = rho, phi = phi, res = res))
}

FSVD_CV <- function(Ly, dat_t, phi, lambda, tran_datset, time_grid){
  
  fit_p <- FSVD_rk(Ly, dat_t, phi, lambda, time_grid, init_num = 500, abs = 10 ^ (-5))
  
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

FSVD_tune <- function(Ly, time_grid, dat_t, phi, tran_datset){
  
  tran_datset <- lapply(1:length(tran_datset), function(k){
    Ly_i <- lapply((1:n), function(i) Ly[[i]][tran_datset[[k]]$mark[[i]]])
    return(list(Lt = tran_datset[[k]]$Lt, Ly = Ly_i, dat_t = tran_datset[[k]]$dat_t, mark = tran_datset[[k]]$mark))
  })
  
  m <- mean(dat_t$sam_num)
  n <- length(Ly)
  time_tun <- exp(seq(log(10^(-6)), log(10 ^ (1)), length.out = 20)) 
  
  loss <- lapply(time_tun, function(lambda){
    FSVD_CV(Ly, dat_t, phi, lambda, tran_datset, time_grid)
  })
  
  mark <- which.min(sapply(1:length(time_tun), function(i) loss[[i]]))
  lambda <- time_tun[mark]
  
  return(lambda)
}

tun_sel <- function(t, y, basis, lambda){
  smooth.basisPar(t, y, basis, Lfdobj = int2Lfd(2), lambda = lambda)$gcv
}

## FSVD implementation
FSVD <- function(Ly, Lt, R_max, R_pre, num_sel,
                 time_grid = seq(0, 1, length.out = 101)
){
  
  dat_t <- tran_dat(Lt, time_grid)
  n <- length(Ly)
  m <- mean(sapply(1:n, function(i) length(Ly[[i]])))
  
  tran_datset <- lapply(1:5, function(k){
    mark <- lapply(1:n, function(i){
      if(k <= length(Lt[[i]])){
        setdiff(1:length(Lt[[i]]), seq(k, length(Lt[[i]]), 5))
      }else{
        1:length(Lt[[i]])
      }
    })
    Lt_i <- lapply((1:n), function(i) Lt[[i]][mark[[i]]])
    dat_t_i <- tran_dat(Lt_i, time_grid)
    return(list(Lt = Lt_i, dat_t = dat_t_i, mark = mark))
  })
  
  time_grid_mat <- unique(unlist(Lt))
  time_grid_mat_mark <- sapply(1:n, function(i){
    sapply(1:length(Lt[[i]]), function(k) which(Lt[[i]][k] == time_grid_mat))
  })
  time_grid_mat_tol_mark <- sapply(1:length(time_grid_mat), function(k) which.min(abs(time_grid_mat[k] - time_grid)))
  dat_raw <- sapply(1:n, function(i){
    A <- rep(NA, length(time_grid_mat))
    A[time_grid_mat_mark[[i]]] <- Ly[[i]]
    return(A)
  })
  
  Comp_Y <-  fill.nuclear(dat_raw)$X
  
  ## R = 1
  phi <- initize_FSVD(Ly, Lt, time_grid, dat_t, Comp_Y)
  lambda <- FSVD_tune(Ly, time_grid, dat_t, phi = phi, tran_datset)
  fit_FSVD <- FSVD_rk(Ly, dat_t, phi = phi, lambda = lambda, time_grid, init_num = 500, abs = 10 ^ (-5))
  
  R <- 1
  
  Rho <- fit_FSVD$rho
  Phi <- matrix(fit_FSVD$phi, length(time_grid))
  A <- matrix(fit_FSVD$a, n)
  Comp_Y <- Comp_Y - Rho * Phi[time_grid_mat_tol_mark] %*% t(A)
  
  ## R = 2
  phi <- initize_FSVD(fit_FSVD$res, Lt, time_grid, dat_t, Comp_Y)
  fit <- lm(phi ~ Phi + 0)
  phi_init <- fit$residuals
  phi_init <- phi_init / sqrt(sum(phi_init ^ 2) * 0.01)
  
  lambda <- FSVD_tune(fit_FSVD$res, time_grid, dat_t, phi = phi_init, tran_datset)
  fit_FSVD <- FSVD_rk(fit_FSVD$res, dat_t, phi = phi_init, lambda, time_grid, init_num = 500, abs = 10 ^ (-5))
  
  while(R < R_max){
    R <- R + 1
    
    Rho <- c(Rho, fit_FSVD$rho)
    Phi <- cbind(Phi, fit_FSVD$phi)
    A <- cbind(A, fit_FSVD$a)
    Comp_Y <- Comp_Y - Rho[R] * Phi[time_grid_mat_tol_mark,R] %*% t(A[,R])
    
    phi <- initize_FSVD(fit_FSVD$res, Lt, time_grid, dat_t, Comp_Y)
    fit <- lm(phi ~ Phi + 0)
    phi_init <- fit$residuals
    phi_init <- phi_init / sqrt(sum(phi_init ^ 2) * 0.01)
    
    lambda <- FSVD_tune(fit_FSVD$res, time_grid, dat_t, phi = phi_init, tran_datset)
    fit_FSVD <- FSVD_rk(fit_FSVD$res, dat_t, phi = phi_init, lambda, time_grid, init_num = 500, abs = 10 ^ (-5))
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

## Functional clustering via FSVD
FClust_fsvd <- function(Ly, Lt, R_max, R_pre, Clu_num, abs){
  
  n <- length(Ly)
  fit_FSVD <- FSVD(Ly, Lt, R_max, R_pre, num_sel = "FD")
  
  Score <- fit_FSVD$Score[,1:fit_FSVD$R]
  if(Clu_num == F){
    fit_clu_FSVD <- Mclust(data = Score)
    Clu_num <- fit_clu_FSVD$G
    fit_clu_FSVD <- fit_clu_FSVD$classification
  }else{
    fit_clu_FSVD <- Mclust(data = Score, G = Clu_num)$classification
  }
  
  mu <- sapply(1:Clu_num, function(k){
    colMeans(matrix(Score[fit_clu_FSVD == k,], nrow = sum(fit_clu_FSVD == k)))
  })
  
  sigma <- sapply(1:Clu_num, function(k){
    crossprod(matrix(Score[fit_clu_FSVD == k,], nrow = sum(fit_clu_FSVD == k))) / sum(fit_clu_FSVD == k)
  }, simplify = "array")
  
  tau <- sapply(1:Clu_num, function(k){
    mean(unlist(lapply(which(fit_clu_FSVD == k), function(i){
      (Ly[[i]] - fit_FSVD$Intric_basis[fit_FSVD$dat_t$time_mark[[i]],1:fit_FSVD$R] %*%
         fit_FSVD$Score[i,1:fit_FSVD$R]) ^ 2
    })))
  })
  
  pi <- sapply(1:Clu_num, function(k) sum(fit_clu_FSVD == k)) / n
  
  tol_T <- unlist(unlist(Lt))
  
  error <- 1
  ite <- 0
  
  A <- lapply(1:n, function(i){
    fit_FSVD$Intric_basis[fit_FSVD$dat_t$time_mark[[i]],1:fit_FSVD$R]
  })
  
  while((error > abs) & (ite < 100)){
    mu_t <- mu
    
    Sig <- lapply(1:n, function(i){
      lapply(1:Clu_num, function(k){
        A[[i]] %*% sigma[,,k] %*% t(A[[i]]) + diag(tau[k], length(Ly[[i]]))
      })
    })
    
    W <- sapply(1:n, function(i){
      sapply(1:Clu_num, function(k){
        return(mvtnorm::dmvnorm(x = Ly[[i]], mean = c(A[[i]] %*% mu[,k]), sigma = Sig[[i]][[k]], log = T) + log(pi[k]))
      })
    })
    
    W <- sapply(1:n, function(i){
      Pro <- exp(W[,i] - max(W[,i]))
      Pro / sum(Pro)
    })
    
    pi <- rowSums(W) / sum(rowSums(W))
    
    for(k in 1:Clu_num){
      Sig_sq <- lapply(1:n, function(i){
        chol(Sig[[i]][[k]])
      })
      
      Y <- unlist(lapply(1:n, function(i){
        forwardsolve(t(Sig_sq[[i]]), Ly[[i]]) * sqrt(W[k,i])
      }))
      X <- NULL
      for(i in 1:n){
        X <- rbind(X, forwardsolve(t(Sig_sq[[i]]), A[[i]]) * sqrt(W[k,i]))
      }
      
      fit <- lm(Y ~ X + 0)
      mu[,k] <- fit$coefficients
      
      res <- lapply(1:n, function(i){
        Ly[[i]] - A[[i]] %*% mu[,k]
      })
      
      fit <- optim(tau[k], fn = function(res, tau, W, sigma, A, Ly){
        Sig_k <- lapply(1:n, function(i){
          chol(A[[i]] %*% sigma[,,k] %*% t(A[[i]]) + diag(tau, length(Ly[[i]])))
        })
        
        sum(sapply(1:n, function(i) crossprod(forwardsolve(t(Sig_k[[i]]), res[[i]])) * (W[k,i]))) +
          sum(sapply(1:n, function(i) 2 * sum(log(diag(Sig_k[[i]]))) * W[k,i]))
        
      }, res = res, W = W, sigma = sigma, A = A, Ly = Ly, method = "L-BFGS-B", lower = rep(10 ^ (-10), 1))
      
      tau[k] <- fit$par
      
      para <- chol(sigma[,,k] + diag(10 ^ (-8), nrow(sigma[,,k])))
      para <- para[upper.tri(para, diag = T)]
      fit <- optim(para, fn = function(res, para, tau, W, A, Ly){
        sigma <- matrix(0, nrow = nrow(sigma[,,k]), ncol = ncol(sigma[,,k]))
        sigma[upper.tri(sigma, diag = T)] <- para
        
        Sig_k <- lapply(1:n, function(i){
          chol(crossprod(sigma %*% t(A[[i]])) + diag(tau, length(Ly[[i]])))
        })
        
        sum(sapply(1:n, function(i) crossprod(forwardsolve(t(Sig_k[[i]]), res[[i]])) * (W[k,i]))) +
          sum(sapply(1:n, function(i) 2 * sum(log(diag(Sig_k[[i]]))) * W[k,i]))
        
      }, res = res, tau = tau[k], W = W, A = A, Ly = Ly)
      
      AA <- matrix(0, nrow = nrow(sigma[,,k]), ncol = ncol(sigma[,,k]))
      AA[upper.tri(AA, diag = T)] <- fit$par
      sigma[,,k] <- crossprod(AA)
    }
    
    error <- sum(abs(mu - mu_t)) / sum(abs(mu_t))
    ite <- ite + 1
    cat(error, pi, "\n")
  }
  
  fit_clu_FSVD_EM <- sapply(1:n, function(i) which.max(W[,i]))
  
  Score_EM <- t(sapply(1:n, function(i){
    rowSums(sapply(1:Clu_num, function(k){
      (mu[,k] + sigma[,,k] %*% t(A[[i]]) %*% 
         solve(A[[i]] %*% sigma[,,k] %*% t(A[[i]]) + diag(tau[k], length(Ly[[i]]))) %*% 
         (Ly[[i]] - A[[i]] %*% mu[,k])) * W[k,i]
    }))
  }))
  
  return(list(
    fit_FSVD = fit_FSVD,
    fit_clu_FSVD = fit_clu_FSVD,
    fit_clu_FSVD_EM = fit_clu_FSVD_EM,
    sigma = sigma,
    mu = mu,
    Probability_clu = W,
    tau = tau,
    Score_EM = Score_EM
  ))
}

################################################################################
# Simulation function and error calculation
## Simulation for intrinsic basis functions' case
sim_func <- function(basis_num, n, obs_point, rat, norm, seed){
  
  set.seed(seed)
  dat_col <- data_gen_fun(basis_num, n, obs_point, rat, norm)
  dat <- dat_col$obs_dat
  
  Ly <- lapply(1:n, function(i) dat[[i]]$obs)
  Lt <- lapply(1:n, function(i) dat[[i]]$time)
  time_grid <- seq(0, 1, 0.01)
  
  # Plot
  t <- c(sapply(1:n, function(i){dat_col$time_grid}))
  Value <- c(sapply(1:n, function(i){dat_col$fda_full[,i]}))
  m <- c(sapply(1:n, function(i){rep(i, 30)}))
  
  ## Smoothing spline
  fit_smo <- lapply(1:n, function(i){
    fit <- smooth.spline(x = Lt[[i]], y = Ly[[i]], cv = F)
    return(predict(fit, time_grid)$y)
  })
  
  ## FPCA-zero mean
  fit_FPCA0 <- FPCA(Ly, Lt, optns = list(error = T, nRegGrid = 101, 
                                         methodMuCovEst = "smooth",
                                         methodBwCov = "GCV",
                                         methodSelectK = "AIC",
                                         maxK = 2 * basis_num,
                                         methodXi = "CE",
                                         userMu = list(t = seq(0, 1, 0.01), mu = rep(0 , 101))))
  
  ## FPCA
  fit_FPCA <- FPCA(Ly, Lt, optns = list(error = T, nRegGrid = 101, 
                                        methodBwCov = "GCV",
                                        methodMuCovEst = "smooth",
                                        methodBwCov = "GCV",
                                        methodXi = "CE",
                                        methodSelectK = "AIC",
                                        maxK = 2 * basis_num
  ))
  
  ## FSVD
  fit_FSVD <- FSVD(Ly, Lt, R_max = 2 * basis_num, R_pre = F, num_sel = "FD")
  
  return(list(dat_col = dat_col,
              fit_smo = fit_smo,
              fit_FPCA0 = fit_FPCA0,
              fit_FPCA = fit_FPCA,
              fit_FSVD = fit_FSVD))
}

sim_func_iid <- function(basis_num, n, obs_point, rat, norm, seed){
  
  set.seed(seed)
  dat_col <- data_gen_fun_iid(basis_num, n, obs_point, rat, norm)
  dat <- dat_col$obs_dat
  
  Ly <- lapply(1:n, function(i) dat[[i]]$obs)
  Lt <- lapply(1:n, function(i) dat[[i]]$time)
  time_grid <- seq(0, 1, 0.01)
  
  # Plot
  t <- c(sapply(1:n, function(i){dat_col$time_grid}))
  Value <- c(sapply(1:n, function(i){dat_col$fda_full[,i]}))
  m <- c(sapply(1:n, function(i){rep(i, 30)}))
  
  ## Smoothing spline
  fit_smo <- lapply(1:n, function(i){
    fit <- smooth.spline(x = Lt[[i]], y = Ly[[i]], cv = F)
    return(predict(fit, time_grid)$y)
  })
  
  ## FPCA-zero mean
  fit_FPCA0 <- FPCA(Ly, Lt, optns = list(error = T, nRegGrid = 101, 
                                         methodMuCovEst = "smooth",
                                         methodBwCov = "GCV",
                                         methodSelectK = "AIC",
                                         maxK = 2 * basis_num,
                                         methodXi = "CE",
                                         dataType = 'Sparse'
                                         # userMu = list(t = seq(0, 1, 0.01), mu = rep(0 , 101))
  ))
  
  ## FPCA
  fit_FPCA <- fit_FPCA0
  
  ## FSVD
  fit_FSVD <- FSVD(Ly, Lt, R_max = 2 * basis_num, R_pre = F, num_sel = "FD")
  
  return(list(dat_col = dat_col,
              fit_smo = fit_smo,
              fit_FPCA0 = fit_FPCA0,
              fit_FPCA = fit_FPCA,
              fit_FSVD = fit_FSVD))
}

### Error of functional completion
Completion_error <- function(Result){
  
  sim_num <- length(Result)
  Ture <- lapply(1:sim_num, function(i){
    Result[[i]]$dat_col$fda_full
  })
  
  SMO_com <- lapply(1:sim_num, function(i){
    sapply(1:length(Result[[i]]$fit_smo), function(k){
      Result[[i]]$fit_smo[[k]]
    })
  })
  
  FPCA0_com <- lapply(1:sim_num, function(i){
    Result[[i]]$fit_FPCA0$phi %*% t(Result[[i]]$fit_FPCA0$xiEst)
  })
  
  FPCA_com <- lapply(1:sim_num, function(i){
    Result[[i]]$fit_FPCA$phi %*% t(Result[[i]]$fit_FPCA$xiEst) + Result[[i]]$fit_FPCA$mu
  })
  
  FSVD_com <- lapply(1:sim_num, function(i){
    R <- Result[[i]]$fit_FSVD$R
    matrix(Result[[i]]$fit_FSVD$Intric_basis[,1:R], nrow = 101) %*% t(matrix((Result[[i]]$fit_FSVD$Score[,1:R]), nrow = nrow(Result[[i]]$fit_FSVD$Score)))
  })
  
  error_smo <- mean(sapply(1:sim_num, function(i){
    (sum((Ture[[i]] - SMO_com[[i]]) ^ 2) / sum(Ture[[i]] ^ 2) * 100)
  }))
  
  error_FPCA0 <- mean(sapply(1:sim_num, function(i){
    (sum((Ture[[i]] - FPCA0_com[[i]]) ^ 2) / sum(Ture[[i]] ^ 2) * 100)
  }))
  
  error_FPCA <- mean(sapply(1:sim_num, function(i){
    (sum((Ture[[i]] - FPCA_com[[i]]) ^ 2) / sum(Ture[[i]] ^ 2) * 100)
  }))
  
  error_FSVD <- mean(sapply(1:sim_num, function(i){
    (sum((Ture[[i]] - FSVD_com[[i]]) ^ 2) / sum(Ture[[i]] ^ 2) * 100)
  }))
  
  return(list(error_smo = error_smo, 
              error_FPCA0 = error_FPCA0,
              error_FPCA = error_FPCA,
              error_FSVD = error_FSVD))
}

## Error of estimating intrinsic basis function
basis_error <- function(Result, k){
  
  sim_num <- length(Result)
  Ture <- lapply(1:sim_num, function(i){
    Result[[i]]$dat_col$basis[,k]
  })
  
  FPCA0_com <- lapply(1:sim_num, function(i){
    Result[[i]]$fit_FPCA0$phi <- matrix(Result[[i]]$fit_FPCA0$phi, nrow = 101)
    if(k <= ncol(Result[[i]]$fit_FPCA0$phi)){
      Result[[i]]$fit_FPCA0$phi[,k]
    }else{
      rep(0, nrow(Result[[i]]$fit_FPCA0$phi))
    }
  })
  
  FPCA_com <- lapply(1:sim_num, function(i){
    Result[[i]]$fit_FPCA$phi <- matrix(Result[[i]]$fit_FPCA$phi, nrow = 101)
    if(k <= ncol(Result[[i]]$fit_FPCA$phi)){
      Result[[i]]$fit_FPCA$phi[,k]
    }else{
      rep(0, nrow(Result[[i]]$fit_FPCA$phi))
    }
  })
  
  FSVD_com <- lapply(1:sim_num, function(i){
    Result[[i]]$fit_FSVD$Intric_basis[,k]
  })
  
  error_FPCA0 <- mean(sapply(1:sim_num, function(i){
    error <- sqrt(1 - sum(Ture[[i]] * FPCA0_com[[i]]) ^ 2 / sum(Ture[[i]] ^ 2) / sum(FPCA0_com[[i]] ^ 2))
    if(is.finite(error) == F){
      error <- 1
    }
    return(error)
  }))
  
  error_FPCA <- mean(sapply(1:sim_num, function(i){
    error <- sqrt(1 - sum(Ture[[i]] * FPCA_com[[i]]) ^ 2 / sum(Ture[[i]] ^ 2) / sum(FPCA_com[[i]] ^ 2))
    if(is.finite(error) == F){
      error <- 1
    }
    return(error)
  }))
  
  error_FSVD <- mean(sapply(1:sim_num, function(i){
    sqrt(1 - sum(Ture[[i]] * FSVD_com[[i]]) ^ 2 / sum(Ture[[i]] ^ 2) / sum(FSVD_com[[i]] ^ 2))
  }))
  
  return(list(error_FPCA0 = error_FPCA0,
              error_FPCA = error_FPCA,
              error_FSVD = error_FSVD))
}

## Simulation for intrinsic basis vectors' case
sim_func_fac <- function(basis_num, n, obs_point, rat, seed){
  
  set.seed(seed)
  dat_col <- data_gen_fac(basis_num, n, obs_point, rat)
  dat <- dat_col$obs_dat
  
  Ly <- lapply(1:n, function(i) dat[[i]]$obs)
  Lt <- lapply(1:n, function(i) dat[[i]]$time)
  time_grid <- seq(0, 1, 0.01)
  
  # Plot
  t <- c(sapply(1:n, function(i){dat_col$time_grid}))
  Value <- c(sapply(1:n, function(i){dat_col$fda_full[,i]}))
  m <- c(sapply(1:n, function(i){rep(i, 101)}))
  
  ## Imputing missing data
  time_grid_mat <- seq(0, 1, length.out = obs_point)
  Y <- sapply(1:n, function(i){
    sapply(1:obs_point, function(t){
      mark <- sapply(1:length(Lt[[i]]), function(k) (abs(Lt[[i]][k] - time_grid_mat[t]) < 0.2))
      if(sum(mark == T) != 0){
        mean(Ly[[i]][mark])
      }else{
        mark <- which.min(sapply(1:length(Lt[[i]]), function(k) (abs(Lt[[i]][k] - time_grid_mat[t]))))
        Ly[[i]][mark]
      }
    })
  })
  
  ## SVD
  fit_svd <- svd(Y)
  IC <- sapply(1:(basis_num * 2), function(k){
    log(mean((Y - matrix(fit_svd$u[,1:k], nrow = nrow(Y)) %*% diag(c(fit_svd$d[1:k]), nrow = k) %*% t(fit_svd$v[,1:k])) ^ 2)) + k * log(min(n, obs_point)) / min(n, obs_point)
  })
  R <- which.min(IC)
  fit_svd <- list(loading = fit_svd$v, R = R)
  
  ## Factor model of time series
  fit_fac <- HDTSA::Factors(Y = Y, lag.k = 2)
  
  ## FSVD
  fit_FSVD <- FSVD(Ly, Lt, R_max = basis_num * 2, R_pre = F, num_sel = "FM")
  
  return(list(dat_col = dat_col,
              fit_svd = fit_svd,
              fit_fac = fit_fac,
              fit_FSVD = fit_FSVD,
              Y = Y))
}

### Error of estimating intrinsic basis vector
load_error <- function(Result, k){
  
  sim_num <- length(Result)
  Ture <- lapply(1:sim_num, function(i){
    Result[[i]]$dat_col$a[,1:3]
  })
  
  SVD_com <- lapply(1:sim_num, function(i){
    A <- Result[[i]]$fit_svd$loading
    if(Result[[i]]$fit_svd$R >= 3){
      qr.Q(qr(A[,1:3]))
    }else{
      cbind(qr.Q(qr(A[,1:Result[[i]]$fit_svd$R])), matrix(0, nrow(A), 3 - Result[[i]]$fit_svd$R))
    }
  })
  
  FAM_com <- lapply(1:sim_num, function(i){
    A <- matrix(Result[[i]]$fit_fac$loading.mat, n)
    if(ncol(A) >= 3){
      qr.Q(qr(A))
    }else{
      cbind(qr.Q(qr(A)), matrix(0, nrow(A), 3 - ncol(A)))
    }
  })
  
  FSVD_com <- lapply(1:sim_num, function(i){
    A <- Result[[i]]$fit_FSVD$Loading
    if(Result[[i]]$fit_FSVD$R >= 3){
      qr.Q(qr(A[,1:3]))
    }else{
      cbind(qr.Q(qr(A[,1:Result[[i]]$fit_FSVD$R])), matrix(0, nrow(A), 3 - Result[[i]]$fit_FSVD$R))
    }
  })
  
  error_SVD <- mean(sapply(1:sim_num, function(i){
    est <- t(Ture[[i]]) %*% SVD_com[[i]]
    est <- svd(est)
    sum(((Ture[[i]] - SVD_com[[i]] %*% est$v %*% t(est$u))) ^ 2) * 100 / 3
  }))
  
  error_FAM <- mean(sapply(1:sim_num, function(i){
    est <- t(Ture[[i]]) %*% FAM_com[[i]]
    est <- svd(est)
    sum(((Ture[[i]] - FAM_com[[i]] %*% est$v %*% t(est$u))) ^ 2) * 100 / 3
  }))
  
  error_FSVD <- mean(sapply(1:sim_num, function(i){
    est <- t(Ture[[i]]) %*% FSVD_com[[i]]
    est <- svd(est)
    sum(((Ture[[i]] - FSVD_com[[i]] %*% est$v %*% t(est$u))) ^ 2) * 100 / 3
  }))
  
  return(list(error_FAM = error_FAM,
              error_SVD = error_SVD,
              error_FSVD = error_FSVD))
}

## Simulation for functional clustering case
sim_func_clu <- function(basis_num, n, obs_point, rat, seed, norm){
  
  set.seed(seed)
  
  fit_clu_FPCA <- NULL
  fit_clu_FSVD <- NULL
  fit_clu_spline <- NULL
  
  while((is.null(fit_clu_FPCA) | is.null(fit_clu_FSVD) | is.null(fit_clu_spline)) == T){
    dat_col <- data_gen_clu(basis_num, n, obs_point, rat, norm)
    dat <- dat_col$obs_dat
    
    Ly <- lapply(1:n, function(i) dat[[i]]$obs)
    Lt <- lapply(1:n, function(i) dat[[i]]$time)
    time_grid <- seq(0, 1, 0.01)
    
    tryCatch({
      ## Smoothing-spline
      ## Average point
      m <- floor(mean(sapply(1:n, function(i) length(Lt[[i]]))))
      
      data <- list(x = unlist(Ly), time = unlist(Lt), 
                   grid = sort(unique(unlist(Lt))),
                   curve = unlist(lapply(1:n, function(i) rep(i, length(Lt[[i]])))),
                   B = NULL)
      
      f_basis <- create.bspline.basis(c(0, 1), nbasis = m, norder = 4)
      fit_coef_spline <- sapply(1:n, function(i){
        lambda <- optimise(tun_sel, interval = c(0, 100),
                           t = Lt[[i]], y = Ly[[i]], basis = f_basis)$minimum
        smooth.basisPar(Lt[[i]], Ly[[i]], f_basis, Lfdobj = int2Lfd(2), lambda = lambda)$fd$coefs
      })
      
      fit_clu_spline_mol <- smoothing.cluts(data = data,
                                            K = 3, q = m, h = 2,
                                            use.covariates = FALSE,
                                            stand.cov = F, index.cov = NULL,
                                            random = T, svd = T, EMplot = F,
                                            EMstep.tol = 1e-6, Mstep.tol = 1e-4, 
                                            B = t(fit_coef_spline), lambda = 1e-10, n.cores = 2)
      
      Pro <- sapply(1:n, function(i){
        sapply(1:3, function(k){
          B <- fit_clu_spline_mol$design$S[which(data$curve == i),]
          mu <- as.vector(B %*% (fit_clu_spline_mol$pars$lambda.zero + fit_clu_spline_mol$pars$Lambda %*% fit_clu_spline_mol$pars$alpha[k,]))
          sigma <- B %*% fit_clu_spline_mol$pars$Gamma[k,,] %*% t(B) + diag(fit_clu_spline_mol$pars$sig2, length(Ly[[i]]))
          dmvnorm(matrix(Ly[[i]], nrow = 1), mean = mu, sigma = sigma, log = T)
        })
      })
      
      fit_clu_spline <- sapply(1:n, function(i) which.max(Pro[,i]))
      
    }, warning = function(m) {NULL})
    
    tryCatch({
      
      ## FPCA
      fit_clu_FPCA <- FClust(Ly, Lt, 3)$cluster
      
    }, warning = function(m) {NULL})
    
    tryCatch({
      ## FSVD
      fit_clu_FSVD <- FClust_fsvd(Ly, Lt, R_max = basis_num * 2, R_pre = F, Clu_num = 3, abs = 0.01)
      
    }, warning = function(m) {NULL})
  }
  
  return(list(dat_col = dat_col,
              fit_clu_spline =  fit_clu_spline,
              fit_clu_FPCA = fit_clu_FPCA,
              fit_clu_FSVD = fit_clu_FSVD$fit_clu_FSVD,
              fit_clu_FSVD_EM = fit_clu_FSVD$fit_clu_FSVD_EM
  ))
}

### Error of functional clustering
clu_error <- function(Result){
  
  sim_num <- length(Result)
  Ture <- lapply(1:sim_num, function(i){
    Result[[i]]$dat_col$cluter_label
  })
  
  spline_com <- lapply(1:sim_num, function(i){
    Result[[i]]$fit_clu_spline
  })
  
  FPCA_com <- lapply(1:sim_num, function(i){
    Result[[i]]$fit_clu_FPCA
  })
  
  FSVD_com <- lapply(1:sim_num, function(i){
    Result[[i]]$fit_clu_FSVD
  })
  
  FSVD_com_EM <- lapply(1:sim_num, function(i){
    Result[[i]]$fit_clu_FSVD_EM
  })
  
  error_spline <- (sapply(1:sim_num, function(i){
    ARI(Ture[[i]], spline_com[[i]])
  }))
  
  error_FPCA <- (sapply(1:sim_num, function(i){
    ARI(Ture[[i]], FPCA_com[[i]])
  }))
  
  error_FSVD <- (sapply(1:sim_num, function(i){
    ARI(Ture[[i]], FSVD_com[[i]])
  }))
  
  error_FSVD_EM <- (sapply(1:sim_num, function(i){
    ARI(Ture[[i]], FSVD_com_EM[[i]])
  }))
  
  return(list(error_spline = error_spline,
              error_FPCA = error_FPCA,
              error_FSVD = error_FSVD,
              error_FSVD_EM = error_FSVD_EM
  ))
}

################################################################################
# Code for other method
## K-NN imputation code
MedImpute <- function(dat, K, time_grid, h, alpha){
  I <- which(rowSums(is.na(dat)) > 0)
  
  dis_t <- rdist(time_grid)
  W <- matrix(0, nrow = nrow(dat), ncol = ncol(dat))
  Z <- matrix(0, nrow = nrow(dat), ncol = nrow(dat))
  for(i in 1:n){
    W[i,which(is.na(dat[i,]) == F)] <- unlist(dat[i,which(is.na(dat[i,]) == F)])
    Z[i,sample((1:n)[-i], 1)] <- 1
  }
  
  n <- nrow(dat)
  p <- ncol(dat)
  
  C <- sapply(h, function(k){
    A <- 2 ^ (- dis_t / k)
    A <- A / rowSums(A)
  }, simplify = 'array') 
  
  error <- 1
  ite <- 0
  
  while((error > 0.001) & (ite < 100)){
    W_t <- W
    
    ite <- 1 + ite
    for(i in 1:n){
      loss <- sapply((1:n)[-i], function(ii){
        sum(((W[i,] -  W[ii,]) * (1 - alpha)) ^ 2)
      })
      mark <- ((1:n)[-i])[order(loss, decreasing = F)[1:K]]
      
      Z[i,mark] <- 1
      Z[i,-mark] <- 0
    }
    
    for(i in 1:n){
      for(j in 1:p){
        if(is.na(dat[i,j]) == T){
          W[i,j] <- (sum((((1 - alpha[j]) * Z[i,] + alpha[j] * C[i,,j]) * W[,j])) + 
                       sum((1 - alpha[j]) * Z[I,i] + alpha[j] * C[I,i,j])) /
            (K + sum(alpha[j] * C[i,,j]) + sum((1 - alpha[j]) * Z[I,i] + alpha[j] * C[I,i,j]))
        }
      }
    }
    
    error <- sum(abs(W - W_t)) / sum(abs(W))
    # cat(error, "\n")
  }
  
  return(W)
}

## Jame's functional clustering method
source("Jame_FC.R")
