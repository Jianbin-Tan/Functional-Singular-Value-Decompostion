impute_vae <- function(
    X_obs,
    latent_dim   = 5L,
    hidden_dims  = c(64L, 32L),
    epochs       = 100,
    lr           = 1e-3,
    batch_size   = 32L
) {
  library(torch)
  
  # 1) Input validation & prep
  if (!is.matrix(X_obs) || !is.numeric(X_obs)) {
    stop("`X_obs` must be a numeric matrix.")
  }
  X_obs     <- as.matrix(X_obs)
  input_dim <- ncol(X_obs)
  
  # Build mask & initial fill
  mask_mat <- ifelse(is.na(X_obs), 0, 1)
  X_filled <- X_obs
  col_means <- colMeans(X_obs, na.rm = TRUE)
  for (j in seq_len(input_dim)) {
    X_filled[is.na(X_obs[, j]), j] <- col_means[j]
  }
  
  # Torch tensors & dataloader
  X_tensor    <- torch_tensor(X_filled, dtype = torch_float())
  mask_tensor <- torch_tensor(mask_mat, dtype = torch_float())
  ds <- tensor_dataset(X_tensor, mask_tensor)
  dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE)
  
  # 2) VAE module with nn_sequential
  vae_module <- nn_module(
    "VAE",
    initialize = function(input_dim, hidden_dims, latent_dim) {
      # build encoder sequentially
      enc_layers <- list()
      in_dim <- input_dim
      for (h in hidden_dims) {
        enc_layers <- c(enc_layers,
                        nn_linear(in_dim, h),
                        nn_relu())
        in_dim <- h
      }
      self$encoder <- do.call(nn_sequential, enc_layers)
      self$mu_layer     <- nn_linear(in_dim,    latent_dim)
      self$logvar_layer <- nn_linear(in_dim,    latent_dim)
      
      # build decoder sequentially
      dec_layers <- list()
      in_dim <- latent_dim
      for (h in rev(hidden_dims)) {
        dec_layers <- c(dec_layers,
                        nn_linear(in_dim, h),
                        nn_relu())
        in_dim <- h
      }
      # final output layer
      dec_layers <- c(dec_layers,
                      nn_linear(in_dim, input_dim))
      self$decoder <- do.call(nn_sequential, dec_layers)
    },
    
    encode = function(x) {
      h <- self$encoder(x)
      mu     <- self$mu_layer(h)
      logvar <- self$logvar_layer(h)
      list(mu, logvar)
    },
    
    reparameterize = function(mu, logvar) {
      std <- (0.5 * logvar)$exp()
      eps <- torch_randn_like(std)
      mu + eps * std
    },
    
    decode = function(z) {
      self$decoder(z)
    },
    
    forward = function(x) {
      enc     <- self$encode(x)
      mu      <- enc[[1]]
      logvar  <- enc[[2]]
      z       <- self$reparameterize(mu, logvar)
      recon   <- self$decode(z)
      list(recon, mu, logvar)
    }
  )
  
  model     <- vae_module(input_dim, hidden_dims, latent_dim)
  optimizer <- optim_adam(model$parameters, lr = lr)
  
  # 3) Training loop
  for (epoch in seq_len(epochs)) {
    model$train()
    total_loss <- 0
    batch_cnt  <- 0
    
    coro::loop(for (batch in dl) {
      x_batch    <- batch[[1]]
      mask_batch <- batch[[2]]
      
      optimizer$zero_grad()
      out    <- model(x_batch)
      recon  <- out[[1]]; mu <- out[[2]]; logvar <- out[[3]]
      
      # masked MSE
      se       <- (recon - x_batch)$pow(2) * mask_batch
      mse_loss <- se$sum(dim = 2) / (mask_batch$sum(dim = 2) + 1e-8)
      
      # KL divergence
      kl_loss <- -0.5 * (1 + logvar - mu$pow(2) - logvar$exp())$sum(dim = 2)
      
      loss <- (mse_loss + kl_loss)$mean()
      loss$backward()
      optimizer$step()
      
      total_loss <- total_loss + loss$item()
      batch_cnt  <- batch_cnt + 1
    })
    
    cat(sprintf("Epoch %02d/%02d — avg loss: %.4f\n",
                epoch, epochs, total_loss / batch_cnt))
  }
  
  # 4) Reconstruct & impute
  model$eval()
  with_no_grad({
    full_out   <- model(X_tensor)
    recon_full <- full_out[[1]]
  })
  recon_mat <- as.matrix(recon_full)
  
  X_imputed <- X_obs
  X_imputed[is.na(X_obs)] <- recon_mat[is.na(X_obs)]
  X_imputed
}

impute_vaa <- function(
    X_obs,
    latent_dim   = 5L,
    hidden_dims  = c(64L, 32L),
    k            = 5L,        # number of neighbors for initialization
    epochs       = 100,
    lr           = 1e-3,
    batch_size   = 32L
) {
  library(torch)
  
  # 1) Input validation & prep
  if (!is.matrix(X_obs) || !is.numeric(X_obs)) {
    stop("`X_obs` must be a numeric matrix.")
  }
  X_obs     <- as.matrix(X_obs)
  n_samples <- nrow(X_obs)
  input_dim <- ncol(X_obs)
  
  # Build mask
  mask_mat <- ifelse(is.na(X_obs), 0, 1)
  
  # 1a) Initial mean‐fill (for neighbor distance calc)
  col_means <- colMeans(X_obs, na.rm = TRUE)
  X_filled  <- X_obs
  for (j in seq_len(input_dim)) {
    X_filled[is.na(X_obs[, j]), j] <- col_means[j]
  }
  
  # 1b) Neighborhood‐aware initialization
  X_init <- X_filled
  for (i in seq_len(n_samples)) {
    obs_idx <- which(mask_mat[i, ] == 1)
    
    # If a sample is completely missing, fall back to col means
    if (length(obs_idx) == 0) {
      X_init[i, ] <- col_means
      next
    }
    
    # compute weighted‐mask distance to all other samples
    diffsq <- sweep(X_filled, 2, X_filled[i, ], "-")^2
    # only use features actually observed in sample i
    weighted <- diffsq * matrix(mask_mat[i, ], nrow = n_samples, ncol = input_dim, byrow = TRUE)
    dists   <- rowSums(weighted)
    dists[i] <- Inf  # exclude self
    
    # find k nearest neighbors
    nn_idx <- order(dists)[seq_len(min(k, n_samples - 1))]
    
    # for each missing feature, replace with neighbor‐mean
    miss_idx <- which(mask_mat[i, ] == 0)
    if (length(miss_idx) > 0) {
      nn_means <- colMeans(X_filled[nn_idx, miss_idx, drop = FALSE], na.rm = TRUE)
      X_init[i, miss_idx] <- nn_means
    }
  }
  
  # Torch tensors & dataloader
  X_tensor    <- torch_tensor(X_init,  dtype = torch_float())
  mask_tensor <- torch_tensor(mask_mat, dtype = torch_float())
  ds <- tensor_dataset(X_tensor, mask_tensor)
  dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE)
  
  # 2) VAE module (unchanged)
  vae_module <- nn_module(
    "VAE",
    initialize = function(input_dim, hidden_dims, latent_dim) {
      enc_layers <- list(); in_dim <- input_dim
      for (h in hidden_dims) {
        enc_layers <- c(enc_layers, nn_linear(in_dim, h), nn_relu())
        in_dim <- h
      }
      self$encoder      <- do.call(nn_sequential, enc_layers)
      self$mu_layer     <- nn_linear(in_dim, latent_dim)
      self$logvar_layer <- nn_linear(in_dim, latent_dim)
      
      dec_layers <- list(); in_dim <- latent_dim
      for (h in rev(hidden_dims)) {
        dec_layers <- c(dec_layers, nn_linear(in_dim, h), nn_relu())
        in_dim <- h
      }
      dec_layers <- c(dec_layers, nn_linear(in_dim, input_dim))
      self$decoder <- do.call(nn_sequential, dec_layers)
    },
    
    encode = function(x) {
      h      <- self$encoder(x)
      mu     <- self$mu_layer(h)
      logvar <- self$logvar_layer(h)
      list(mu, logvar)
    },
    
    reparameterize = function(mu, logvar) {
      std <- (0.5 * logvar)$exp()
      eps <- torch_randn_like(std)
      mu + eps * std
    },
    
    decode = function(z) {
      self$decoder(z)
    },
    
    forward = function(x) {
      enc     <- self$encode(x)
      mu      <- enc[[1]]; logvar <- enc[[2]]
      z       <- self$reparameterize(mu, logvar)
      recon   <- self$decode(z)
      list(recon, mu, logvar)
    }
  )
  
  model     <- vae_module(input_dim, hidden_dims, latent_dim)
  optimizer <- optim_adam(model$parameters, lr = lr)
  
  # 3) Training loop (unchanged)
  for (epoch in seq_len(epochs)) {
    model$train(); total_loss <- 0; batch_cnt <- 0
    coro::loop(for (batch in dl) {
      x_batch    <- batch[[1]]; mask_batch <- batch[[2]]
      optimizer$zero_grad()
      out    <- model(x_batch)
      recon  <- out[[1]]; mu <- out[[2]]; logvar <- out[[3]]
      
      se       <- (recon - x_batch)$pow(2) * mask_batch
      mse_loss <- se$sum(dim = 2) / (mask_batch$sum(dim = 2) + 1e-8)
      kl_loss  <- -0.5 * (1 + logvar - mu$pow(2) - logvar$exp())$sum(dim = 2)
      loss     <- (mse_loss + kl_loss)$mean()
      
      loss$backward(); optimizer$step()
      total_loss <- total_loss + loss$item(); batch_cnt <- batch_cnt + 1
    })
    cat(sprintf("Epoch %02d/%02d — avg loss: %.4f\n",
                epoch, epochs, total_loss / batch_cnt))
  }
  
  # 4) Reconstruct & final imputation
  model$eval()
  with_no_grad({
    full_out   <- model(X_tensor)
    recon_full <- full_out[[1]]
  })
  recon_mat <- as.matrix(recon_full)
  
  X_imputed <- X_obs
  X_imputed[is.na(X_obs)] <- recon_mat[is.na(X_obs)]
  X_imputed
}

