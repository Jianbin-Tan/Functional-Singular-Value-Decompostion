# Set address
rm(list=ls())
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Load functions and packages
source("Function.R")

# Load data
dat <- read.csv(file = "Data/lab_drg_870_872_Nov_iv.csv")

## Preprocessing and implementation functions
subject_id <- unique(dat$SUBJECT_ID)
feature_id <- unique(dat$FEATURE_NAME)
time_grid <- seq(0, 1, length.out = 101)

imp_func <- function(j, subject_id, feature_id, dat, time_grid){
  dat_sub <- dat[(dat$SUBJECT_ID == subject_id[j]),]
  feature_id_sub <- unique(dat_sub$FEATURE_NAME)
  n_sub <- length(feature_id_sub)
  
  Ly <- list()
  Lt <- list()
  for(i in 1:n_sub){
    sel_dat <- dat_sub$VALUE[dat_sub$FEATURE_NAME == feature_id_sub[i]]
    sel_time <- dat_sub$RECORD_MIN[dat_sub$FEATURE_NAME == feature_id_sub[i]]
    
    sel_time <- (sel_time - min(sel_time)) / (max(sel_time) - min(sel_time))
    mark <- which(is.na(sel_dat) == F)
    Ly[[i]] <- sel_dat[mark]
    Lt[[i]] <- sel_time[mark]
  }
  
  Lt <- lapply(1:n_sub, function(i) sapply(1:length(Lt[[i]]), function(k) time_grid[which.min(abs(Lt[[i]][k] - time_grid))]))
  for(i in 1:n_sub){
    uni_t <- unique(Lt[[i]])
    
    Ly[[i]] <- sapply(1:length(uni_t), function(l) mean(Ly[[i]][which(Lt[[i]] == uni_t[[l]])]))
    Lt[[i]] <- uni_t
  }
  
  mark <- which(sapply(1:n_sub, function(i) (length(Ly[[i]]) >= 5) & length(which(Ly[[i]] != 0)) >= 5))
  
  if(length(mark) > 1){
    Lt <- lapply(mark, function(i) Lt[[i]][which(Ly[[i]] != 0)])
    Ly <- lapply(mark, function(i) Ly[[i]][which(Ly[[i]] != 0)] / sqrt(mean(Ly[[i]][which(Ly[[i]] != 0)] ^ 2)))
    n_sub <- length(Ly)
    
    Lt_tol <- sort(unique(unlist(Lt)))
    
    if(length(Lt_tol) >= 5){
      
      Ly_fit <- list()
      Lt_fit <- list()
      Ly_pred <- list()
      Lt_pred <- list()
      
      for(i in 1:n_sub){
        mark <- which(1:length(Ly[[i]]) %% 2 == (i %% 2))
        
        if(length(mark) < 4){
          mark <- sample(1:length(Ly[[i]]), 4)
        }
        
        Ly_fit[[i]] <- Ly[[i]][mark]
        
        Lt_fit[[i]] <- Lt[[i]][mark]
        
        Ly_pred[[i]] <- Ly[[i]][-mark]
        
        Lt_pred[[i]] <- Lt[[i]][-mark]
      }
      
      result <- rep(NA, 4)
      ## FSVD
      tryCatch({
        fit_FSVD <- FSVD(Ly = Ly_fit, Lt = Lt_fit, R_max = 10, R_pre = F, num_sel = "FM", time_grid = time_grid)
        fit_FSVD <- sapply(1:n_sub, function(i) (fit_FSVD$Fac_serial %*% c(fit_FSVD$Loading[i,])))
        
        result[1] <- mean(unlist(lapply(1:n_sub, function(i){
          mark <- sapply(1:length(Lt_pred[[i]]), function(k) which.min(abs(Lt_pred[[i]][k] - time_grid)))
          (fit_FSVD[mark,i] - Ly_pred[[i]])
        })) ^ 2)
      }, error = function(e) NULL)
      
      
      ## VAE
      time_grid_mat <- sort(unique(unlist(Lt_fit)))
      dat_raw <- sapply(1:n_sub, function(i){
        mark <- sapply(1:length(Lt_fit[[i]]), function(k) which(Lt_fit[[i]][k] == time_grid_mat))
        A <- rep(NA, length(time_grid_mat))
        A[mark] <- Ly_fit[[i]]
        return(A)
      })
      
      tryCatch({
        VAE_fit <- impute_vae(dat_raw)
        
        result[2] <- mean(unlist(lapply(1:n_sub, function(i){
          mark <- sapply(1:length(Lt_pred[[i]]), function(k) which.min(abs(Lt_pred[[i]][k] - time_grid_mat)))
          (VAE_fit[mark,i] - Ly_pred[[i]])
        })) ^ 2)
      }, error = function(e) NULL)
      
      
      ## Matrix completion
      tryCatch({
        SVD_imp <-  fill.nuclear(dat_raw)$X
        
        result[3] <- mean(unlist(lapply(1:n_sub, function(i){
          mark <- sapply(1:length(Lt_pred[[i]]), function(k) which.min(abs(Lt_pred[[i]][k] - time_grid_mat)))
          (SVD_imp[mark,i] - Ly_pred[[i]])
        })) ^ 2)
        
      }, error = function(e) NULL)
      
      
      ## KNN
      tryCatch({
        dat_time_serie <- as.data.frame(dat_raw)
        KNN_imp <- MedImpute(dat = dat_time_serie, K = 5, time_grid = time_grid_mat, h = rep(1, ncol(dat_time_serie)), alpha = rep(0, ncol(dat_time_serie)))
        
        result[4] <- mean(unlist(lapply(1:n_sub, function(i){
          mark <- sapply(1:length(Lt_pred[[i]]), function(k) which.min(abs(Lt_pred[[i]][k] - time_grid_mat)))
          (KNN_imp[mark,i] - Ly_pred[[i]])
        })) ^ 2)
      }, error = function(e) NULL)
      
    }else{
      result <- rep(NA, 4)
    }
    
  }else{
    result <- rep(NA, 4)
  }
  
  return(result)
}

# Implementation
library(snowfall)

# Simulation
sfInit(parallel = T, cpus = 40)
sfExport("subject_id", "feature_id", "dat", "time_grid", "imp_func")
sfSource("Function.R")
Result <- sfLapply(1:length(subject_id), imp_func,
                   subject_id = subject_id, feature_id = feature_id, 
                   dat = dat, time_grid = time_grid)

sfStop()

save(Result, file = paste0("Result/Result_EHR.rda"), version = 2)

# Plot
data_plot <- NULL
for(k in 1:length(Result)){
  
  res <- Result[[k]]
  data_plot <- rbind(data_plot, data.frame(value = res, 
                                           method = c("FSVD", "VAE", "Matrix completion", "K-NN")))
}

data_plot$method <- factor(data_plot$method, levels = c("Matrix completion", "VAE", "K-NN", "FSVD"))

ggplot(data_plot) +
  geom_jitter(
    aes(x = method, y = value, color = method),
    width = 0.29, height = 0, size = 1, alpha = 0.6
  ) +
  geom_boxplot(
    aes(x = method, y = value, fill = method),
    width = 0.6, size = 0.4, outlier.shape = NA, alpha = 0.85
  ) +
  geom_hline(yintercept = median(sapply(1:length(Result), function(i) Result[[i]][1]), na.rm = T), linetype = 2) +
  labs(x = "", y = "", title = "") +
  scale_y_log10(limits = c(10^(-4), 10)) +
  scale_fill_manual(values = c("#F0A780", "#2c7fb8", "#f768a1", "#e34a33")) +
  scale_color_manual(values = c("#F0A780", "#2c7fb8", "#f768a1", "#e34a33"), guide = "none") +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, size = 14),
    axis.title.y = element_text(size = 14),
    panel.border = element_blank(),
    plot.title = element_text(size = 14, hjust = 0, face = "bold")
  )

ggsave(paste0("Figure/", "boxplot", ".pdf"), width = 8, height = 4, dpi = 300)
