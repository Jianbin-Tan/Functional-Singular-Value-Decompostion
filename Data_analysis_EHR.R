# Set address
rm(list=ls())
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Load functions and packages
source("Function.R")
load("Data/dat_ehr.rda")
library(keras)
library(tsDyn)
library(gridExtra)
library(filling)
library(mtsdi)

# Data illustration
feature_id <- unique(dat$Feature)
n <- length(feature_id)

Lt_null <- lapply(1:n, function(i){
  dat$Time[which(dat$Feature == feature_id[i])]
})

Ly_null <- lapply(1:n, function(i){
  dat$Value[which(dat$Feature == feature_id[i])]
})

mean_t <- c(min(unlist(Lt_null)), max(unlist(Lt_null)))
Lt <- lapply(1:n, function(i){
  (Lt_null[[i]] - mean_t[1]) / (mean_t[2] - mean_t[1])
})

time_grid <- seq(0, 1, length.out = 101)

time_mark <- lapply(1:n, function(i){
  sapply(1:length(Lt[[i]]), function(k) which.min(abs(time_grid - Lt[[i]][k])))
})

mean_y <- sapply(1:n, function(i) sqrt(mean(Ly_null[[i]] ^ 2)))
Ly <- list()
for(i in 1:n){
  num <- unique(time_mark[[i]])
  Lt[[i]] <- sapply(1:length(num), function(j) time_grid[num[j]])
  Ly[[i]] <- sapply(1:length(num), function(j){
    mean((Ly_null[[i]] / mean_y[i])[which(time_mark[[i]] == num[j])])
  }) 
  time_mark[[i]] <- sapply(1:length(num), function(j) num)
}

dat_plot <- data.frame(Time = (unlist(Lt) * (mean_t[2] - mean_t[1]) + mean_t[1]) / 60,
                       Value = unlist(Ly), Feature = c(unlist(lapply(1:n, function(i) rep(feature_id[i], length(Ly[[i]]))))))

p_ori <- ggplot(dat_plot) + 
  geom_point(aes(x = Time, y = Value), color = "orange") + 
  facet_wrap(.~Feature, scales = "free_y", ncol = 4) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "Time (min)", y = "Value",
       title = "",
       colour = "", fill = "", linetype = "") 
p_ori

ggsave(paste0("Figure/", "real_data", ".pdf"), width = 9.5, height = 5, dpi = 300)

## FSVD
fit_FSVD <- FSVD(Ly, Lt, R_max = 10, R_pre = F, num_sel = "FM", time_grid = time_grid)

dat_plot_fsvd <- data.frame(
  Time = rep((time_grid * (mean_t[2] - mean_t[1]) + mean_t[1]) / 60, n),
  Value = c(sapply(1:n, function(i){
    val <- rep(NA, length(time_grid))
    val[time_mark[[i]]] <- Ly[[i]]
    return(val)
  })), 
  Feature = c(sapply(1:n, function(i) rep(feature_id[i], length(time_grid)))),
  Estimation = c(sapply(1:n, function(i) (fit_FSVD$Fac_serial[,1:fit_FSVD$R] %*% c(fit_FSVD$Loading[i,1:fit_FSVD$R]))))
)

p_1 <- ggplot(dat_plot_fsvd) + 
  geom_point(aes(x = Time, y = Value), color = "orange") + 
  geom_line(aes(x = Time, y = Estimation), size = .6, color = "blue") + 
  facet_wrap(.~Feature, scales = "free_y", ncol = 3) +
  labs(x = "Time (min)", y = "Value",
       title = "",
       colour = "", fill = "", linetype = "") + 
  # scale_color_manual(values = c("blue")) +
  # scale_linetype_manual(values = c(2, 1)) + 
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))  
p_1
ggsave(paste0("Figure/", "fit_dat_FSVD", ".pdf"), width = 7, height = 6, dpi = 300)

## Smoothing spline
fit_smo <- lapply(1:n, function(i){
  fit <- smooth.spline(x = Lt[[i]], y = Ly[[i]], cv = F, all.knots = T)
  y <- predict(fit, x = time_grid)$y
  return(as.numeric(y))
})

dat_plot_smo <- data.frame(
  Time = rep((time_grid * (mean_t[2] - mean_t[1]) + mean_t[1]) / 60, n),
  Value = c(sapply(1:n, function(i){
    val <- rep(NA, length(time_grid))
    val[time_mark[[i]]] <- Ly[[i]]
    return(val)
  })), 
  Feature = c(sapply(1:n, function(i) rep(feature_id[i], length(time_grid)))),
  Estimation = c(sapply(1:n, function(i) fit_smo[[i]]))
)

p_2 <- ggplot(dat_plot_smo) + 
  geom_point(aes(x = Time, y = Value), color = "orange") + 
  geom_line(aes(x = Time, y = Estimation), size = 0.5, color = "blue") + 
  facet_wrap(.~Feature, scales = "free_y", ncol = 3) +
  labs(x = "Time (min)", y = "Value",
       title = "",
       colour = "", fill = "", linetype = "") + 
  # scale_color_manual(values = c("blue")) +
  # scale_linetype_manual(values = c(2, 1)) + 
  # theme_bw(base_family = "Times") +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))  
p_2
ggsave(paste0("Figure/", "fit_dat_smo", ".pdf"), width = 7, height = 6, dpi = 300)

## Matrix completion
time_grid_mat <- unique(unlist(Lt))
dat_raw <- sapply(1:n, function(i){
  mark <- sapply(1:length(Lt[[i]]), function(k) which(Lt[[i]][k] == time_grid_mat))
  A <- rep(NA, length(time_grid_mat))
  A[mark] <- Ly[[i]]
  return(A)
})

SVD_imp <-  fill.nuclear(dat_raw)

dat_plot_mat <- data.frame(
  Time = rep((rep(time_grid * (mean_t[2] - mean_t[1]) + mean_t[1], n) / 60), 1),
  Value = rep(c(sapply(1:n, function(i){
    val <- rep(NA, length(time_grid))
    val[time_mark[[i]]] <- Ly[[i]]
    return(val)
  })), 1), 
  Feature = rep(c(sapply(1:n, function(i) rep(feature_id[i], length(time_grid)))), 1),
  Estimation = c(c(sapply(1:n, function(i){
    val <- rep(NA, length(time_grid))
    time_miss <- time_grid_mat[is.na(dat_raw[,i])]
    mark <- sapply(1:length(time_miss), function(k) which(time_miss[k] == time_grid))
    val[mark] <- SVD_imp$X[is.na(dat_raw[,i]),i]
    return(val)
  })))
)

p_3 <- ggplot(dat_plot_mat) + 
  geom_point(aes(x = Time, y = Value), color = "orange") + 
  geom_point(aes(x = Time, y = Estimation), size = 0.7, color = "blue") + 
  facet_wrap(.~Feature, scales = "free_y", ncol = 3) +
  labs(x = "Time (min)", y = "Value",
       title = "",
       colour = "", fill = "", linetype = "") + 
  # scale_color_manual(values = c("blue")) +
  # scale_linetype_manual(values = c(2, 1)) + 
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))  
p_3
ggsave(paste0("Figure/", "fit_dat_SVD", ".pdf"), width = 7, height = 6, dpi = 300)

## Imputation by KNN
dat_time_serie <- as.data.frame(dat_raw)

KNN_imp <- MedImpute(dat_time_serie, K = 5, time_grid = time_grid_mat, h = rep(1, ncol(dat_time_serie)), alpha = rep(0, ncol(dat_time_serie)))

dat_plot_knn <- data.frame(
  Time = rep((rep(time_grid * (mean_t[2] - mean_t[1]) + mean_t[1], n)) / 60, 1),
  Value = rep(c(sapply(1:n, function(i){
    val <- rep(NA, length(time_grid))
    val[time_mark[[i]]] <- Ly[[i]]
    return(val)
  })), 1), 
  Feature = rep(c(sapply(1:n, function(i) rep(feature_id[i], length(time_grid)))), 1),
  Estimation = c(c(sapply(1:n, function(i){
    val <- rep(NA, length(time_grid))
    time_miss <- time_grid_mat[is.na(dat_time_serie[,i])]
    mark <- sapply(1:length(time_miss), function(k) which(time_miss[k] == time_grid))
    val[mark] <- KNN_imp[is.na(dat_time_serie[,i]),i]
    return(val)
  })))
)

p_4 <- ggplot(dat_plot_knn) + 
  geom_point(aes(x = Time, y = Value), color = "orange") + 
  geom_point(aes(x = Time, y = Estimation), size = 0.7, color = "blue") + 
  facet_wrap(.~Feature, scales = "free_y", ncol = 3) +
  labs(x = "Time (min)", y = "Value",
       title = "",
       colour = "", fill = "", linetype = "") + 
  # scale_color_manual(values = c("blue")) +
  # scale_linetype_manual(values = c(2, 1)) + 
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))  
p_4
ggsave(paste0("Figure/", "fit_dat_KNN", ".pdf"), width = 7, height = 6, dpi = 300)

## Total plot
dat_plot_tol <- rbind(dat_plot_fsvd, dat_plot_mat, dat_plot_smo, dat_plot_knn)

num <- nrow(dat_plot_fsvd)
dat_plot_tol <- data.frame(dat_plot_tol, Method = c(rep("FSVD", num),
                                            rep("Matrix completion", num),
                                            rep("Smoothing spline", num),
                                            rep("K-NN", num)
))
dat_plot_tol <- data.frame(dat_plot_tol, Curve = c(rep(1, num),
                                                  rep(NA, num),
                                                  rep(1, num),
                                                  rep(NA, num)),
                       Point = c(rep(NA, num),
                                 rep(1, num),
                                 rep(NA, num),
                                 rep(1, num)))

dat_plot_tol$Method <- factor(dat_plot_tol$Method, levels = c("Matrix completion", 'Smoothing spline', 'K-NN', "FSVD"))

## Figure A
dat_plot <- dat_plot_tol[(dat_plot_tol$Feature == feature_id[1]) |
         (dat_plot_tol$Feature == feature_id[2]) |
         (dat_plot_tol$Feature == feature_id[3]) |
           (dat_plot_tol$Feature == feature_id[4]) |
           (dat_plot_tol$Feature == feature_id[5]) |
           (dat_plot_tol$Feature == feature_id[6])
           ,]
dat_plot$Estimation_Curve <- dat_plot$Curve * dat_plot$Estimation
dat_plot$Estimation <- dat_plot$Estimation * dat_plot$Point

dat_plot$Estimation[dat_plot$Method == "K-NN"][dat_plot$Estimation[dat_plot$Method == "K-NN"] < 0.5] <- NA

label <- (dat_plot$Method == "Smoothing spline") & (dat_plot$Feature == "Arterial Blood Pressure systolic")
dat_plot$Estimation_Curve[label][dat_plot$Estimation_Curve[label] < 0.5] <- NA

label <- (dat_plot$Method == "Smoothing spline") & (dat_plot$Feature == "Respiratory Rate")
dat_plot$Estimation_Curve[label][dat_plot$Estimation_Curve[label] > 1.6] <- NA

label <- (dat_plot$Method == "Matrix completion") & (dat_plot$Feature != "Base Excess")
dat_plot$Estimation[label][dat_plot$Estimation[label] < 0.4] <- NA

# dat_plot$Estimation[dat_plot$Method == "Matrix completion"][dat_plot$Estimation[dat_plot$Method == "Matrix completion"] < 0.4] <- NA

p_fit <- ggplot(dat_plot) + 
  geom_point(aes(x = Time, y = Value), color = "orange", size = 1.2) + 
  geom_line(aes(x = Time, y = Estimation_Curve), color = "blue", size = 0.6) +
  geom_point(aes(x = Time, y = Estimation), color = "blue", size = 0.6) + 
  facet_wrap(.~  Feature + Method, ncol = 4, scale = "free_y") +
  labs(x = "Time (min)", y = "Value",
       title = "(A)",
       colour = "", fill = "", linetype = "") + 
  # scale_color_manual(values = c("blue", rep("blue", 4))) +
  # scale_linetype_manual(values = c(1, 1)) + 
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        legend.position = "none",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  geom_blank(aes(y = x.value), data = data.frame(
    Feature = c(t(matrix(rep(unique(dat_plot$Feature), 2), ncol = 2))),
    x.value = c(rep(c(0.5, 1.6), 1),
                rep(c(0.5, 1.5), 1),
                rep(c(0.3, 1.6), 1),
                rep(c(0.5, 1.5), 1),
                rep(c(0.4, 1.6), 1),
                rep(c(-3, 3), 1)
                )
    ))
p_fit

ggsave(paste0("Figure/", "fit_dat_tol", ".pdf"), width = 9, height = 8, dpi = 300)

## Figure B
dat_plot <- dat_plot_tol[(dat_plot_tol$Feature == feature_id[7]) |
                           (dat_plot_tol$Feature == feature_id[8]) |
                           (dat_plot_tol$Feature == feature_id[9]) |
                           (dat_plot_tol$Feature == feature_id[10]) |
                           (dat_plot_tol$Feature == feature_id[11]) |
                           (dat_plot_tol$Feature == feature_id[12])
                         ,]
dat_plot$Estimation_Curve <- dat_plot$Curve * dat_plot$Estimation
dat_plot$Estimation <- dat_plot$Estimation * dat_plot$Point

# dat_plot$Estimation[dat_plot$Method == "K-NN"][dat_plot$Estimation[dat_plot$Method == "K-NN"] < 0.5] <- NA
# 
label <- (dat_plot$Method == "Smoothing spline") & (dat_plot$Feature == "Glucose")
dat_plot$Estimation_Curve[label][dat_plot$Estimation_Curve[label] < 0.05] <- NA

label <- (dat_plot$Method == "Smoothing spline") & (dat_plot$Feature == "Lactate")
dat_plot$Estimation_Curve[label][dat_plot$Estimation_Curve[label] < 0.4] <- NA

p_fit_2 <- ggplot(dat_plot) + 
  geom_point(aes(x = Time, y = Value), color = "orange", size = 1.2) + 
  geom_line(aes(x = Time, y = Estimation_Curve), color = "blue", size = 0.6) +
  geom_point(aes(x = Time, y = Estimation), color = "blue", size = 0.6) + 
  facet_wrap(.~  Feature + Method, ncol = 4, scale = "free_y") +
  labs(x = "Time (min)", y = "Value",
       title = "(B)",
       colour = "", fill = "", linetype = "") + 
  # scale_color_manual(values = c("blue", rep("blue", 4))) +
  # scale_linetype_manual(values = c(1, 1)) + 
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 10),
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  geom_blank(aes(y = x.value), data = data.frame(
    Feature = c(t(matrix(rep(unique(dat_plot$Feature), 2), ncol = 2))),
    x.value = c(rep(c(0.05, 1.6), 1),
                rep(c(0.2, 1.5), 1),
                rep(c(0.5, 2), 1),
                rep(c(0.2, 1.6), 1),
                rep(c(0.4, 1.5), 1),
                rep(c(0.2, 1.6), 1)
    )
  ))
p_fit_2

ggsave(paste0("Figure/", "fit_dat_tol", ".pdf"), width = 9, height = 8, dpi = 300)

grid.arrange(p_fit, p_fit_2, nrow = 1)
ppp <- arrangeGrob(p_fit, p_fit_2, nrow = 1) 

ggsave(paste0("Figure/", "fit_dat_tol", ".pdf"), ppp, width = 15, height = 9, dpi = 300)

## MSE
mean((dat_plot_mat$Value[is.finite(dat_plot_mat$Value) == T] - dat_plot_mat$Estimation[is.finite(dat_plot_mat$Value) == T]) ^ 2) 
mean((dat_plot_smo$Value[is.finite(dat_plot_fsvd$Value) == T] - dat_plot_smo$Estimation[is.finite(dat_plot_smo$Value) == T]) ^ 2) 
mean((dat_plot_knn$Value[is.finite(dat_plot_knn$Value) == T] - dat_plot_knn$Estimation[is.finite(dat_plot_knn$Value) == T]) ^ 2) 
mean((dat_plot_fsvd$Value[is.finite(dat_plot_fsvd$Value) == T] - dat_plot_fsvd$Estimation[is.finite(dat_plot_fsvd$Value) == T]) ^ 2) 

## Factor models
fit_FSVD$R <- 3
R <- fit_FSVD$R
fit_loa <- fit_FSVD$Loading[,1:R]
fit_fac <- fit_FSVD$Fac_serial[,1:R]

sig <- c(sign(t(fit_fac) %*% rep(1, nrow(fit_fac))))
fit_loa <- t(t(fit_loa) * sig)
fit_fac <- t(t(fit_fac) * sig)

dat_plot_fac <- data.frame(
  Time = rep((time_grid * (mean_t[2] - mean_t[1]) + mean_t[1]) / 60, fit_FSVD$R),
  Value = c(fit_fac),
  Component = c(sapply(1:fit_FSVD$R, function(i) rep(paste0("Factor ", i), length(time_grid))))
)

dat_plot_fac$Component <- as.factor(dat_plot_fac$Component)

p <- ggplot(dat_plot_fac) + 
  geom_line(aes(x = Time, y = Value), color = "#3182bd", size = 0.8) + 
  facet_wrap(.~Component, nrow = 1) +
  labs(x = "Time (min)", y = "Value",
       title = "",
       colour = "", fill = "", linetype = "") + 
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))  
p

# ggsave(paste0("code/Figure/", "fit_fac", ".pdf"), width = 10, height = 3, dpi = 300)

dat_plot_loa <- data.frame(
  Feature = rep(feature_id, fit_FSVD$R),
  Value = c(fit_loa),
  sig = c(fit_loa > 0),
  Component = c(sapply(1:fit_FSVD$R, function(i) rep(paste0("Factor ", i), n)))
)

dat_plot_loa$Component <- as.factor(dat_plot_loa$Component)

pp <- ggplot(dat_plot_loa) + 
  geom_bar(aes(x = Feature, y = Value, fill = sig), stat = "identity", width = 0.7) + 
  facet_wrap(.~Component, nrow = 1) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        # text = element_text(size = 10),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  # scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  scale_fill_manual(values = c("#D6AFB9", "#7E9BB7")) + 
  labs(x = "", y = "Factor loading",
       title = "",
       colour = "", fill = "", linetype = "") +
  coord_flip() +
  scale_x_discrete(limits = rev(feature_id))
pp

A <- matrix(c(NA, 1, 1, 1, 1, 1, 
              NA, 1, 1, 1, 1, 1, 
              NA, 1, 1, 1, 1, 1, 
              2, 2, 2, 2, 2, 2,
              2, 2, 2, 2, 2, 2,
              2, 2, 2, 2, 2, 2, 
              2, 2, 2, 2, 2, 2), byrow = T, nrow = 7)
grid.arrange(p, pp, ncol = 1, layout_matrix = A)
ppp <- arrangeGrob(p, pp, ncol = 1, layout_matrix = A) 
  
ggsave(paste0("Figure/", "fit_fac_loa", ".pdf"), ppp, width = 8.6, height = 6, dpi = 300)
