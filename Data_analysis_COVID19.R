# set working directory to be where the current script is located
rm(list=ls())
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Load functions and data
source("Function.R")
dat <- read.csv(file = "Data/time_series_covid19_confirmed_global.csv")
Population_dat <- read.csv(file = "Data/UID_ISO_FIPS_LookUp_Table.csv")
n <- nrow(dat)

# Data processing
region_name <- c("Luxembourg", "Iceland", "Ireland", "Switzerland", "Qatar",
                 "Belgium", "Spain", "Portugal", "Estonia", "Austria", "Italy",
                 "Norway", "Panama", "Netherlands", "Denmark", "Sweden", "Israel",
                 "Slovenia", "France", "Serbia", "Finland", "Bahrain", "Germany", 
                 "Czechia", "Croatia", "Chile", "Belarus", "United Kingdom", "Peru", "Iran",
                 "Romania", "Albania", "US", "Slovakia", "Canada", "Greece", "Saudi Arabia",
                 "Bulgaria", "Poland", "Kuwait", "United Arab Emirates", "Costa Rica", "Australia", "Russia",
                 "Georgia", "Lebanon", "Brazil", "Korea, South", "Singapore", "South Africa",
                 "Argentina", "Mexico", "China", "Algeria", "Pakistan", "Philippines",
                 "Malaysia", "Iraq", "Egypt", "Indonesia", "Taiwan*", "Thailand", "India", "Japan")

mark <- sapply(1:length(region_name), function(k){
  mark <- which(dat$Country.Region == region_name[k]) 
  if(length(mark) > 1){
    mark <- mark[which(dat$Province.State[mark] == "")]
    if(length(mark) == 0){
      mark <- which(dat$Country.Region == region_name[k])[1]
    }
  }
  return(mark)
})
dat <- dat[mark,]

Population <- sapply(1:n, function(i){
  Population_dat[which(Population_dat$Country_Region == region_name[i])[1],]$Population
})

## Data transformation
Ly <- list()
Lt <- list()
n <- nrow(dat)

for(i in 1:n){
  num <- unlist(dat[i,5:ncol(dat)])
  num <- num[is.na(num) == F]
  num <- log10(num[num >= 20] / Population[i] * 10 ^ 6)[1:67]
  
  Ly[[i]] <- c(unique(num))
  Lt[[i]] <- c(sapply(1:length(Ly[[i]]), function(k){
    min(which(Ly[[i]][k] == num))
  }) - 1)
}

mean_t <- c(min(unlist(Lt)), max(unlist(Lt)))

Lt <- lapply(1:n, function(i){
  (Lt[[i]] - mean_t[1]) / mean_t[2]
})

## Data plot
dat_plot <- cbind(value = unlist(Ly), time = unlist(Lt) * mean_t[2], Region = unlist(lapply(1:n, function(i){
  rep(i, length(Lt[[i]]))
})))
dat_plot <- as.data.frame(dat_plot)

select_coun <- c("Malaysia", "US", "Thailand", "Luxembourg", "Qatar", "Taiwan*")
mark <- sapply(1:length(select_coun), function(k) which(region_name == select_coun[k]))

p_1 <- ggplot() + 
  geom_line(data = dat_plot, aes(x = time, y = value, group = Region), color = "gray") +
  geom_line(aes(x = Lt[[mark[1]]] * mean_t[2], y = Ly[[mark[1]]])) +
  geom_point(aes(x = Lt[[mark[1]]] * mean_t[2], y = Ly[[mark[1]]], color = "Observed case count"), size = 0.5) +
  geom_text(aes(x = 35, y = 0.8, label = "Malaysia")) +
  geom_line(aes(x = Lt[[mark[2]]] * mean_t[2], y = Ly[[mark[2]]])) +
  geom_point(aes(x = Lt[[mark[2]]] * mean_t[2], y = Ly[[mark[2]]], color = "Observed case count"), size = 0.5) +
  geom_text(aes(x = 12, y = 0.8, label = "United States")) +
  geom_line(aes(x = Lt[[mark[3]]] * mean_t[2], y = Ly[[mark[3]]])) +
  geom_point(aes(x = Lt[[mark[3]]] * mean_t[2], y = Ly[[mark[3]]], color = "Observed case count"), size = 0.5) +
  geom_text(aes(x = 50, y = 0.3, label = "Thailand")) +
  geom_line(aes(x = Lt[[mark[4]]] * mean_t[2], y = Ly[[mark[4]]])) +
  geom_point(aes(x = Lt[[mark[4]]] * mean_t[2], y = Ly[[mark[4]]], color = "Observed case count"), size = 0.5) +
  geom_text(aes(x = 46, y = 4, label = "Luxembourg")) +
  # geom_line(aes(x = Lt[[mark[5]]] * mean_t[2], y = Ly[[mark[5]]])) +
  # geom_text(aes(x = 13, y = 2.6, label = "Qatar")) +
  geom_line(aes(x = Lt[[mark[6]]] * mean_t[2], y = Ly[[mark[6]]])) +
  geom_point(aes(x = Lt[[mark[6]]] * mean_t[2], y = Ly[[mark[6]]], color = "Observed case count"), size = 0.5) +
  geom_text(aes(x = 45, y = 2.2, label = "Taiwan")) +
  labs(x = "Day since 20+ cases", y = "Total casess per million (log scale)",
       title = "(A)",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        element_text(size = 30),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  scale_color_manual(values = "black") +
  ylim(c(-2, 5))
p_1

# CV for functional completion
time_grid <- seq(0, 1, length.out = 101)

time_mark <- lapply(1:n, function(i){
  sapply(1:length(Lt[[i]]), function(k) which.min(abs(time_grid - Lt[[i]][k])))
})

cv <- 5
sam_mark <- lapply(1:n, function(i){
  num <- length(Ly[[i]])
  g <- floor(num / cv)
  r <- num - g * cv
  if(r == 0){
    mark <- rep(1:cv, g)
  }else{
    mark <- c(rep(1:cv, g), 1:r)
  }
  return(mark)
})

CV_error <- sapply(1:cv, function(cv_num){
  Ly_tran <- lapply(1:n, function(i){
    Ly[[i]][sam_mark[[i]] != cv_num]
  })
  
  Lt_tran <- lapply(1:n, function(i){
    Lt[[i]][sam_mark[[i]] != cv_num]
  })
  
  fit_FPCA <- FPCA(Ly_tran, Lt_tran, optns = list(error = T, nRegGrid = 101, 
                                                  methodBwCov = "GCV",
                                                  methodMuCovEst = "smooth",
                                                  # FVEthreshold = 0.95,
                                                  methodBwCov = "GCV",
                                                  methodXi = "CE",
                                                  methodSelectK = 3
  ))
  
  fit_FSVD <- FSVD(Ly_tran, Lt_tran, R_max = 4, R_pre = 4, num_sel = "FD")
  
  ## Error
  err_fpca <-  sapply(1:n, function(i) 
    sum((((fit_FPCA$phi %*% fit_FPCA$xiEst[i,] + fit_FPCA$mu)[time_mark[[i]]])[sam_mark[[i]] == cv_num] - 
           Ly[[i]][sam_mark[[i]] == cv_num]) ^ 2)
  )
  
  err_fsvd <-  sapply(1:n, function(i) 
    sum(((((fit_FSVD$Intric_basis %*% (fit_FSVD$Score[i,]))[time_mark[[i]]])[sam_mark[[i]] == cv_num]) - 
           Ly[[i]][sam_mark[[i]] == cv_num]) ^ 2)
  )
  
  return(cbind(err_fpca = err_fpca,
               err_fsvd = err_fsvd))
}, simplify = "array")

error <- apply(CV_error, c(2), mean)

(error[1] - error[2]) / error[1]

## FSVD implementation
fit_FSVD <- FSVD(Ly, Lt, R_max = 4, R_pre = NULL, num_sel = "FD")

## Plot
### Intrinsic basis function
fit_FSVD$Intric_basis <- sapply(1:4, function(k)  fit_FSVD$Intric_basis[,k] * ifelse(fit_FSVD$Intric_basis[2,k] > fit_FSVD$Intric_basis[1,k], 1, -1))
dat_plot <- data.frame(
  time = rep(seq(0, 66, length.out = 101), 4),
  value = c(fit_FSVD$Intric_basis),
  mark = c(rep(paste0("1st IBF"), 101), 
           rep(paste0("2nd IBF"), 101),
           rep(paste0("3rd IBF"), 101),
           rep(paste0("4th IBF"), 101))
)

dat_plot$mark <- as.factor(dat_plot$mark)

p_2 <- ggplot(dat_plot) + 
  geom_line(aes(x = time, y = value, group = mark, color = mark), size = 0.7) +
  labs(x = "Day since 20+ cases", y = "Value of functions",
       title = "(B)",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  scale_x_continuous(labels = c(0, 20, 40, 60)) +
  # scale_x_continuous(breaks = seq(0, 1, length.out = 7)) +
  ylim(c(-2.5, 2.5)) +
  annotate("text", x = 0, y = 2, label = c(paste0("Cross-validation error: ", round(error[2], 3))), hjust = 0, vjust = 0, size = 3.5) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  # scale_linetype_manual(values = c(1, 2, 3, 4)) +
  scale_color_manual(values = c('red','#fc8d62','#8da0cb','#e78ac3')) 
p_2

# FPCA implementation
fit_FPCA <- FPCA(Ly, Lt, optns = list(error = T, nRegGrid = 101, 
                                      methodBwCov = "GCV",
                                      methodMuCovEst = "smooth",
                                      methodBwCov = "GCV",
                                      methodXi = "CE",
                                      methodSelectK = "AIC"
))

fit_FPCA$phi <- sapply(1:3, function(k)  fit_FPCA$phi[,k] * ifelse(fit_FPCA$phi[2,k] > fit_FPCA$phi[1,k], 1, -1))

dat_plot <- data.frame(
  time = rep(seq(0, 66, length.out = 101), 4),
  value = c(fit_FPCA$mu / sqrt(sum(fit_FPCA$mu ^ 2 * 0.01)), c(fit_FPCA$phi)),
  mark = c(rep(paste0(" MF"), 101), 
           rep(paste0("1st EF"), 101),
           rep(paste0("2nd EF"), 101),
           rep(paste0("3rd EF"), 101))
)

dat_plot$mark <- as.factor(dat_plot$mark)

p_3 <- ggplot(dat_plot) + 
  geom_line(aes(x = time, y = value, group = mark, color = mark), size = 0.7) +
  labs(x = "Day since 20+ cases", y = "Value of functions",
       title = "(C)",
       colour = "", fill = "", linetype = "") +
  annotate("text", x = 0, y = 2, label = c(paste0("Cross-validation error: ", round(error[1], 3))), hjust = 0, vjust = 0, size = 3.5) +
  theme_bw(base_family = "Times") +
  scale_x_continuous(labels = c(0, 20, 40, 60)) +
  ylim(c(-2.5, 2.5)) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  # scale_linetype_manual(values = c(1, 2, 3, 4)) +
  scale_color_manual(values = c('red','#fc8d62','#8da0cb','#e78ac3')) 
p_3

## Combine figures
A <- matrix(c(1, 2, 3,
              1, 2, 3,
              1, 2, 3,
              1, 2, 3,
              1, 2, 3,
              1, 2, 3
              ), byrow = T, nrow = 6)

gridExtra::grid.arrange(p_1, p_2, p_3, ncol = 3, layout_matrix = A)
p <- gridExtra::arrangeGrob(p_1, p_2, p_3, ncol = 3, layout_matrix = A) 

ggsave(paste0("Figure/", "curve_analysis", ".pdf"), p, width = 11, height = 3.5, dpi = 300)

# Functional clustering via FSVD
fit_clu_FSVD <- FClust_fsvd(Ly, Lt, R_max = 4, R_pre = NULL, Clu_num = F, abs = 0.001)

dat_plot <- data.frame(
  X = - fit_clu_FSVD$Score_EM[,1],
  Y = fit_clu_FSVD$Score_EM[,2],
  label = fit_clu_FSVD$fit_clu_FSVD_EM,
  Country = c("Luxembourg", "Iceland", "Ireland", "Switzerland", "Qatar",
              "Belgium", "Spain", "Portugal", "Estonia", "Austria", "Italy",
              "Norway", "Panama", "Netherlands", "Denmark", "Sweden", "Israel",
              "Slovenia", "France", "Serbia", "Finland", "Bahrain", "Germany", 
              "Czechia", "Croatia", "Chile", "Belarus", "United Kingdom", "Peru", "Iran",
              "Romania", "Albania", "United States", "Slovakia", "Canada", "Greece", "Saudi Arabia",
              "Bulgaria", "Poland", "Kuwait", "United Arab Emirates", "Costa Rica", "Australia", "Russia",
              "Georgia", "Lebanon", "Brazil", "Korea South", "Singapore", "South Africa",
              "Argentina", "Mexico", "China", "Algeria", "Pakistan", "Philippines",
              "Malaysia", "Iraq", "Egypt", "Indonesia", "Taiwan", "Thailand", "India", "Japan")
)

dat_plot$label <- factor(dat_plot$label)

p_1 <- ggplot(dat_plot) + 
  geom_point(aes(x = X, y = Y, color = label), size = 3) + 
  ggrepel::geom_text_repel(aes(x = X, y = Y, label = Country), size = 3,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 5)) + 
  labs(x = "Values of the 1st projected scores", y = "Values of the 2nd projected scores",
       title = "",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  # scale_x_continuous(labels = c(0, 20, 40, 60)) +
  # xlim(c(-1, 1)) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) + 
  scale_color_manual(values = c('#fc8d62', 'blue')) 
p_1
# ggsave(paste0("code/Figure/", "Map_clu", ".pdf"), width = 10, height = 10, dpi = 300)

dat_plot <- data.frame(
  X = dat$Long,
  Y = dat$Lat,
  label = fit_clu_FSVD$fit_clu_FSVD_EM,
  Country = c("Luxembourg", "Iceland", "Ireland", "Switzerland", "Qatar",
              "Belgium", "Spain", "Portugal", "Estonia", "Austria", "Italy",
              "Norway", "Panama", "Netherlands", "Denmark", "Sweden", "Israel",
              "Slovenia", "France", "Serbia", "Finland", "Bahrain", "Germany", 
              "Czechia", "Croatia", "Chile", "Belarus", "United Kingdom", "Peru", "Iran",
              "Romania", "Albania", "United States", "Slovakia", "Canada", "Greece", "Saudi Arabia",
              "Bulgaria", "Poland", "Kuwait", "United Arab Emirates", "Costa Rica", "Australia", "Russia",
              "Georgia", "Lebanon", "Brazil", "Korea South", "Singapore", "South Africa",
              "Argentina", "Mexico", "China", "Algeria", "Pakistan", "Philippines",
              "Malaysia", "Iraq", "Egypt", "Indonesia", "Taiwan", "Thailand", "India", "Japan")
)

dat_plot$label <- factor(dat_plot$label)
library(maps)
world <- map_data("world")

p_2 <- ggplot() + 
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray") + 
  # coord_fixed(1.3) +
  geom_point(data = dat_plot, aes(x = X, y = Y, color = label), size = 1.3) + 
  ggrepel::geom_text_repel(data = dat_plot, aes(x = X, y = Y, label = Country),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 10), size = 2.6) +
  labs(x = "", y = "",
       title = "(A)",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  # scale_x_continuous(labels = c(0, 20, 40, 60)) +
  ylim(c(-60, 100)) +
  xlim(c(-150, 190)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none",
        axis.line=element_blank(),
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  scale_color_manual(values = c('#fc8d62', 'blue')) 
p_2

ggsave(paste0("Figure/", "Map_clu", ".pdf"), p_2, width = 7, height = 5, dpi = 300)

## Clusters' mean function
time_grid <- seq(0, 1, length.out = 101)

mean_1 <- fit_clu_FSVD$fit_FSVD$Intric_basis %*% fit_clu_FSVD$mu[,1]
mean_2 <- fit_clu_FSVD$fit_FSVD$Intric_basis %*% fit_clu_FSVD$mu[,2]

dat_plot <- data.frame(time = rep(time_grid, 2) * 67, 
                       val = c(mean_1, mean_2),
                       label = c(rep("Cluster 2", length(time_grid)), rep("Cluster 1", length(time_grid)))
)

p_3 <- ggplot(dat_plot) + 
  geom_line(aes(x = time, y = val, color = label), size = 1) +
  labs(x = "Day since 20+ cases", y = "Value of functions",
       title = "(B)",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  scale_x_continuous(breaks = c(0, 20, 40, 60)) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  scale_color_manual(values = c('blue', '#fc8d62')) 
p_3

gridExtra::grid.arrange(p_2, p_3, ncol = 2, layout_matrix = matrix(c(1, 1, 1, 1, 2, 2, 2, 2, 2), nrow = 1))
p <- gridExtra::grid.arrange(p_2, p_3, ncol = 2, layout_matrix = matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2), nrow = 1))

ggsave(paste0("Figure/", "clu_map_mean", ".pdf"), p, width = 10, height = 3.5, dpi = 300)


