# set working directory to be where the current script is located
rm(list=ls())
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# Load main function
source("Function.R")

###############################################################################
# Setting 1
## Intrinsic basis functions' case
### Heterogeneous cases
basis_num <- 3 # Number of basis functions used
rat <- 0.05 # Noise level

sample_mark <- c(50, 100, 150) # Sample size
point_mark <- c(6, 8, 10) # Mean number of observed time points
for(n in sample_mark){
  for(obs_point in point_mark){
    source("Simulation_fun.R")
  }
}

## Estimation errors of basis functions
bas_error <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    NULL
  })
})

for(i in 1:3){
  for(k in 1:3){
    n <- sample_mark[i]
    obs_point <- point_mark[k]
    load(paste0("Result/Result_", n, "_", obs_point, ".rda"))
    bas_error[[i]][[k]] <- sapply(1:basis_num, function(k) basis_error(Result, k))
  }
}

lapply(1:3, function(i){
  a <- cbind(bas_error[[i]][[1]], bas_error[[i]][[2]], bas_error[[i]][[3]])
  rownames(a) <- c("FPCA(0)", "FPCA", "FSVD")
  a <- a[c(2,1,3),]
  a <- xtable(a)
  return(a)
})

## Estimation error of functional completion
com_error <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    NULL
  })
})

for(i in 1:3){
  for(k in 1:3){
    n <- sample_mark[i]
    obs_point <- point_mark[k]
    load(paste0("Result/Result_", n, "_", obs_point, ".rda"))
    com_error[[i]][[k]] <- Completion_error(Result)
  }
}

lapply(1:3, function(i){
  a <- cbind(com_error[[i]][[1]], com_error[[i]][[2]], com_error[[i]][[3]])
  rownames(a) <- c("Smooth spline", "FPCA(0)", "FPCA", "FSVD")
  a <- a[c(3,2,1,4),]
  a <- xtable(a)
  return(a)
})

### i.i.d. cases
basis_num <- 3 # Number of basis functions used
rat <- 0.05 # Noise level

sample_mark <- c(50, 100, 150) # Sample size
point_mark <- c(6, 8, 10) # Mean number of observed time points
for(n in sample_mark){
  for(obs_point in point_mark){
    source("Simulation_fun_iid.R")
  }
}

## Estimation errors of basis functions
bas_error <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    NULL
  })
})

for(i in 1:3){
  for(k in 1:3){
    n <- sample_mark[i]
    obs_point <- point_mark[k]
    load(paste0("Result/Result_iid_", n, "_", obs_point, ".rda"))
    bas_error[[i]][[k]] <- sapply(1:basis_num, function(k) basis_error(Result, k))
  }
}

lapply(1:3, function(i){
  a <- cbind(bas_error[[i]][[1]], bas_error[[i]][[2]], bas_error[[i]][[3]])
  rownames(a) <- c("FPCA(0)", "FPCA", "FSVD")
  a <- a[c(2,1,3),]
  a <- xtable(a)
  return(a)
})

## Estimation error of functional completion
com_error <- lapply(1:3, function(i){
  lapply(1:3, function(k){
    NULL
  })
})

for(i in 1:3){
  for(k in 1:3){
    n <- sample_mark[i]
    obs_point <- point_mark[k]
    load(paste0("Result/Result_iid_", n, "_", obs_point, ".rda"))
    com_error[[i]][[k]] <- Completion_error(Result)
  }
}

lapply(1:3, function(i){
  a <- cbind(com_error[[i]][[1]], com_error[[i]][[2]], com_error[[i]][[3]])
  rownames(a) <- c("Smooth spline", "FPCA(0)", "FPCA", "FSVD")
  a <- a[c(3,2,1,4),]
  a <- xtable(a)
  return(a)
})

################################################################################
# Setting 2
## Functional clustering case
basis_num <- 3
rat <- 0.05

sample_mark <- c(50, 100, 150)
point_mark <- c(6, 8, 10)
for(n in sample_mark){
  for(obs_point in point_mark){
    source("Simulation_clu.R")
  }
}

## Estimation error of functional clustering
clu_Error <- NULL

for(i in 1:3){
  for(k in 1:3){
    n <- sample_mark[i]
    obs_point <- point_mark[k]
    load(paste0("Result/Result_", n, "_", obs_point, "_clu", ".rda"))
    res <- clu_error(Result)
    clu_Error <- rbind(clu_Error, 
                       data.frame(value = unlist(res),
                                  Method = rep(c(rep("Smoothing-clustering", 100), rep("FPCA-clustering", 100),
                                                 rep("FSVD-clustering", 100), rep("FSVD-EM-clustering", 100)), 1),
                                  n = rep(n, 400),
                                  p = rep(paste0("{", obs_point - 2, ",...,", obs_point + 2, "}"), 400)))
  }
}

clu_Error$Method <- factor(clu_Error$Method, levels = c("Smoothing-clustering",
                                                        "FPCA-clustering",
                                                        "FSVD-clustering",
                                                        "FSVD-EM-clustering"))

ggplot(clu_Error) +
  geom_boxplot(aes(x = Method, y = value, fill = Method),
               position = "dodge", size = 0.3) +
  facet_wrap(n ~  p) +
  labs(x = "Method", y = "ARI",
       title = "",
       colour = "", fill = "", linetype = "") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        axis.text.x=element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5)) +
  # scale_fill_manual(values = c("blue", "orange")) 
  scale_fill_manual(values = c("#2c7fb8", "#fdae61", "red", "#b30000")) +
  scale_color_manual(values = c("#2c7fb8", "#fdae61", "red", "#b30000"))
# scale_color_manual(values = c("blue", "orange", "red", "#b30000"))

ggsave(paste0("Figure/", "sim_clu", ".pdf"), width = 6.5, height = 5.5, dpi = 300)

###############################################################################
# Setting 3
## Model setting
basis_num <- 3
rat <- 0.05

sample_mark <- c(50, 100, 150)
point_mark <- c(6, 10, 14)

for(n in sample_mark){
  for(obs_point in point_mark){
    source("Simulation_fac.R")
  }
}

dat_plot <- NULL
for(i in 1:3){
  for(k in 1:3){
    n <- sample_mark[i]
    obs_point <- point_mark[k]
    load(paste0("Result/Result_fac_", n, "_", obs_point, ".rda"))
    p <- data.frame(error = unlist(load_error(Result)),
                    Method = c(" FAM", "  SVD", "FSVD"),
                    n = rep(n, 3),
                    J = rep(paste0("{", obs_point - 2, ",...,", obs_point + 2, "}"), 3))
    dat_plot <- rbind(dat_plot, p)
  }
}

colnames(dat_plot) <- c("NMSE", "Method", "n", "J")
dat_plot$J <- factor(dat_plot$J, level = c("{4,...,8}", "{8,...,12}", "{12,...,16}"))

ggplot(dat_plot) + 
  geom_bar(aes(x = Method, y = NMSE, fill = Method),
           stat = "identity",
           position = "dodge",
           width = 0.7
           # alpha = 0.8
  ) +
  facet_wrap(n ~  J, nrow = 3) +
  scale_fill_manual(values = c("#F0A780", "#96B6D8", "#e34a33")) +
  labs(x = "Method", y = "NMSE (%)",
       title = "",
       colour = "", fill = "", linetype = "") +
  # scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = c("0%", "25%", "50%",
  #                                                                     "75%", "100%")) +
  # scale_x_discrete(breaks = c("0", "0.15", "0.3"), labels = c("0%", "15%", "30%")) +
  theme_bw() +
  # ggthemes::theme_tufte() +
  theme(text=element_text(size=15),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        axis.text.x=element_blank(),
        # text = element_text(family = "STHeiti"),
        plot.title = element_text(hjust = 0.5)) 

ggsave(paste0("Figure/", "sim_fac_error", ".pdf"), width = 6, height = 5.5, dpi = 300)
