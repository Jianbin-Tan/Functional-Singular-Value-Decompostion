# set working directory to be where the current script is located
rm(list=ls())
mydir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(mydir)

# COVID-19 Data
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

continent <- c(
  "Europe", "Europe", "Europe", "Europe", "Asia", "Europe", "Europe", "Europe",
  "Europe", "Europe", "Europe", "Europe", "North America", "Europe", "Europe", "Europe",
  "Asia", "Europe", "Europe", "Europe", "Europe", "Asia", "Europe", "Europe",
  "Europe", "South America", "Europe", "Europe", "South America", "Asia", "Europe", "Europe",
  "North America", "Europe", "North America", "Europe", "Asia", "Europe", "Europe", "Asia",
  "Asia", "North America", "Oceania", "Europe", "Asia", "Asia", "South America", "Asia",
  "Asia", "Africa", "South America", "North America", "Asia", "Africa", "Asia", "Asia",
  "Asia", "Asia", "Africa", "Asia", "Asia", "Asia", "Asia", "Asia"
)

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
})),
cont = unlist(lapply(1:n, function(i){
  rep(continent[i], length(Lt[[i]]))
}))
)
dat_plot <- as.data.frame(dat_plot)
dat_plot$value <- as.numeric(dat_plot$value)
dat_plot$time <- as.numeric(dat_plot$time)

p_1 <- ggplot(dat_plot) + 
  geom_line(aes(x = time, y = value, group = Region), size = 0.3) +
  geom_vline(xintercept = seq(0, 1, length.out = 13) * 66, linetype = 2) + 
  geom_point(aes(x = time, y = value, group = Region, color = cont), size = 1.2) +
  # facet_wrap(.~cont, scales = "free_y", ncol = 3) +
  labs(x = "Day since 20+ cases", y = "Total cases per million (log scale)",
       title = "(A) Regional COVID-19 Case Trajectories",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),  # <-- This removes all major grid lines
        legend.position = "none",
        panel.border = element_blank(),
        element_text(size = 30),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) 
  # scale_x_continuous(breaks = floor(seq(0, 66, length.out = 13))) + 
  # scale_color_manual(values = "black") +
  # ylim(c(-2, 5))
p_1

# Data exploration
cont <- unique(continent)
time_group <- seq(0, 1, length.out = 13)
data_plot <- data.frame()

for(i in 1:length(cont)){
  for(j in 1:12){
    mark <- which(continent == cont[i])
    sample <- unlist(lapply(mark, function(k) Ly[[k]][(Lt[[k]] >= time_group[j]) & (Lt[[k]] < time_group[j + 1])]))
    
    data_plot <- rbind(data_plot,
                       data.frame(value = sample,
                                  cont = rep(cont[i], length(sample)),
                                  time = rep(paste0("Time period ", j), length(sample)))
                       )
  }
}

data_plot$time <- factor(data_plot$time, levels = sapply(1:12, function(j) paste0("Time period ", j)))

p_2 <- ggplot(data_plot) + 
  geom_boxplot(aes(x = cont, y = value, color = cont), outliers = FALSE) +
  facet_wrap(.~time, scales = "free_y", ncol = 4) +
  labs(x = NULL, y = NULL,
       title = "(C) Samples from Different Continents",
       colour = "", fill = "", linetype = "") +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),  # <-- removes major grid lines
        legend.position = "top",
        panel.border = element_blank(),
        element_text(size = 30),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +   # <-- removes x-axis tick labels
  guides(color = guide_legend(nrow = 1))
p_2

# EHR data
load("Data/dat_ehr.rda")

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

feature_id[4] <- "Blood Pressure"
dat_plot <- data.frame(Time = (unlist(Lt) * (mean_t[2] - mean_t[1]) + mean_t[1]) / 60,
                       Value = unlist(Ly), Feature = c(unlist(lapply(1:n, function(i) rep(feature_id[i], length(Ly[[i]]))))))

p_3 <- ggplot(dat_plot) + 
  geom_line(aes(x = Time, y = Value, group = Feature), size = 0.3) +
  geom_vline(xintercept = seq(0, 1, length.out = 4) * 600, linetype = 2) + 
  geom_point(aes(x = Time, y = Value, group = Feature, color = Feature), size = 1.2) +
  # facet_wrap(.~Feature, scales = "free_y", ncol = 4) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),  # <-- Removes major grid lines
        legend.position = "none",
        panel.border = element_blank(),
        element_text(size = 30),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)) +
  labs(x = "Time (min)", y = "Value of Clinical Features",
       title = "(B) Clinical Longitudinal Data of a Patient",
       colour = "", fill = "", linetype = "")
p_3

# Data exploration
dat_plot <- data.frame()
for(i in 1:length(Ly)){
  num_1 <- length(Ly[[i]][Lt[[i]] <= 1/3])
  num_3 <- length(Ly[[i]][Lt[[i]] > 2/3])
  num_2 <- length(Ly[[i]]) - num_1 - num_3
  dat_plot <- rbind(dat_plot, rbind(data.frame(Time = rep("Time period 1", num_1), Value = Ly[[i]][Lt[[i]] <= 1/3], Feature =  rep(feature_id[i], num_1)),
                                    data.frame(Time = rep("Time period 2", num_2), Value = Ly[[i]][(Lt[[i]] > 1/3)&(Lt[[i]] <= 2/3)], Feature =  rep(feature_id[i], num_2)),
                                    data.frame(Time = rep("Time period 3", num_3), Value = Ly[[i]][Lt[[i]] > 2/3], Feature =  rep(feature_id[i], num_3))
                                    ))
}

p_4 <- ggplot(dat_plot) + 
  geom_boxplot(aes(x = Time, y = Value, color = Time), outliers = FALSE) +
  # geom_point(aes(x = Time, y = Value, group = Feature), size = 0.5) + 
  facet_wrap(.~Feature, scales = "free_y", ncol = 4) +
  theme_bw(base_family = "Times") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),  # <-- This removes major grid lines
        legend.position = "top",
        panel.border = element_blank(),
        element_text(size = 30),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "Time period", y = "",
       title = "(D) Samples from Different Time Periods",
       colour = "", fill = "", linetype = "")
p_4

p <- grid.arrange(p_1, p_3, p_2, p_4, ncol= 2, layout_matrix = matrix(c(1, 1, 1, 3, 3, 3, 3, 2, 2, 2, 4, 4, 4, 4), 7, 2))

ggsave(paste0("Figure/", "data_combine", ".pdf"), p, width = 13, height = 8, dpi = 300)
