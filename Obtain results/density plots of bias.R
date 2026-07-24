#CLPM
# 1. Load Required Libraries
require(foreach)
require(ggplot2)
require(tidyr)
require(dplyr)

# 2. Define True Values and Target Columns
true_value_FT <- c(.3,.3,.1,.1,1.1383929,1.1383929,.4092262,1,1,.3)
true_value_MVN <- c(.3,.3,.1,.1,1.1383929,1.1383929,.4092262,1,1,.3)
param_cols <- c("XXest", "YYest", "XYest", "YXest", "VXest", "VYest", "CXYest", "VeXest", "VeYest", "CeXYest")

# 3. Read Data and Calculate Bias
raw_bias_list <- list()

for (TYPE in c("MVN", "FT")) {
  # Get the true values for the current TYPE
  true_vals <- get(paste0("true_value_", TYPE))
  
  for (i in c(1, 2, 4, 8)) {
    # Read the data
    filepath <- (paste0("~/Library/CloudStorage/OneDrive-UCB-O365/research/project/Simulation comparison/Power results/CLPM/",TYPE,"_Estimates_", i,".txt"))
    dat <- read.csv(filepath)
    
    # Subset and rename columns 
    dat_S<-dat[,2:ncol(dat)]
    names(dat_S)<-c("CLPM_chisq","XXest", "YYest", "XYest", "YXest", "VXest","VYest", "CXYest", "VeXest", "VeYest", "CeXYest","CLPMnoXX_chisq","XXest_noXX", "YYest_noXX", "XYest_noXX", "YXest_noXX", "VXest_noXX","VYest_noXX", "CXYest_noXX", "VeXest_noXX", "VeYest_noXX", "CeXYest_noXX","CLPMnoYY_chisq","XXest_noYY", "YYest_noYY", "XYest_noYY", "YXest_noYY", "VXest_noYY","VYest_noYY", "CXYest_noYY", "VeXest_noYY", "VeYest_noYY", "CeXYest_noYY","CLPMnoXY_chisq","XXest_noXY", "YYest_noXY", "XYest_noXY", "YXest_noXY", "VXest_noXY","VYest_noXY", "CXYest_noXY", "VeXest_noXY", "VeYest_noXY", "CeXYest_noXY","CLPMnoYX_chisq","XX_noYXest", "YYest_noYX", "XYest_noYX", "YXest_noYX", "VXest_noYX","VYest_noYX", "CXYest_noYX", "VeXest_noYX", "VeYest_noYX", "CeXYest_noYX","CLPMnoCXY_chisq","XXest_noCXY", "YYest_noCXY", "XYest_noCXY", "YXest_noCXY", "VXest_noCXY","VYest_noCXY", "CXYest_noCXY", "VeXest_noCXY", "VeYest_noCXY", "CeXYest_noCXY","CLPMnoCeXY_chisq","XXest_noCeXY", "YYest_noCeXY", "XYest_noCeXY", "YXest_noCeXY", "VXest_noCeXY","VYest_noCeXY", "CXYest_noCeXY", "VeXest_noCeXY", "VeYest_noCeXY", "CeXYest_noCeXY")
    dat2 <- as.data.frame(lapply(dat_S, function(x) as.numeric(as.character(x))))
    
    # Extract estimates and calculate bias (Estimate - True Value)
    estimates_df <- dat2[, param_cols]
    bias_df <- sweep(estimates_df, 2, true_vals, "-")
    
    # Tag with TYPE and File_I, then add to list
    bias_df$TYPE <- TYPE
    bias_df$File_I <- i
    raw_bias_list[[length(raw_bias_list) + 1]] <- bias_df
  }
}

# Combine all files into one dataframe
all_bias_data <- bind_rows(raw_bias_list)

# 4. Reshape and Reorder Data
# Pivot to long format for ggplot
long_bias_data <- pivot_longer(all_bias_data, 
                               cols = all_of(param_cols), 
                               names_to = "Parameter", 
                               values_to = "Bias")

# Set the factor levels to switch the order of XYest and YXest in the plot
new_param_order <- c("XXest", "YYest", "YXest", "XYest", "VXest", "VYest", "CXYest", "VeXest", "VeYest", "CeXYest")
long_bias_data$Parameter <- factor(long_bias_data$Parameter, levels = new_param_order)

# 5. Define Plotmath Labels for the Facets (Updated for Bold)
param_labels <- c(
  "XXest"   = "bold(a[xx])",
  "YYest"   = "bold(a[yy])",
  "XYest"   = "bold(a[xy])",
  "YXest"   = "bold(a[yx])",
  "VXest"   = "bold(V[x])",
  "VYest"   = "bold(V[y])",
  "CXYest"  = "bold(Cov[xy])",
  "VeXest"  = "bold(Re[x])",
  "VeYest"  = "bold(Re[y])",
  "CeXYest" = "bold(CovRe[xy])"
)

# 6. Generate the Plot
bias_density_plot <- ggplot(long_bias_data, aes(x = Bias, fill = TYPE, color = TYPE)) +
  geom_density(alpha = 0.5) + 
  
  # Apply labels and 3-column layout
  facet_wrap(~ Parameter, scales = "free", ncol = 3, 
             labeller = as_labeller(param_labels, default = label_parsed)) + 
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14),       # Facet labels (Parameter names)
    legend.text = element_text(size = 12),      # Legend labels (process-based, etc.)
    legend.title = element_text(size = 14)      # Legend title (Simulation)
  ) +
  # Set titles and axis labels
  labs(title = "",
       x = "",
       y = "Density",
       fill = "Simulation",
       color = "Simulation") + 
  
  # Customize legend: order, colors, and names
  scale_fill_manual(
    values = c("FT" = "#ff7f0e", "MVN" = "#1f77b4"),
    breaks = c("FT", "MVN"), 
    labels = c("FT" = "process-based", "MVN" = "covariance-based")
  ) +
  scale_color_manual(
    values = c("FT" = "#b45f06", "MVN" = "#0b5394"),
    breaks = c("FT", "MVN"),
    labels = c("FT" = "process-based", "MVN" = "covariance-based")
  ) +
  
  # Add reference line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6)

# 7. Print the Plot
print(bias_density_plot)

#CoT
# 1. Load Required Libraries
require(foreach)
require(ggplot2)
require(tidyr)
require(dplyr)

# 2. Define True Values and Target Columns
true_value_GE <- c(.212042, .2000046, .3151295, .2000117, .1999925, .1094478)
true_value_MVN <- c(.212175, .2, .3150366, .2, .2, 0.1095111)
param_cols <- c("VAest","Dest","Xest","Sest","Eest","MUest")
# 3. Read Data and Calculate Bias
raw_bias_list <- list()

for (TYPE in c("MVN", "GE")) {
  # Get the true values for the current TYPE
  true_vals <- get(paste0("true_value_", TYPE))
  
  for (i in c(1,2,4,8)) {
    # Read the data
    file_suffix <- ifelse(TYPE == "GE", "_Estimates_noabound_nomissing_", "_Estimates_noabound_")
    filepath <- (paste0("~/Library/CloudStorage/OneDrive-UCB-O365/research/project/Simulation comparison/Power results/UniCoTP/", TYPE, file_suffix, i))
    dat <- read.csv(filepath)
    
    # Subset and rename columns exactly as you had them
    dat_S<-dat[,2:ncol(dat)]
    names(dat_S)<-c(
      "CoTP_ll","CoTP_status","Aest","Sest","Eest","Mest","Dest","Test","MUest","Qest","Xest","West",
      "CoTPnoA_ll","CoTPnoA_status","Aest_noA","Sest_noA","Eest_noA","Mest_noA","Dest_noA","Test_noA","MUest_noA",
      "CoTPnoS_ll","CoTPnoS_status","Aest_noS","Sest_noS","Eest_noS","Mest_noS","Dest_noS","Test_noS","MUest_noS",
      "CoTPnoM_ll","CoTPnoM_status","Aest_noM","Sest_noM","Eest_noM","Mest_noM","Dest_noM","Test_noM","MUest_noM",
      "CoTPnoD_ll","CoTPnoD_status","Aest_noD","Sest_noD","Eest_noD","Mest_noD","Dest_noD","Test_noD","MUest_noD",
      "CoTPnoE_ll","CoTPnoE_status","Aest_noE","Sest_noE","Eest_noE","Mest_noE","Dest_noE","Test_noE","MUest_noE",
      "CoTPnoMU_ll","CoTPnoMU_status","Aest_noMU","Sest_noMU","Eest_noMU","Mest_noMU","Dest_noMU","Test_noMU","MUest_noMU"
    )
    
    dat2a<-subset(dat_S,CoTP_status<=1)
    dat2<-subset(dat2a,round(CoTP_ll,10)<=(round(CoTPnoA_ll,10)+.00001)|CoTPnoA_status>1)
    dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoS_ll,10)+.00001)|CoTPnoS_status>1)
    dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoM_ll,10)+.00001)|CoTPnoM_status>1)
    dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoD_ll,10)+.00001)|CoTPnoD_status>1)
    dat2<-subset(dat2,round(CoTP_ll,10)<=(round(CoTPnoMU_ll,10)+.00001)|CoTPnoMU_status>1)
    dat2 <- as.data.frame(lapply(dat2, function(x) as.numeric(as.character(x))))
    dat2$VAest<-dat2$Aest*dat2$Aest*dat2$Qest
    
    # Extract estimates and calculate bias (Estimate - True Value)
    estimates_df <- dat2[, param_cols]
    bias_df <- sweep(estimates_df, 2, true_vals, "-")
    
    # Tag with TYPE and File_I, then add to list
    bias_df$TYPE <- TYPE
    bias_df$File_I <- i
    raw_bias_list[[length(raw_bias_list) + 1]] <- bias_df
  }
}

# Combine all files into one dataframe
all_bias_data <- bind_rows(raw_bias_list)

# 4. Reshape and Reorder Data
# Pivot to long format for ggplot
long_bias_data <- pivot_longer(all_bias_data, 
                               cols = all_of(param_cols), 
                               names_to = "Parameter", 
                               values_to = "Bias")

# Set the factor levels
long_bias_data$Parameter <- factor(long_bias_data$Parameter, levels = param_cols)

# 5. Define Plotmath Labels for the Facets (Updated for Bold)
param_labels <- c(
  "VAest"   = "bold(V[A])",
  "Dest"    = "bold(V[D])",
  "Xest"    = "bold(V[F])",
  "Sest"    = "bold(V[S])",
  "Eest"    = "bold(V[E])",
  "MUest"   = "bold(mu)"
)

# 6. Generate the Plot
bias_density_plot <- ggplot(long_bias_data, aes(x = Bias, fill = TYPE, color = TYPE)) +
  geom_density(alpha = 0.5) + 
  
  # Apply labels and 3-column layout
  facet_wrap(~ Parameter, scales = "free", ncol = 3, 
             labeller = as_labeller(param_labels, default = label_parsed)) + 
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14),       # Facet labels (Parameter names)
    legend.text = element_text(size = 12),      # Legend labels (process-based, etc.)
    legend.title = element_text(size = 14)      # Legend title (Simulation)
  ) +
  # Set titles and axis labels
  labs(title = "",
       x = "",
       y = "Density",
       fill = "Simulation",
       color = "Simulation") + 
  
  # Customize legend: order, colors, and names
  scale_fill_manual(
    values = c("GE" = "#ff7f0e", "MVN" = "#1f77b4"),
    breaks = c("GE", "MVN"), 
    labels = c("GE" = "process-based", "MVN" = "covariance-based")
  ) +
  scale_color_manual(
    values = c("GE" = "#b45f06", "MVN" = "#0b5394"),
    breaks = c("GE", "MVN"),
    labels = c("GE" = "process-based", "MVN" = "covariance-based")
  ) +
  
  # Add reference line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6)

# 7. Print the Plot
print(bias_density_plot)
