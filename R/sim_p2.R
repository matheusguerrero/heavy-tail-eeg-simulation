rm(list = ls())


source("./R/set_packages.R")
source("./R/fun_main.R")
source("./R/fun_plot.R")


#### Common parameters ####  
sample_size = 10
sampling_rate = 256
window_duration = 1


#### Patient 02 #### 

# Setting up patient parameters (patient scenario)
patient_02 <- list(
  pt_id = 2,
  seizure = 1,
  sample_size = sample_size,
  sampling_rate = sampling_rate,
  window_duration = window_duration,
  n_channels = 4,
  n_clusters = 3,
  clusters = c(1, 2, 2, 3),
  fixed_weights = list(
    W_1 = get_fixed_weights(
      before = c(0.00, 0.50, 3.00, 0.50, 2.00),
      during = c(1.00, 1.00, 8.00, 0.50, 12.0),
      after  = c(1.50, 0.00, 2.00, 0.00, 2.00),
      row_sumupto1 = FALSE
    ),
    W_2 = get_fixed_weights(
      before = c(0.50, 0.50, 2.00, 0.50, 2.50),
      during = c(1.00, 0.50, 6.00, 1.00, 10.0),
      after  = c(0.00, 0.00, 1.00, 0.00, 2.00),
      row_sumupto1 = FALSE
    ),
    W_3 = get_fixed_weights(
      before = c(2.00, 0.50, 1.00, 0.50, 2.00),
      during = c(0.50, 0.50, 3.00, 0.00, 11.0),
      after =  c(0.00, 0.50, 1.00, 0.50, 1.00),
      row_sumupto1 = FALSE
    )
  ),
  seizure_start = c(2.5, 2.6, 2.6),
  seizure_end   = c(6.0, 6.2, 5.9),
  error_covariance = rbind(
    c(NA, 0.10, 0.10, 0.10),
    c(NA,   NA, 0.85, 0.10), 
    c(NA,   NA,   NA, 0.20),  
    c(NA,   NA,   NA,   NA)
  ),
  error_variances = rep(1, 4),
  tail_type = "exponential",
  tail_par = c(1.1, 2.2, 0.5),
  tail_band_weight = c(0.3, 0.9, 1.1, 0.75, 1.7),
  tail_factor_place = "out_z"
)  

# Simulating signals
p2 <- sim_patient_eeg(
  patient_02$pt_id, 
  patient_02$seizure, 
  patient_02$sample_size,
  patient_02$sampling_rate,
  patient_02$window_duration,
  patient_02$n_channels,
  patient_02$n_clusters, 
  patient_02$clusters, 
  patient_02$fixed_weights, 
  patient_02$seizure_start, 
  patient_02$seizure_end, 
  patient_02$error_covariance,
  patient_02$error_variances,
  patient_02$tail_type, 
  patient_02$tail_par,
  patient_02$tail_band_weight, 
  patient_02$tail_factor_place
)  


# Checking EEGs
rbind(head(p2$eeg, n = 3), tail(p2$eeg, n = 3))


# Plotting signals and periodogram
plot_signal(p2)


# Decomposing signals into the five canonical frequency bands and plotting
p2_bands <- get_eeg_bands(p2)
bands <- c("delta", "theta", "alpha", "beta", "gamma")
for (i in 1:5) {
  print(plot_signal(p2, p2_bands[[i]], bands[i])) 
}


# Checking tail dependence via chi measure
plot_chi_win(p2, "chi", upp_quantile = .925)


# Study of the GPD parameters
p2_gpd <- get_gpd(p2)
plot_gpd_shape(p2_gpd)
plot_gpd_scale(p2_gpd)


# We can do everything at once. Note: it may take a while.
p2_full <- sim_patient(patient_02)