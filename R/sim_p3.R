rm(list = ls())


source("./R/set_packages.R")
source("./R/fun_main.R")
source("./R/fun_plot.R")


#### Common parameters ####  
sample_size = 10
sampling_rate = 256
window_duration = 1


#### Patient 03 #### 

# Setting up patient parameters (patient scenario)
patient_03 <- list(
  pt_id = 3,
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
  tail_band_weight = c(0, 0, 1, 0, 1),
  tail_factor_place = "in_z"
)  

# Simulating signals
p3 <- sim_patient_eeg(
  patient_03$pt_id, 
  patient_03$seizure, 
  patient_03$sample_size,
  patient_03$sampling_rate,
  patient_03$window_duration,
  patient_03$n_channels,
  patient_03$n_clusters, 
  patient_03$clusters, 
  patient_03$fixed_weights, 
  patient_03$seizure_start, 
  patient_03$seizure_end, 
  patient_03$error_covariance,
  patient_03$error_variances,
  patient_03$tail_type, 
  patient_03$tail_par,
  patient_03$tail_band_weight, 
  patient_03$tail_factor_place
)  


# Checking EEGs
rbind(head(p3$eeg, n = 3), tail(p3$eeg, n = 3))


# Plotting signals and periodogram
plot_signal(p3)


# Decomposing signals into the five canonical frequency bands and plotting
p3_bands <- get_eeg_bands(p3)
band <- c("delta", "theta", "alpha", "beta", "gamma")
for (i in 1:5) {
  print(plot_signal(p3, p3_bands[[i]], band[i])) 
}


# Checking tail dependence via chi measure
plot_chi_win(p3, "chi", upp_quantile = .9)


# Study of the GPD parameters
p3_gpd <- get_gpd(p3)
plot_gpd_shape(p3_gpd)
plot_gpd_scale(p3_gpd)


# We can do everything at once. Note: it may take a while (around 3.5 min).
start_time <- Sys.time()
p3_full <- sim_patient(patient_03)
end_time <- Sys.time()
end_time - start_time


# We can also plot everything at once. Results are saved into the plot sub-directory.
plot_all(p3_full)