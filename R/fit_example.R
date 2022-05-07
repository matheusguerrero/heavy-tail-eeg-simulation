rm(list = ls())

source("./R/set_packages.R")
source("./R/fun_main.R")
source("./R/fun_plot.R")


#### Common parameters ####
sample_size = 100
sampling_rate = 256
window_duration = 1
n <- sample_size * sampling_rate


# Latent process: common for all patients. (Is this realistic?)
Z = get_latent_process(sample_size, sampling_rate)
dim(Z)
n == nrow(Z)
head(Z, n = 2)


# Plotting the latent process and its spectral decomposition.
signal <- rowSums(Z)
bands = c("signal", "delta", "theta", "alpha", "beta", "gamma")
for (i in seq_along(bands)) plot_bands(1:length(signal), signal, bands[i], sampling_rate)


# White noise to be used for all patients a and b.
WNa <- rnorm(n); WNb <- rnorm(n)


#### Non seizure patient ####


# Non-seizure weights.

# Patient a
(Wa_1 = get_fixed_weights(c(0.50, 6.00, 3.00, 1.00, 8.00), row_sumupto1 = TRUE))
sum(Wa_1)
Aa_1 = get_timevarying_weights(sample_size, sampling_rate, fixed_weights = Wa_1)
dim(Aa_1)
head(Aa_1, n = 2)


# Patient b
(Wb_1 = get_fixed_weights(c(9.00, 0.05, 2.00, 8.00, 3.00), row_sumupto1 = TRUE))
sum(Wb_1)
Ab_1 = get_timevarying_weights(sample_size, sampling_rate, fixed_weights = Wb_1)
dim(Ab_1)
head(Ab_1, n = 2)


# Non-seizure EEG without heavy-tail factor.


# Patient a
eeg_a_1 <- get_auxiliary_eeg(Z, Aa_1, tail_factor = NULL, tail_factor_place = NULL)
signal_a_1 <- as.vector(scale(eeg_a_1 + WNa, T, F))
plot_patient(signal_a_1, sample_size, sampling_rate, "Non-seizure EEG - Patient a \n Without heavy-tail factor")


# Patient b
eeg_b_1 <- get_auxiliary_eeg(Z, Ab_1, tail_factor = NULL, tail_factor_place = NULL)
signal_b_1 <- as.vector(scale(eeg_b_1 + WNb, T, F))
plot_patient(signal_b_1, sample_size, sampling_rate, "Non-seizure EEG - Patient b \n Without heavy-tail factor")


# Checking extremal dependence
plot(texmex::chi(cbind(signal_a_1, signal_b_1)))


#### Seizure patients ####


# Time occurrence (begin and ending) of seizure in seconds.
seizure_start = 40; seizure_end = 88


# Seizure weights.

# Patient a
(Wa_2 = get_fixed_weights(
before = c(0.00, 0.50, 3.00, 0.50, 5.00),
during = c(1.00, 1.00, 8.00, 0.50, 12.0),
after  = c(1.50, 0.00, 2.00, 0.00, 6.00),
row_sumupto1 = FALSE # IMPORTANT: If TRUE, we can't have the seizure behavior.
))
rowSums(Wa_2)
Aa_2 = get_timevarying_weights(sample_size, sampling_rate, seizure_start, seizure_end, fixed_weights = Wa_2)
dim(Aa_2)
rbind(head(Aa_2, n = 2), tail(Aa_2, n = 2))


# Seizure EEG without heavy-tail factor.
eeg_a_2 <- get_auxiliary_eeg(Z, Aa_2, tail_factor = NULL, tail_factor_place = NULL)
signal_a_2 <- as.vector(scale(eeg_a_2 + WNa, T, F))
plot_patient(signal_a_2, sample_size, sampling_rate, "Seizure EEG - Patient a \n Without heavy-tail factor")


# Checking GPD's shape parameter for seizure and non-seizure moments.
get_evi(eeg_a_2, sampling_rate, seizure_start, seizure_end)


# Exponential heavy-tail factor.
exp_rate <- 1.8
delta_exp <- rexp(n, exp_rate)
theta_exp <- rexp(n, exp_rate)
alpha_exp <- rexp(n, exp_rate)
beta_exp <- rexp(n, exp_rate)
gamma_exp <- rexp(n, exp_rate)
Factor <- matrix(c(delta_exp, theta_exp, alpha_exp, beta_exp, gamma_exp), nrow = n, ncol = 5)
colnames(Factor) <- c("delta", "theta", "alpha", "beta", "gamma")
dim(Factor)
head(Factor, n = 2)


# Seizure EEG with heavy-tail factor outside the latent process Z.
eeg_a_3 <- get_auxiliary_eeg(Z, Aa_2, Factor, tail_factor_place = "out_z")
signal_a_3 <- as.vector(scale(eeg_a_3 + WNa, T, F))
plot_patient(signal_a_3, sample_size, sampling_rate, "Seizure EEG - Patient a \n With heavy-tail factor outside Z")


# Checking GPD's shape parameter for seizure and non-seizure moments.
get_evi(eeg_a_3, sampling_rate, seizure_start, seizure_end)

# Seizure EEG with weighted heavy-tail factor inside the latent process Z.
tail_bands_weights <- c(1.0, 0.5, 5.0, 0.5, 8.0); w <- tail_bands_weights
weighted_Factor <- matrix(c(w[1] * delta_exp, w[2] * theta_exp, w[3] * alpha_exp, w[4] * beta_exp, w[5] * gamma_exp), nrow = n, ncol = 5)
colnames(weighted_Factor) <- colnames(Factor)
dim(weighted_Factor)
head(weighted_Factor, n = 2)

eeg_a_4 <- get_auxiliary_eeg(Z, Aa_2, weighted_Factor, tail_factor_place = "out_z")
signal_a_4 <- as.vector(scale(eeg_a_4 + WNa, T, F))
plot_patient(signal_a_4, sample_size, sampling_rate, "Seizure EEG - Patient a \n With weighted heavy-tail factor outside Z")
### IMPORTANT: Maybe this is the best way. At least the simulated signals are more natural.


# Checking GPD's shape parameter for seizure and non-seizure moments.
get_evi(eeg_a_4, sampling_rate, seizure_start, seizure_end)

# Seizure EEG with heavy-tail factor inside the latent process Z.
eeg_a_5 <- get_auxiliary_eeg(Z, Aa_2, Factor, tail_factor_place = "in_z")
signal_a_5 <- as.vector(scale(eeg_a_5 + WNa, T, F))
plot_patient(signal_a_5, sample_size, sampling_rate, "Seizure EEG - Patient a \n With heavy-tail factor inside Z")


# Checking GPD's shape parameter for seizure and non-seizure moments.
get_evi(eeg_a_5, sampling_rate, seizure_start, seizure_end)


# Seizure EEG with weighted heavy-tail factor inside the latent process Z.
eeg_a_6 <- get_auxiliary_eeg(Z, Aa_2, weighted_Factor, tail_factor_place = "in_z")
signal_a_6 <- as.vector(scale(eeg_a_6 + WNa, T, F))
plot_patient(signal_a_6, sample_size, sampling_rate, "Seizure EEG - Patient a \n With weighted heavy-tail factor inside Z")


# Checking GPD's shape parameter for seizure and non-seizure moments.
get_evi(eeg_a_6, sampling_rate, seizure_start, seizure_end)


# Plotting the different signals together.
x <- (1:n)/sampling_rate
par(mfrow = c(6, 1))
plot_single_signal(x, signal_a_1, main = "Non-seizure EEG \n Without heavy-tail factor")
plot_single_signal(x, signal_a_2, main = "Seizure EEG \n Without heavy-tail factor")
plot_single_signal(x, signal_a_3, main = "Seizure EEG \n With heavy-tail factor outside Z")
plot_single_signal(x, signal_a_4, main = "Seizure EEG \n With weighted heavy-tail factor outside Z")
plot_single_signal(x, signal_a_5, main = "Seizure EEG \n With heavy-tail factor inside Z")
plot_single_signal(x, signal_a_6, main = "Seizure EEG \n With weighted heavy-tail factor inside Z")
par(mfrow = c(1, 1))