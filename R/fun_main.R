#### Helper functions ####


#' Function: expand.grid.unique
#'
#' Create a data frame with only unique combinations of the supplied vector.
#' 
#' @param x A vector.
#' @param y A vector.
#' @param include.equals A logical indicating if only unique combinations are to be displayed (FALSE).
#' 
#' @return A data frame containing one row for each combination of the supplied vectors.
expand.grid.unique <- function(x, y, include.equals = FALSE){
  
  x <- unique(x)
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
  
}


#' Function: is_pos_semi_def
#' 
#' Verify if a matrix is positive semidefinite.
#' 
#' @param m A matrix (symmetric). Note that \code{m} will be forced to be symmetric.
#' @param m_diag A vector. The main diagonal of \code{m}.
#' @return A logical indicating whether \code{m} is positive semidefinite or not.
is_pos_semi_def <- function(m, m_diag) {
  
  m[lower.tri(m)] <- m[upper.tri(m)]
  diag(m) <- m_diag
  m <- as.matrix(Matrix:::forceSymmetric(m))
  
  return(matrixcalc:::is.positive.semi.definite(m))
  
}    


#### Main functions ####


#' Function: sim_ar2bands
#'
#' Simulate an AR(2) process for a given canonical frequency band.
#' 
#' @param sample_size An integer indicating the number of data points (in seconds) to be generated.
#' @param band A string indicating which canonical band to be simulated.
#' @param sampling_rate An integer indicating the number of data points per second.
#' 
#' @return A vector with the simulated AR(2) process.
sim_ar2bands <- function(
    sample_size, # in seconds
    sampling_rate = 256,
    band = c("delta", "theta", "alpha", "beta", "gamma")
) {
  
  ## IMPORTANT ##
  ## Change it as you wish. It controls the sharpness of the periodogram peaks.
  sharpness <- c(delta = .5, theta = .5, alpha = .5, beta = .5, gamma = .5)
  
  n <- sample_size * sampling_rate
  
  match.arg(band)
  lwHz <- c(0, 4, 8, 16, 32); upHz <- c(4, 8, 16, 32, 64)
  ## IMPORTANT ##
  ## Note that the periodogram peaks occur in the middle range of each frequency band.
  peak_location = lwHz + (upHz - lwHz) / 2
  names(peak_location) <- c("delta", "theta", "alpha", "beta", "gamma")
  
  
  psi <- 2 * peak_location[band] / sampling_rate
  par <- c(
    exp(-sharpness[band]) * 2 * cos(2 * pi * psi), 
    -exp(-2 * sharpness[band])
  )
  names(par) <- c("phi1", "phi2")
  
  signal <- arima.sim(
    n = n, 
    list(ar = par), 
    sd = 1, 
    n.start = 10000)
  
  return(as.vector(signal))
  
}


#' Function: get_latent_process
#'
#' Simulate a latent process \emph{Z} composed of the five canonical frequency bands generated using AR(2) processes.
#' 
#' @param sample_size An integer indicating the number of data points to be generated.
#' @param band A string indicating which canonical band to be simulated.
#' @param sampling_rate An integer indicating the number of samples per second.
#' 
#' @return A vector with the simulated latent process \emph{Z}
get_latent_process <- function(sample_size, sampling_rate) {
  
  delta <- sim_ar2bands(sample_size, sampling_rate, "delta")
  theta <- sim_ar2bands(sample_size, sampling_rate, "theta")
  alpha <- sim_ar2bands(sample_size, sampling_rate, "alpha")
  beta <- sim_ar2bands(sample_size, sampling_rate, "beta")
  gamma <- sim_ar2bands(sample_size, sampling_rate, "gamma")
  
  Z <- as.matrix(cbind(delta, theta, alpha, beta, gamma))
  colnames(Z) <- c("delta", "theta", "alpha", "beta", "gamma")
  
  return(Z)
  
}


#' Function: get_fixed_weights
#'
#' Organize a weight matrix/vector \emph{A}, which controls the influence of each frequency band in the latent process \emph{Z}.
#' 
#' @note These are the fixed weights of each frequency bands. For the weights varying in time check [get_timevarying_weights()].
#' 
#' @param weights A vector of length 5 indicating the influence of each frequency band in the case of NO seizure.
#' @param before A vector of length 5 indicating the influence of each frequency band BEFORE seizure.
#' @param during A vector of length 5 indicating the influence of each frequency band DURING seizure.
#' @param after A vector of length 5 indicating the influence of each frequency band AFTER seizure.
#' @param row_sumupto1 A logical indicating if the influence of the five frequency bands should sum up to 1 (row-wise sum).
#' 
#' @return A vector with the simulated latent process \emph{Z}.
#' 
#' @example 
#' W_1 = get_fixed_weights(
#'    before = c(0.00, 0.50, 3.00, 0.50, 2.00),
#'    during = c(1.00, 1.00, 8.00, 0.50, 12.0),
#'    after  = c(1.50, 0.00, 2.00, 0.00, 2.00),
#'    row_sumupto1 = TRUE
#' )
#' or
#' W_2 = get_fixed_weights(c(0.00, 0.50, 3.00, 0.50, 2.00))
get_fixed_weights <- function(
    weights = NULL, 
    before = NULL, 
    during = NULL, 
    after = NULL,
    row_sumupto1 = FALSE
) {
  
  if (is.null(weights)) {
    weights <- rbind(before, during, after)
    colnames(weights) <- c("delta", "theta", "alpha", "beta", "gamma")
    rownames(weights) <- c("before", "during", "after")
    if (row_sumupto1 == TRUE) {
      weights[1, ] <- weights[1, ] / sum(weights[1, ])
      weights[2, ] <- weights[2, ] / sum(weights[2, ])
      weights[3, ] <- weights[3, ] / sum(weights[3, ])
    }
  } else {
    names(weights) <- c("delta", "theta", "alpha", "beta", "gamma")
    if (row_sumupto1 == TRUE) {weights <- weights / sum(weights)}
  }
  
  return(weights)
  
} 


#' Function: get_timevarying_weights
#'
#' Replicate the weights of each frequency band (over the latent process \emph{Z}) accordingly with the time occurrence of a seizure. 
#' 
#' @note There are three distinct moments: before, during, and after the seizure. 
#' If there is no seizure, the weights are kept constant throughout the simulated time series.
#'
#' @param sample_size An integer indicating the number of data points to be generated.
#' @param sampling_rate An integer indicating the number of samples per second.
#' @param seizure_start An integer indicating the seizure onset.
#' @param seizure_end An integer indicating the end of the seizure process.
#' @param fixed_weights A matrix (or vector) of fixed weights for each frequency bands. The output of [get_fixed_weights()].
#' 
#' @return A matrix (dim = (n = sample_size x sampling_rate, 5)) of weights varying in time but kept constant within each seizure moment.
get_timevarying_weights <- function(
    sample_size, sampling_rate, 
    seizure_start = NULL, seizure_end = NULL, 
    fixed_weights 
) {
  
  n <- sample_size * sampling_rate
  
  if (is.null(seizure_start)) {
    
    ## NO seizure case.
    A <- as.matrix(
      cbind(
        rep(fixed_weights["delta"], n), 
        rep(fixed_weights["theta"], n), 
        rep(fixed_weights["alpha"], n), 
        rep(fixed_weights["beta"], n), 
        rep(fixed_weights["gamma"], n)
      )
    )
    colnames(A) <- c("delta", "theta", "alpha", "beta", "gamma")
    rownames(A) <- NULL
    
  } else {
    
    ## Seizure
    ## Establishing the moments before, during and after seizure.
    t_before <- ceiling(seizure_start * sampling_rate) - 1
    t_during <- ceiling((seizure_end - seizure_start) * sampling_rate) + 1
    t_after <- n - (t_before + t_during) 
    t_vec <- c(t_before, t_during, t_after)
    names(t_vec) <- c("before", "during", "after")
    
    ## Auxiliary to properly replicate the weights along the seizure moments.
    w_band <- function(w, t) {
      return(c(rep(w[1], t[1]), rep(w[2], t[2]), rep(w[3], t[3])))
    }
    
    ## Creating weight matrix
    A_delta <- w_band(fixed_weights[, "delta"], t_vec)
    A_theta <- w_band(fixed_weights[, "theta"], t_vec)
    A_alpha <- w_band(fixed_weights[, "alpha"], t_vec)
    A_beta <- w_band(fixed_weights[, "beta"], t_vec)
    A_gamma <- w_band(fixed_weights[, "gamma"], t_vec)
    
    A <- as.matrix(cbind(A_delta, A_theta, A_alpha, A_beta, A_gamma))
    colnames(A) <- c("delta", "theta", "alpha", "beta", "gamma")
    
  }
  
  return(A)
  
}  


#' Function: get_auxiliary_eeg
#'
#' Simulate an EEG signal composed of a latent process \emph{Z} weighted by a matrix \emph{A} with the addition of a heavy-tail component.
#'
#' @note You can pass any heavy-tail time series as the tail_factor. (Or, you could decide not to include it.)
#' If included, you can choose where the factor is inserted in the simulated signal: 
#' "inside" Z, sharing the same weight matrix A; 
#' or, "outside" Z, not been influenced by A.
#'
#' @param latent_z_process A matrix (nrow = n) with the time series of the five frequency bands composing the latent process \emph{Z};
#' output of [get_latent_process()].
#' @param weight_matrix A matrix (nrow = n) \emph{A} with the time varying weights given by the function [get_timevarying_weights()].
#' @param tail_factor A matrix (nrow = n) with the heavy-\code{tail factor} \emph{F} decomposed into the five frequency bands.
#' @param tail_factor_place A string indicating where the \code{tail_factor} \emph{F} should be placed.
#' 
#' @return A vector of the simulated EEG signal with length n.
get_auxiliary_eeg <- function(
    latent_z_process, 
    weight_matrix,
    tail_factor = NULL, 
    tail_factor_place = c(NULL, "in_z", "out_z")
) {
  
  match.arg(tail_factor_place)
  
  n <- nrow(latent_z_process)
  if (n != nrow(weight_matrix)) stop("Non-conformable elements.")
  
  X <- vector("numeric", length = n)
  if (is.null(tail_factor_place)) {
    for (i in 1:n) X[i] <- weight_matrix[i, ] %*% latent_z_process[i, ]
  } else if (tail_factor_place == "in_z") {
    for (i in 1:n) X[i] <- weight_matrix[i, ] %*% (latent_z_process[i, ] + tail_factor[i, ])
  } else {
    for (i in 1:n) X[i] <- weight_matrix[i, ] %*% latent_z_process[i, ] + sum(tail_factor[i, ])
  }
  
  return(X)
  
}


#' Function: sim_patient_eeg
#'
#' Simulate a seizure (or non-seizure) EEG signal with (or without) a heavy-tail factor. 
#' The signal is simulated based on a latent process composed by five AR(2) processes.
#'
#' @param pt_id An integer; unique patient identifier.
#' @param seizure A logical; 1 for seizure, 0 for non-seizure.
#' @param sample_size An integer indicating the number of data points (in seconds) to be generated.
#' @param sampling_rate An integer indicating the number of data points per second.
#' @param window_duration An integer; duration (in seconds) of the non-overlapping sliding window.
#' @param n_channels An integer; the quantity of EEG channels to be simulated.
#' @param n_clusters An integer indicating the number of distinct clusters to be built.
#' @param clusters A vector of length \code{n_clusters} indicating in which cluster a channel belongs. 
#' If \code{n_channels} = 5, and \code{n_clusters} = 3, an example of \code{clusters} is c(1, 3, 3, 1, 2). 
#' Channels 1 and 4 are in cluster 1. Channels 2 and 3 in cluster 3. And, channel 5 is in cluster 2.
#' @param fixed_weights A matrix or list of matrices depending on how many clusters do you have.
#' It is the weight matrix \emph{A}, controlling the influence of each frequency band in the latent process \emph{Z}.
#' @param seizure_start An integer indicating the seizure onset (in seconds).
#' @param seizure_end An integer indicating the ending of the seizure moment (in seconds).
#' @param error_covariance An upper triangular matrix (nrow = ncol = n_channels) giving the pairwise "bulk" correlation between the simulated channels.
#' @param error_variances A vector of length n_channels giving the variance of each simulated EEG channel.
#' @param tail_type A string indicating the distribution of the heavy-tail factor. Either NULL or "exponential".
#' @param tail_par A vector of length n_clusters with the parameter of the heavy-tail factor (cluster-wise). Set to NULL if no such factor is used.
#' @param tail_band_weight A vector of lenght 5 indicating the influence of each frequency band in the heavy-tail factor.
#' @param tail_factor_place A string indicating where to place the heavy-tail factor. 
#' Can be none (NULL), inside \emph{Z} ("in_z"), or outside \emph{Z} ("out_z").
#' 
#' @return A list with \emph{eeg}, a data frame of the simulated EEG channels, and \emph{setup} all function inputs.
sim_patient_eeg <- function(
    pt_id, 
    seizure, 
    sample_size,
    sampling_rate,
    window_duration,
    n_channels,
    n_clusters, 
    clusters, 
    fixed_weights, 
    seizure_start, 
    seizure_end, 
    error_covariance,
    error_variances,
    tail_type = c(NULL, "exponential"),
    tail_par = NULL,
    tail_band_weight = c(1, 1, 1, 1, 1),
    tail_factor_place = c(NULL, "in_z", "out_z")
) {
  
  
  match.arg(tail_type)
  match.arg(tail_factor_place)
  
  # set.seed(pt_id + 10122020)
  n <- sample_size * sampling_rate
  xtime <- seq(0, sample_size, length.out = n)
  
  ## Generating the latent process Z
  Z <- get_latent_process(sample_size, sampling_rate)
  
  ## Generating the weight matrix A
  A <- vector("list", length = n_clusters)
  for (i in 1:n_clusters) {
    A[[i]] <- get_timevarying_weights(
      sample_size,
      sampling_rate,
      seizure_start = seizure_start[i], 
      seizure_end = seizure_end[i], 
      fixed_weights = fixed_weights[[i]])
  }
  
  ## Normal Errors
  error_covariance[lower.tri(error_covariance)] <- error_covariance[upper.tri(error_covariance)]
  rownames(error_covariance) = colnames(error_covariance) <- paste0("Ch_", 1:n_channels)
  diag(error_covariance) <- error_variances
  error_covariance <- as.matrix(Matrix:::forceSymmetric(error_covariance))
  Error <- mvtnorm:::rmvnorm(n = n, mean = rep(0, n_channels), sigma = error_covariance)
  
  ## Tail factors
  w <- tail_band_weight
  Factor <- vector("list", length = n_clusters)
  if (is.null(tail_type)) {
    for(i in 1:n_clusters) Factor[[i]] <- NULL
  } else {
    for(i in 1:n_clusters) {
      Factor[[i]] <- matrix(
        c(w[1] * rexp(n, tail_par[i]),
          w[2] * rexp(n, tail_par[i]),
          w[3] * rexp(n, tail_par[i]),
          w[4] * rexp(n, tail_par[i]),
          w[5] * rexp(n, tail_par[i]))
        , nrow = n, ncol = 5)
      colnames(Factor[[i]]) <- c("delta", "theta", "alpha", "beta", "gamma")
    }
  } 
  
  
  # Simulated signals
  Chns <- vector("list", length = n_channels) 
  for (i in 1:n_channels) {
    Chns[[i]] <- as.vector(
      scale(
        get_auxiliary_eeg(
          Z, A[[clusters[i]]], Factor[[clusters[i]]], tail_factor_place
        ) + Error[, i]
        , center = TRUE, scale = TRUE
      )
    )
  }
  
  dt_chn <- as.data.frame(matrix(unlist(Chns), ncol = n_channels, nrow = n))
  names(dt_chn) <- paste0("Channel_", 1:n_channels)
                          
  dt_eeg <- cbind(data.frame(
    Patient = pt_id,
    Seizure = seizure,
    Time = xtime), 
    dt_chn
  )  
  
  setup <- list(
    pt_id = pt_id, 
    seizure = seizure, 
    sample_size = sample_size,
    sampling_rate = sampling_rate,
    window_duration = window_duration,
    n_channels = n_channels,
    n_clusters = n_clusters, 
    clusters = clusters, 
    fixed_weights = fixed_weights, 
    seizure_start = seizure_start, 
    seizure_end = seizure_end, 
    error_covariance = error_covariance,
    error_variances = error_variances,
    tail_type = tail_type,
    tail_par = tail_par,
    tail_band_weight = tail_band_weight,
    tail_factor_place = tail_factor_place
  )
  
  return(list(eeg = dt_eeg,
              setup = setup))
  
}  


#' Function: get_eeg_bands
#'
#' Decompose the simulated EEG signals into the five canonical frequency bands using a Butterworth filter.
#'
#' @param patient A list. The output of the function [sim_patient_eeg()]
#' 
#' @return A list of data frames. Each element is the simulated EEG signal in a canonical frequency band
get_eeg_bands <- function(patient) {
  
  
  sampling_rate <- patient$setup$sampling_rate
  
  signals <- patient$eeg[, -c(1:3)]
  
  my_filter <- function(x, bf) {
    as.vector(signal:::filter(bf, x))
  }  
  
  bands <- c("delta", "theta", "alpha", "beta", "gamma")
  lwHz <- c(0, 4, 8, 16, 32); upHz <- c(4, 8, 16, 32, 64)
  names(lwHz) = names(upHz) <- bands
  
  dt_bands <- vector("list", length = 5)
  names(dt_bands) <- bands
  for(i in seq_along(bands)) {
    startband <- lwHz[bands[i]]; endband <- upHz[bands[i]]
    bf <- signal:::butter(3, c(startband / sampling_rate, endband / sampling_rate))
    
    dt_aux <- apply(signals, 2, my_filter, bf = bf)
    dt_bands[[i]] <- cbind(patient$eeg[, 1:3], dt_aux)
  }
  
  return(dt_bands)
  
}


#### EVT based functions ####


#' Function: get_evi
#'
#' Estimate the GPD extreme value index (evi) of seizure and non-seizure moments from a simulated EEG signal using [get_auxiliary_eeg()] 
#'
#' @param signal A vector. An EEG signal, the output of [get_auxiliary_eeg()].  
#' @param sampling_rate An integer indicating the number of data points per second.
#' @param seizure_start An integer indicating the seizure onset (in seconds).
#' @param seizure_end An integer indicating the ending of the seizure moment (in seconds).
#' 
#' @return A data frame with the estimated GPD evi for seizure and non-seizure moments.
get_evi <- function(signal, sampling_rate, seizure_start, seizure_end) {
  
  k_seq <- c(50, 75, 100, 125, 150)
  
  get_par <- function(signal) {
    evi <- vector("numeric", length = length(k_seq))
    for (k in seq_along(k_seq)) {
      evi[k] <- eva:::gpdFit(signal
                             , nextremes = k_seq[k]
                             , method = c("mle")
                             , start = NULL
                             , opt = "Nelder-Mead"
                             , maxit = 10000
      )$par.ests[2]
    } 
    return(round(evi, 4))
  }
  
  seizure <- (seizure_start*sampling_rate):(seizure_end*sampling_rate)
  
  dt <- data.frame(
    Exceedances = k_seq,
    Seizure = get_par(signal[seizure]),
    `Non_seizure` = get_par(signal[-seizure])
  )
  
  return(dt)
  
} 



#' Function: get_gpd
#'
#' Estimate the GPD extreme value index (evi) and scale from the EEG signals simulated using the function [sim_patient_eeg()].
#'
#' @param patient A list. The output of the function [sim_patient_eeg()]
#' @param patient_band A list. The output of the function [get_eeg_bands()] for the given \code{band}.
#' @param band A string. The canonical frequency band for which you want to estimate the GPD parameters.
#' 
#' @return A data frame with the estimated GPD evi and scale the EEG signal.
get_gpd <- function(
    patient, 
    patient_band = NULL, 
    band = c(NULL, "delta", "theta", "alpha", "beta", "gamma")
    
) {
  
  match.arg(band)
  
  window_duration <- patient$setup$window_duration
  sampling_rate <- patient$setup$sampling_rate
  sample_size <- patient$setup$sample_size
  n <- sampling_rate * sample_size
  
  event <- ifelse(
    patient$setup$seizure == 1, 
    "Seizure",
    "Non-seizure")
  
  pt_id <- patient$setup$pt_id
  
  if (is.null(patient_band)) {
    dt_split <- patient$eeg[, -c(1:3)]
    input <- "Simulated signal"
  } else {
    dt_split <- patient_band[, -c(1:3)]
    input <- paste0("Filtered ", band, " band")
  }
  
  window_n <- sample_size / window_duration
  window_size <- window_duration * sampling_rate
  blk <- gl(n = window_n, k = window_size, length = n)
  slices <- split(dt_split, blk)
  
  n_chns <- patient$setup$n_channels
  k_seq <- seq(40, 120, by = 1)
  dt_gpd <- vector("list", length = window_n)
  for (i in 1:window_n) {
    dt_chn <- vector("list", length = n_chns)
    for (j in 1:n_chns) {
      dt_par <- vector("list", length = length(k_seq))
      for (k in seq_along(k_seq)) {
        gpd_par <-  tryCatch(
          eva:::gpdFit(slices[[i]][, j]
                       , nextremes = k_seq[k]
                       , method = c("mle")
                       , start = NULL
                       , opt = "Nelder-Mead"
                       , maxit = 10000
          )$par.ests, 
          error = function(e) {gpd_par <- NA}
        )
        dt_par[[k]] <- data.frame(
          Patient = pt_id,
          Seizure = event,
          Channel = paste0("Chn-", j),
          Scale   = as.double(gpd_par[1]),
          Shape   = as.double(gpd_par[2]),
          k = k_seq[k],
          Window = i,
          Band = input
        )
      }
      dt_chn[[j]] <- do.call("rbind", dt_par)
    }
    dt_gpd[[i]] <- do.call("rbind", dt_chn)
  }
  dt_gpd <- do.call("rbind", dt_gpd)
  
  alpha <- 0.05
  sig <- alpha / (2 * n_chns)
  z <- qnorm(1 - sig)
  dt_gpd$Shp_lwb <- dt_gpd$Shape - z / sqrt(dt_gpd$k) * (1 + dt_gpd$Shape)
  dt_gpd$Shp_upb <- dt_gpd$Shape + z / sqrt(dt_gpd$k) * (1 + dt_gpd$Shape)
  dt_gpd$Scl_lwb <- dt_gpd$Scale - z / sqrt(dt_gpd$k) * (1 + dt_gpd$Scale)
  dt_gpd$Scl_upb <- dt_gpd$Scale + z / sqrt(dt_gpd$k) * (1 + dt_gpd$Scale)
  dt_gpd$Scl_lwb[dt_gpd$Scl_lwb < 0] <- 0
  
  return(dt_gpd)
  
}    


#' Function: get_gpd_all
#'
#' Estimate the GPD extreme value index (evi) and scale from the EEG signals simulated using the function [sim_patient_eeg()].
#'
#' @param patient A list. The output of the function [sim_patient_eeg()].
#' @param patient_band A list. The output of the function [get_eeg_bands()] for all signals and theirs spectral decomposition.
#' 
#' @return A data frame with the estimated GPD evi and scale for the EEG signals and theirs spectral decomposition.
get_gpd_all <- function(
    patient,
    patient_band
) {
  
  n_chns <- patient$setup$n_channels
  dt_gpd <- vector("list", n_chns)
  band = c("delta", "theta", "alpha", "beta", "gamma")
  dt_gpd[[1]] <- get_gpd(patient)
  
  for (i in seq_along(band)) dt_gpd[[i + 1]] <- get_gpd(patient, patient_band[[i]], band[i])
  
  dt_gpd <- do.call("rbind", dt_gpd)
  
}  


#### Simulating a complete scenario ####


sim_patient <- function(
    patient_scenario, 
    output_folder = "./output/"
) {
  
  p <- patient_scenario
  
  output_folder = ifelse(
    is.null(p$tail_factor_place), 
    output_folder,
    paste0(output_folder, p$tail_factor_place, "/")
  )
  
  p_signal <- sim_patient_eeg(
    p$pt_id, 
    p$seizure, 
    p$sample_size, 
    p$sampling_rate, 
    p$window_duration,
    p$n_channels, 
    p$n_clusters, 
    p$clusters, 
    p$fixed_weights, 
    p$seizure_start, 
    p$seizure_end, 
    p$error_covariance, 
    p$error_variances, 
    p$tail_type, 
    p$tail_par, 
    p$tail_band_weight, 
    p$tail_factor_place
  )
  
  p_bands <- get_eeg_bands(p_signal)
  
  p_gpd <- get_gpd_all(p_signal, p_bands)
  
  patient <- list(signal = p_signal, bands = p_bands, gpd = p_gpd)
  write_rds(patient, 
            paste0(output_folder, "patient-", p$pt_id, "-", p$seizure, ".RDS")
  )
  
  return(patient)
  
}