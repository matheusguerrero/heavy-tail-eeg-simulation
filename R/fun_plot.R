# Customized plot theme
my_theme <- theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold")
  )

# Helper fuction for better ggplot2 axis ticks
number_ticks <- function(n) {function(limits) pretty(limits, n)}


#' Function: get_periodogram
#'
#' Create a data frame with the estimated periodogram of a given signal.
#' 
#' @param signal A vector with a univariate time series.
#' @param id A string indicating a label for \code{signal}.
#' 
#' @return A data frame containing the estimated periodogram of \code{signal}.
get_periodogram <- function(signal, id) {
  
  perdgram <- TSA::periodogram(signal, plot = FALSE)
  dt_perdgram <- data.frame(Channel = id, Freq = perdgram$freq, Spec = perdgram$spec)
  return(dt_perdgram)
  
}


plot_bands <- function(
    x, y, 
    band = c("signal", "delta", "theta", "alpha", "beta", "gamma"),
    sampling_rate
) {
  
  match.arg(band)
  
  if (band != "signal") {
    lwHz <- c(0, 4, 8, 16, 32); upHz <- c(4, 8, 16, 32, 64)
    names(lwHz) =  names(upHz)  <- c("delta", "theta", "alpha", "beta", "gamma")
    startband <- lwHz[band]; endband <- upHz[band]
    
    bf <- signal:::butter(3, c(4 * startband / sampling_rate, 4 * endband / sampling_rate))
    y <- signal:::filter(bf, y)
  }
  
  Band <- c(signal = "Simulated signal", delta = "Delta band", theta = "Theta band", 
            alpha  = "Alpha band",  beta  = "Beta band",  gamma = "Gamma band")
  my_title <- paste0(Band[band], " | Freq: ", sampling_rate, " Hz \n 1 second sample")
  
  par(mfrow = c(3, 1))
  plot(x = x[1:sampling_rate], y = y[1:sampling_rate], type = "l", main = my_title, xlab = "Time", ylab = " ")
  plot(x = x, y = y, type = "l", main = "10 seconds sample", xlab = "Time [s]", ylab = " ")
  periodogram(y, main = "Periodogram", ylab = " ")
  par(mfrow = c(1, 1))
  
}  

plot_single_signal <- function(x, signal, ...) {
  plot(x, signal, type = "l", 
    xlab = "Time [s]", ylab = "EEG Amplitude", xaxs="i", xaxt = "n",...
  )
  axis(1, at = c(20, 40, 60, 80))
}       

plot_spectrogram <- function(signal, nfft, sampling_rate, window, overlap,...) {
  
  # create spectrogram
  spec = signal::specgram(x = signal, n = nfft, Fs = sampling_rate, window = window, overlap = overlap)
  
  # discard phase information
  P = abs(spec$S)
  #normalize
  P = P/max(P)
  # config time axis
  t = spec$t
  
  # plot spectrogram
  oce::imagep(x = t,
              y = spec$f,
              z = t(P),
              ylim = c(0, 64),
              col = oce.colorsViridis,
              ylab = 'Frequency [Hz]',
              xlab = 'Time [s]',
              drawPalette = T,
              decimate = F,
              axes = F,
              ...)
  axis(1, at = c(20, 40, 60, 80))
  axis(2, at = c(4, 8, 16, 32, 64), labels = c(4, 8, 16, 32, 64))
  box()
}


plot_patient <- function(signal, sample_size, sampling_rate, title) {
  
  n <- sample_size * sampling_rate
  
  par(mfrow=c(2,1), oma = c(0, 0, 3, 0))
  plot_single_signal(x = (1:n)/sampling_rate, signal)
  plot_spectrogram(signal, nfft = 1024, sampling_rate, window = 256, overlap = 128)
  mtext(title, outer = TRUE, cex = 1.5, font = 2, adj = 0)
  par(mfrow = c(1, 1))
  
}


plot_signal <- function(
    patient, 
    dt = NULL, 
    band = c(NULL, "delta", "theta", "alpha", "beta", "gamma")
) {
  
  match.arg(band)
  
  pt_id <- patient$setup$pt_id
  
  event <- ifelse(
    patient$setup$seizure == 1, 
    "Seizure",
    "Non-seizure")
  
  if (is.null(dt)) {
    dt_aux <- patient$eeg[, -c(1:3)]
    dt <- patient$eeg
    input <- "Simulated signal"
  } else {
    dt_aux <- dt[, -c(1:3)]
    input <- paste0("Filtered ", band, " band")
  }
  
  n_chns <- patient$setup$n_channels
  dt_perdgram <- vector("list", length = n_chns)
  for(i in 1:n_chns) {
    dt_perdgram[[i]] <- get_periodogram(dt_aux[, 1], paste0("Channel_", i)) 
  }
  dt_perdgram <- do.call("rbind", dt_perdgram)
  dt_perdgram$Patient <- pt_id
  dt_perdgram$Seizure <- patient$setup$seizure
  
  
  dt %>% 
    pivot_longer(-c(Patient, Time, Seizure), 
                 names_to = "Channel", 
                 values_to = "Amplitude") %>% 
    ggplot(aes(x = Time, y = Amplitude)) +
    geom_line() +
    scale_x_continuous(breaks = number_ticks(5)) +
    scale_y_continuous(breaks = number_ticks(3)) +
    labs(x = "Time [s]", y = "EEG Amplitude") +
    facet_wrap(~Channel, nrow = 1) +
    my_theme -> plot_signals
  
  dt_perdgram %>% 
    ggplot(aes(x = Freq, y = Spec)) +
    geom_segment(aes(x = Freq, xend = Freq, y = 0, yend = Spec), size = .25) +
    scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4)) +
    scale_y_continuous(breaks = number_ticks(3)) +
    labs(x = "Frequency [Hz]", y = "Spectrum") +
    facet_wrap(~Channel, nrow = 1) +
    my_theme +
    theme(strip.text.x = element_blank()) -> plot_perdgram
  
  my_plots <- plot_grid(plot_signals, plot_perdgram, ncol = 1, nrow = 2)
  title <- ggdraw() + 
    draw_label("Simulated EEG Signals w/ Periodogram | Freq: 256 Hz", 
               fontface='bold', x = 0, hjust = 0, size = 18) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  subtitle <- ggdraw() + 
    draw_label(paste0("Patient ", pt_id, " | ", event, " | ", input), 
               fontface='bold', x = 0, hjust = 0, size = 16) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  plot_signals <- plot_grid(title, subtitle, my_plots, ncol = 1, rel_heights = c(0.05, 0.05, 1))
  
  return(plot_signals)
  
} 


plot_chi <- function(
    patient, 
    type = c("chi", "chibar"),
    dt = NULL, 
    band = c(NULL, "delta", "theta", "alpha", "beta", "gamma")
) {
  
  
  match.arg(type)
  match.arg(band)
  
  if (is.null(dt)) {
    dt_aux <- patient$eeg[, -c(1:3)]
    input <- "Simulated signal"
  } else {
    dt_aux <- dt[, -c(1:3)]
    input <- paste0("Filtered ", band, " band")
  }
  
  pt_id <- patient$setup$pt_id
  seizure <- patient$setup$seizure
  
  
  pairs <- expand.grid.unique(1:6, 1:6)
  dt_chi <- vector("list", length = nrow(pairs))
  dt_chibar <- vector("list", length = nrow(pairs))
  for(i in 1:nrow(pairs)) {
    
    chi_i <- texmex:::chi(
      cbind(
        abs(dt_aux[, pairs[i, 1]]), 
        abs(dt_aux[, pairs[i, 2]])
      )
    )
    
    dt_chi[[i]] <- data.frame(chi_i$chi, quantile = chi_i$quantile)
    dt_chibar[[i]] <- data.frame(chi_i$chibar, quantile = chi_i$quantile)
    dt_chi[[i]]$Pair = dt_chibar[[i]]$Pair <- paste0("Ch", pairs[i, 1], " x Ch", pairs[i, 2])
    
  }
  dt_chi <- do.call("rbind", dt_chi)
  dt_chibar <- do.call("rbind", dt_chibar)
  dt_chi$Patient = dt_chibar$Patient <- pt_id
  dt_chi$Seizure = dt_chibar$Seizure <- seizure
  names(dt_chibar)[1:3] = names(dt_chi)[1:3]
  
  if (type == "chi") {
    dt_plot <- dt_chi
    mytitle <- "Tail Correlation Coefficient"
    y_label <- expression(chi(u))
  } else {
    dt_plot <- dt_chibar
    mytitle <- "Complementary Tail Correlation Coefficient"
    y_label <- expression(bar(chi)(u))
  }
  
  
  event <- ifelse(seizure == 1, "Seizure", "Non-seizure")
  
  
  dt_plot %>% 
    ggplot(aes(x = quantile, y = chi)) +
    geom_line(size = .80) +
    geom_ribbon(aes(ymin = chilow, ymax = chiupp), 
                alpha = 0.1, linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = c(-1, 0, 1), 
               linetype = "dashed", size = .65, color = "darkgrey") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0.25, 0.5, 0.75), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.25, 1.25), breaks = c(-1, -.5, 0, .5, 1), expand = c(0, 0)) +
    facet_wrap(~Pair) +
    ggtitle(mytitle,
            subtitle = paste0("Patient ", pt_id, " | ", event, " | ", input)) +
    labs(x = "Quantile", y = y_label) +
    my_theme + 
    theme(panel.grid.minor = element_blank()) -> chi_plot
  
  return(chi_plot)
  
}   


plot_chi_win <- function(
    patient, 
    type = c("chi", "chibar"),
    window_duration, 
    upp_quantile = .8,
    dt = NULL, 
    band = c(NULL, "delta", "theta", "alpha", "beta", "gamma")
) {
  
  
  match.arg(type)
  match.arg(band)
  
  if (is.null(dt)) {
    dt_split <- patient$eeg[, -c(1:3)]
    input <- "Simulated signal"
  } else {
    dt_split <- dt[, -c(1:3)]
    input <- paste0("Filtered ", band, " band")
  }
  
  sampling_rate <- 256
  n <- sampling_rate * patient$setup$sample_size
  
  
  window_n <- patient$setup$sample_size / window_duration
  window_size <- window_duration * sampling_rate
  blk <- gl(n = window_n, k = window_size, length = n)
  slices <- split(dt_split, blk)
  
  
  pairs <- expand.grid.unique(1:6, 1:6)
  dt_chi <- vector("list", length = length(slices))
  dt_chibar <- vector("list", length = length(slices))
  for (j in 1:length(slices)) {
    chi_slice <- vector("list", length = nrow(pairs))
    chibar_slice <- vector("list", length = nrow(pairs))
    for(i in 1:nrow(pairs)) {
      chi_i <- texmex:::chi(
        cbind(
          abs(slices[[j]][, pairs[i, 1]]), 
          abs(slices[[j]][, pairs[i, 2]])
        ), nq = 1, qlim = c(upp_quantile, upp_quantile)
      )
      chi_slice[[i]] <- data.frame(
        Chilow = chi_i$chi[1], 
        Chi = chi_i$chi[2], 
        Chiupp = chi_i$chi[3], 
        Quantile = upp_quantile,
        Pair = paste0("Ch", pairs[i, 1], " x Ch", pairs[i, 2]),
        Patient = patient$setup$pt_id,
        Seizure = patient$setup$seizure,
        Min = (j - 1), 
        Max = j
      )
      chibar_slice[[i]] <- data.frame(
        Chilow = chi_i$chibar[1], 
        Chi = chi_i$chibar[2], 
        Chiupp = chi_i$chibar[3], 
        Quantile = upp_quantile,
        Pair = paste0("Ch", pairs[i, 1], " x Ch", pairs[i, 2]),
        Patient = patient$setup$pt_id,
        Seizure = patient$setup$seizure,
        Min = (j - 1), 
        Max = j
      )
      rownames(chi_slice[[i]]) = rownames(chibar_slice[[i]]) <- NULL
    }
    dt_chi[[j]] <- do.call("rbind", chi_slice)
    dt_chibar[[j]] <- do.call("rbind", chibar_slice)
  }
  dt_chi <- do.call("rbind", dt_chi)
  dt_chibar <- do.call("rbind", dt_chibar)
  
  
  if (type == "chi") {
    dt <- dt_chi
    mytitle <- "Tail Correlation Coefficient for the different slide windows"
    y_label <- bquote(chi(.(upp_quantile)))
  } else {
    dt <- dt_chibar
    mytitle <- "Complementary Tail Correlation Coefficient for the different slide windows"
    y_label <- bquote(bar(chi)(.(upp_quantile)))
  }
  
  
  event <- ifelse(
    patient$setup$seizure == 1, 
    "Seizure",
    "Non-seizure")
  
  pt_id <- patient$setup$pt_id
  
  
  cols <- c("chi" = "black", "ci" = "black")
  cols_labs <- c("chi" = "Chi Measure", "ci" = "95% C.I.")
  lines_labs <- c("chi" = "solid", "ci" = "dashed")
  
  dt_tail <- dt %>% 
    group_by(Pair) %>% 
    arrange(Max) %>%
    dplyr::filter(row_number() == n()) %>% 
    mutate(Min = Min + 1,
           Max = Max + 1)
  dt <- rbind(dt, dt_tail)
  
  
  edges <- unique(dt$Min)
  scaleFUN <- function(x) sprintf("%.0f", x)
  dt %>% 
    ggplot() +
    geom_step(aes(x = Min, y = Chi, color = "chi", linetype = "chi"), 
              size = .8) +
    geom_step(aes(x = Min, y = Chilow, color = "ci", linetype = "ci"), size = .8) +
    geom_step(aes(x = Min, y = Chiupp, color = "ci", linetype = "ci"), size = .8) +
    geom_hline(yintercept = c(-1, 0, 1), 
               linetype = "dashed", size = .65, color = "darkgrey") +
    scale_x_continuous(breaks = edges, labels = scaleFUN) +
    scale_y_continuous(breaks = c(-1, -.5, 0, .5, 1),
                       labels = c(-1, "", 0, "", 1)) +
    coord_cartesian(ylim = c(-1, 1)) +
    facet_wrap(~Pair) +
    scale_color_manual(name = "my_legend", values = cols, labels = cols_labs) +
    scale_linetype_manual(name = "my_legend", values = lines_labs, labels = cols_labs) +
    ggtitle(mytitle,
            subtitle = paste0("Patient ", pt_id, " | ", event, " | ", input)) +
    labs(x = "Time [s]", y = y_label) +
    my_theme +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(1.2,"cm")
    ) -> plot_chi
  
  return(plot_chi)
  
}


plot_gpd_shape <- function(dt_gpd) {
  
  pt_id <- unique(dt_gpd$Patient)
  event <- unique(dt_gpd$Seizure)
  input <- unique(dt_gpd$Band)
  
  dt_gpd %>% 
    ggplot(aes(x = k, y = Shape)) +
    geom_line() +
    geom_ribbon(aes(ymin = Shp_lwb, 
                    ymax = Shp_upb), 
                alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", size = .65, color = "darkgrey") +
    scale_x_continuous(breaks = c(40, 60, 80, 100, 120),
                       labels = c("", 60, "", 100, "")) +
    scale_y_continuous(breaks = number_ticks(3)) +
    facet_grid(Channel~Window) +
    ggtitle("GPD Extreme value index estimates for various thresholds for the different slide windows",
            subtitle = paste0("Patient ", pt_id, " | ", event, " | ", input)) +
    labs(x = "Number of exceedances", 
         y = expression(hat(gamma))) +
    my_theme -> shape
  
  return(shape)
  
}  


plot_gpd_scale <- function(dt_gpd) {
  
  pt_id <- unique(dt_gpd$Patient)
  event <- unique(dt_gpd$Seizure)
  input <- unique(dt_gpd$Band)
  
  dt_gpd %>% 
    ggplot(aes(x = k, y = Scale)) +
    geom_line() +
    geom_ribbon(aes(ymin = Scl_lwb, 
                    ymax = Scl_upb), 
                alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", size = .65, color = "darkgrey") +
    scale_x_continuous(breaks = c(40, 60, 80, 100, 120),
                       labels = c("", 60, "", 100, "")) +
    scale_y_continuous(breaks = number_ticks(3)) +
    facet_grid(Channel~Window) +
    ggtitle("GPD scale estimates for various thresholds for the different slide windows",
            subtitle = paste0("Patient ", pt_id, " | ", event, " | ", input)) +
    labs(x = "Number of exceedances", 
         y = expression(hat(alpha))) +
    my_theme -> scale
  
  return(scale)
  
} 


plot_signal_all <- function(
    patient, 
    patient_bands,
    output_folder) {
  
  pt <- patient$setup$pt_id
  sz <- patient$setup$seizure
  
  dir <- paste0(output_folder, "eeg-", pt, "-", sz, "-0.PNG")
  png(dir, width = 1080, height = 720)
  plot_signal(patient) -> plot
  print(plot)
  dev.off()
  
  bands <- c("delta", "theta", "alpha", "beta", "gamma")
  for (i in 1:5) {
    dir <- paste0(output_folder, "eeg-", pt, "-", sz, "-", i, ".PNG")
    png(dir, width = 1080, height = 720)
    plot_signal(patient, patient_bands[[i]], bands[i]) -> plot
    print(plot)
    dev.off()
  }  
  
}


plot_chi_all <- function(
    patient, 
    patient_bands,
    window_duration,
    upp_quantile = .85,
    output_folder) {
  
  pt <- patient$setup$pt_id
  sz <- patient$setup$seizure
  
  # Chi plot - Simulated signal
  dir <- paste0(output_folder, "chi-", pt, "-", sz, "-0.PNG")
  png(dir, width = 1080, height = 720)
  plot_chi(patient, "chi") -> plot
  print(plot)
  dev.off()
  
  # Chi plot - Frequency bands
  bands <- c("delta", "theta", "alpha", "beta", "gamma")
  for(i in 1:5) {
    dir <- paste0(output_folder, "chi-", pt, "-", sz, "-", i, ".PNG")
    png(dir, width = 1080, height = 720)
    plot_chi(patient, "chi", patient_bands[[i]], bands[i]) -> plot
    print(plot)
    dev.off()
  }
  
  # Chibar plot - Simulated signal
  dir <- paste0(output_folder, "chibar-", pt, "-", sz, "-0.PNG")
  png(dir, width = 1080, height = 720)
  plot_chi(patient, "chibar") -> plot
  print(plot)
  dev.off()
  
  # Chibar plot - Frequency bands
  for(i in 1:5) {
    dir <- paste0(output_folder, "chibar-", pt, "-", sz, "-", i, ".PNG")
    png(dir, width = 1080, height = 720)
    plot_chi(patient, "chibar", patient_bands[[i]], bands[i]) -> plot
    print(plot)
    dev.off()
  }
  
  # Chi plot window-wise - Simulated signal
  dir <- paste0(output_folder, "chi-w-", pt, "-", sz, "-0.PNG")
  png(dir, width = 1080, height = 720)
  plot_chi_win(patient, "chi", window_duration, upp_quantile) -> plot
  print(plot)
  dev.off()
  
  # Chi plot window-wise - Frequency bands
  for(i in 1:5) {
    dir <- paste0(output_folder, "chi-w-", pt, "-", sz, "-", i, ".PNG")
    png(dir, width = 1080, height = 720)
    plot_chi_win(patient, "chi", window_duration, upp_quantile, patient_bands[[i]], bands[i]) -> plot
    print(plot)
    dev.off()
  }
  
  # Chibar plot window-wise - Simulated signal
  dir <- paste0(output_folder, "chibar-w-", pt, "-", sz, "-0.PNG")
  png(dir, width = 1080, height = 720)
  plot_chi_win(patient, "chibar", window_duration, upp_quantile) -> plot
  print(plot)
  dev.off()
  
  # Chi plot window-wise - Frequency bands
  for(i in 1:5) {
    dir <- paste0(output_folder, "chibar-w-", pt, "-", sz, "-", i, ".PNG")
    png(dir, width = 1080, height = 720)
    plot_chi_win(patient, "chibar", window_duration, upp_quantile, patient_bands[[i]], bands[i]) -> plot
    print(plot)
    dev.off()
  }  
  
}


plot_gpd_all <- function(
    patient, 
    patient_bands,
    gpd_bands,
    output_folder
) {
  
  pt <- patient$setup$pt_id
  sz <- patient$setup$seizure
  
  dt_gpd <- gpd_bands %>% 
    dplyr::filter(Band == "Simulated signal")
  
  dir <- paste0(output_folder, "gpd-shp-", pt, "-", sz, "-0.PNG")
  png(dir, width = 1080, height = 720)
  plot_gpd_shape(dt_gpd) -> plot
  print(plot)
  dev.off()
  
  dir <- paste0(output_folder, "gpd-scl-", pt, "-", sz, "-0.PNG")
  png(dir, width = 1080, height = 720)
  plot_gpd_scale(dt_gpd) -> plot
  print(plot)
  dev.off()
  
  bands <- c("delta", "theta", "alpha", "beta", "gamma")
  for (i in 1:5) {
    
    dt_gpd <- gpd_bands %>% 
      dplyr::filter(Band == paste0("Filtered ", bands[i], " band"))
    
    dir <- paste0(output_folder, "gpd-shp-", pt, "-", sz, "-", i, ".PNG")
    png(dir, width = 1080, height = 720)
    plot_gpd_shape(dt_gpd) -> plot
    print(plot)
    dev.off()
    
    dir <- paste0(output_folder, "gpd-scl-", pt, "-", sz, "-", i, ".PNG")
    png(dir, width = 1080, height = 720)
    plot_gpd_scale(dt_gpd) -> plot
    print(plot)
    dev.off()
  }  
  
} 


plot_all <- function(
    patient_scenario, 
    patient_simulated, 
    output_folder_plots = "./plot/"
) {
  
  p <- patient_scenario
  p_simu <- patient_simulated
  
  
  output_folder_plots = ifelse(
    is.null(p$tail_factor_place), 
    output_folder_plots,
    paste0("./plot/eeg-simulation/", p$tail_factor_place, "/")
  )
  
  plot_signal_all(p_simu$signal, p_simu$bands, output_folder_plots)
  plot_chi_all(p_simu$signal, p_simu$bands, p$window_duration, upp_quantile = .5, output_folder_plots)
  plot_gpd_all(p_simu$signal, p_simu$bands, p_simu$gpd, output_folder_plots)
  
}