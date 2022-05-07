rm(list = ls())

source("./R/set_packages.R")
source("./R/fun_main.R")
source("./R/fun_cluster.R")

get_extdata <- function(df_eeg, threshold) {
  
  # Function to transform data to Fréchet scale
  Frechettrans <- function(x) 1 / (1 - ecdf(x)(x) * length(x) / (length(x) + 1)) 
  
  # Function to compute the row's norm
  norm_vec <- function(x) sqrt(sum(x ^ 2)) 
  
  # Transforming data to Fréchet scale and computing the norm
  df_eeg <- apply(df_eeg, 2, Frechettrans)
  norms <- apply(df_eeg, 1, norm_vec)
  
  # Finding the exceedances
  index_selection <- (norms > quantile(norms, threshold))
  
  # Filtering extremes in the dataset
  eeg.ext <- df_eeg[index_selection, ]
  norms <- apply(eeg.ext, 1, norm_vec)
  eeg.ext <- eeg.ext / norms 
  
  return(eeg.ext)
  
}

p1 <- read_rds("./output/in_z/patient-1-1.RDS")
p1_4cluster <- get_extdata(abs(p1$signal$eeg[, -c(1, 2, 3)]), threshold = 0.1)

clusterMeans(p1_4cluster,k = 3)
clusterPC(p1_4cluster,k = 3)