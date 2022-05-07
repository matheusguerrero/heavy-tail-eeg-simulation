rm(list = ls())

## Reading each channel and storing in a matrix
## Note: We are loading only 21 channels. 
##       The ones we have in the figure.
number_of_channels <- 21
eeg_raw_matrix <- matrix(NA, ncol = number_of_channels, nrow = 50000, byrow = T)
for (i in 1:number_of_channels) {
  if (i < 10) {
    path_name <- paste0("./input/ch0", i, ".inp")
  } else {
    path_name <- paste0("./input/ch", i, ".inp")
  }
  eeg_raw_matrix[, i] <- as.matrix(
    read.table(path_name, header = FALSE))
} 

## Passing the correct names of each EEG channel
channels_names <- c("Fp1", "Fp2", "F3", "F4" , "C3",  "C4", "P3", 
                    "P4",  "O1",  "O2", "Sp1", "Sp2", "F7", "F8", 
                    "T3",  "T4",  "T5", "T6",  "Fz",  "Cz", "Pz")
colnames(eeg_raw_matrix) <- channels_names

## Removing Sp1 and Sp2
eeg_matrix <- eeg_raw_matrix[, -c(11,12)]

## Centering each channel
eeg_matrix <- scale(eeg_matrix, center = T, scale = T)

## Basic setup
eeq_freq <- 100
eeg_seizure <- 35000
eeg_end <- 50000
eeg_start <- eeg_seizure - (eeg_end - eeg_seizure) + 1

eeg_reduced <- eeg_matrix[eeg_start:eeg_end, ]
seizure_at <- 15000


## Saving data in a .RDS file
saveRDS(eeg_matrix, file = "./input/original-eeg-data.RDS")
saveRDS(eeg_reduced, file = "./input/original-eeg-reduced.RDS")