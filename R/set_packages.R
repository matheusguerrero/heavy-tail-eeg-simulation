#### Loading packages ####
pcks <- c("tidyverse", "ggthemes", "ggfortify", "gridExtra", 
          "cowplot", "wesanderson", "png", "grid",
          "lazyeval", "ggsci", "RColorBrewer", "scales",
          "eva", "texmex", "TSA", "signal", "ramify", "skmeans", "oce")
if (!require("librarian")) install.packages("librarian")
librarian::shelf(pcks, quiet = TRUE)


##### Basic helper functions ####
## Helper function to pass to the filter() function in the dyplr package.
`%notin%` <- Negate(`%in%`)


#### Basic setup ####
eeq_freq <- 100
eeg_seizure <- 35000
eeg_end <- 50000