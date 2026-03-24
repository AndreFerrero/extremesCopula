###################################################################
# Script that takes the initial data text file and outputs a two  #
# dimensional data matrix with consecutive pairs. Once these have #
# been obtained then the data are plotted to obtain figure 2 in   #
# Section 5.3                                                     #
###################################################################

# Read in data from text file
hw.data <- read.csv("winter2016/data_winter2016.csv", header = T)

years <- seq(1937, 2013, by = 1)
new.date <- numeric(0)

# The date is given in the form yyyymmdd where the example 19370101
# would imply 1st January 1937. The aim is to look solely at the
# summer months and as such it is necessary to look at months from
# yyyy0601-yyyy0831. To do this I create a new vector 'newdate' by
# removing the year*10000. In this way only pick the rows that are
# associated with the summer months.

for (j in 1:length(hw.data$DATE)) {
  for (i in 1:(length(years) - 1)) {
    if (hw.data$DATE[j] >= (years[i] * 10000) && hw.data$DATE[j] < (years[i + 1] * 10000)) {
      new.date <- c(new.date, (hw.data$DATE[j] - (years[i] * 10000)))
    }
  }
}

hw.data.mat <- cbind(hw.data, new.date)

JJA.data <- hw.data.mat[which(hw.data.mat$new.date > 531 & hw.data.mat$new.date < 901), ]

# Choose to use data commencing on 1946 to avoid the missing data
# at the start of the series
JJA.data <- JJA.data[JJA.data$DATE > 19460000, ]

# Obtain two dimensional data by taking consecutive pairs
# Remove any overlapping rows where comparing
rem.nos <- seq(1, length(JJA.data$TX * 0.1) + 1, by = 92)

xdat <- c(NA, JJA.data$TX * 0.1) # Add empty data point to shift data for
ydat <- c(JJA.data$TX * 0.1, NA) # comparison of times t and t+1
data <- cbind(xdat, ydat)
j.data <- data[-rem.nos, ] # start and end of season in different years

rem.xdat <- which(j.data[, 1] < -900) # Remove any missing values which are
j.data <- j.data[-rem.xdat, ] # defined as -999.9
rem.xdat <- which(j.data[, 2] < -900)
j.data <- j.data[-rem.xdat, ]

JJA.data.clean <- JJA.data[JJA.data$TX * 0.1 > -100, ]

data <- JJA.data.clean$TX * 0.1

# Remove any unnecessary variables from the workspace
rm(xdat, ydat, data, rem.xdat, hw.data.mat)
