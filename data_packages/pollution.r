library(texmex)
library(extRemes)

data(summer)
head(summer)

data(winter)
head(winter)

par(mfrow = c(1, 2))
hist(summer$SO2, breaks = 20, main = "Histogram of SO2", xlab = "Winter SO2 concentration")
hist(winter$SO2, breaks = 20, main = "Histogram of SO2", xlab = "Summer SO2 concentration")

hist(winter$NO2, breaks = 20, main = "Histogram of NO2", xlab = "Winter NO2 concentration")
hist(summer$NO2, breaks = 20, main = "Histogram of NO2", xlab = "Summer NO2 concentration")

hist(winter$O3, breaks = 20, main = "Histogram of O3", xlab = "Winter O3 concentration")
hist(summer$O3, breaks = 20, main = "Histogram of O3", xlab = "Summer O3 concentration")

hist(winter$PM10, breaks = 20, main = "Histogram of PM10", xlab = "Winter PM10 concentration")
hist(summer$PM10, breaks = 20, main = "Histogram of PM10", xlab = "Summer PM10 concentration")

par(mfrow = c(1, 1))
