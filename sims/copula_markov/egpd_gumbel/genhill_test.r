source("code/models/copula_markov/simulate.r")
library(ReIns)

genHillk <- function(x, k_max = length(x)) {
    hill <- ReIns::Hill(x)

    genhill <- ReIns::genHill(x, hill$gamma)

    plot(1:k_max, genhill$gamma[1:k_max], type = "l")

    return(genhill)
}

genhill <- genHillk(egpd_data, k_max = 1000)
genhill$gamma[1000]

genhill_gumbel <- genHillk(egpd_gumbel_data$x, k_max = 2000)
genhill_gumbel$gamma[1500]

paretoqq <- function(x) {
    # Step 1: Sort the data in descending order
    x_sorted <- sort(x, decreasing = TRUE)
    n <- length(x_sorted)
    ranks <- 1:n

    # Step 2: Calculate coordinates
    # X-axis: Theoretical log-quantiles -log(i/(n+1))
    # Y-axis: Log of ordered observations
    plot_x <- -log(ranks / (n + 1))
    plot_y <- log(x_sorted)

    # Step 3: Plot
    plot(plot_x, plot_y, 
        main = "Standard Pareto QQ Plot",
        xlab = "-log(i / (n+1))", ylab = "log(X_i,n)",
        pch = 20, col = "darkblue")
    grid()
    # The slope of the rightmost points is the Hill estimate (gamma)
}

genparetoqq <- function(x) {
    x_sorted <- sort(x, decreasing = TRUE)
    # Function to calculate Hill estimates for all k (1 to n-1)
    calc_hill_all <- function(x_sorted) {
    n <- length(x_sorted)
    log_x <- log(x_sorted)
    # Efficiently calculate Hill for all j using cumulative sums
    cumsum_log <- cumsum(log_x)
    hill_estimates <- (cumsum_log[1:(n-1)] / (1:(n-1))) - log_x[2:n]
    return(hill_estimates)
    }

    # 1. Get Hill estimates for all j
    H_j <- calc_hill_all(x_sorted) 

    # 2. Construct UH statistics
    # UH_j is X_j,n * H_j,n (note: we use x_sorted[1:n-1] to align with H_j)
    UH_j <- x_sorted[1:(n-1)] * H_j

    # 3. Coordinates for UH QQ plot
    # We plot log(UH) against the log-ranks
    plot_uh_x <- -log((1:(n-1)) / (n + 1))
    plot_uh_y <- log(UH_j)

    # 4. Plot
    plot(plot_uh_x, plot_uh_y, 
        main = "Generalized (UH) QQ Plot",
        xlab = "-log(j / (n+1))", ylab = "log(UH_j,n)",
        pch = 20, col = "darkred")
    grid()
    # The slope of this plot represents the Generalized Hill estimate
}

paretoqq(egpd_gumbel_data$x)
genparetoqq(egpd_gumbel_data$x)

paretoqq(egpd_data)
genparetoqq(egpd_data)
