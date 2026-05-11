# =========================================================
# Declustering illustration
#
# Idea:
# - exceedances occur in temporal clusters
# - keep only the cluster maximum
# - discard dependent exceedances nearby
# =========================================================

set.seed(123)

# ---------------------------------------------------------
# Simulate dependent series
# ---------------------------------------------------------
n <- 240

x <- arima.sim(model = list(ar = 0.8), n = n)

# Create clusters of extremes
cluster_centers <- c(45, 95, 150, 205)

for(c in cluster_centers){

  idx <- (c-2):(c+2)

  x[idx] <- x[idx] + c(2,4,6,4,2)
}

t <- 1:n

# ---------------------------------------------------------
# Threshold
# ---------------------------------------------------------
u <- quantile(x, 0.92)

exc <- which(x > u)

# =========================================================
# Simple declustering rule
# If exceedances are within r observations,
# they belong to the same cluster
# =========================================================

r <- 6

clusters <- list()

current_cluster <- exc[1]

for(i in 2:length(exc)){

  if(exc[i] - exc[i-1] <= r){

    current_cluster <- c(current_cluster, exc[i])

  } else {

    clusters[[length(clusters)+1]] <- current_cluster
    current_cluster <- exc[i]
  }
}

clusters[[length(clusters)+1]] <- current_cluster

# ---------------------------------------------------------
# Cluster maxima
# ---------------------------------------------------------
cluster_max_t <- numeric(length(clusters))
cluster_max_x <- numeric(length(clusters))

for(i in seq_along(clusters)){

  idx <- clusters[[i]]

  j <- idx[which.max(x[idx])]

  cluster_max_t[i] <- j
  cluster_max_x[i] <- x[j]
}

# =========================================================
# Plot
# =========================================================

png("slides/figures/declustering.png",
    width = 1200,
    height = 500,
    res = 150)

par(mar = c(4,4,3,1))

plot(t, x,
     type = "l",
     lwd = 1.5,
     xlab = "",
     ylab = "Observations",
     main = "Declustering of Extremes")

# ---------------------------------------------------------
# Threshold
# ---------------------------------------------------------
abline(h = u,
       col = "red",
       lty = 2,
       lwd = 2)

# ---------------------------------------------------------
# Shade clusters
# ---------------------------------------------------------
for(cl in clusters){

  rect(min(cl)-1,
       par("usr")[3],
       max(cl)+1,
       par("usr")[4],
       col = rgb(0,0,0,0.05),
       border = NA)
}

# Replot series
lines(t, x, lwd = 1.5)

# ---------------------------------------------------------
# All exceedances
# ---------------------------------------------------------
points(exc,
       x[exc],
       pch = 19,
       col = "gray50",
       cex = 1)

# ---------------------------------------------------------
# Cluster maxima kept after declustering
# ---------------------------------------------------------
points(cluster_max_t,
       cluster_max_x,
       pch = 19,
       col = "red",
       cex = 1.2)

# ---------------------------------------------------------
# Optional segments from threshold
# ---------------------------------------------------------
segments(exc,
         u,
         exc,
         x[exc],
         col = "gray70",
         lty = 3)

# ---------------------------------------------------------
# Labels
# ---------------------------------------------------------
text(10,
     u + 0.8,
     labels = expression(u),
     col = "red",
     cex = 1.2)

dev.off()