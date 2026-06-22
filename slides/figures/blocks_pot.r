# =========================================================
# Separate EVT illustrations:
# 1. Block Maxima
# 2. Peaks Over Threshold (POT)
# =========================================================

set.seed(123)

# ---------------------------------------------------------
# Simulate data
# ---------------------------------------------------------
n <- 240
x <- arima.sim(model = list(ar = 0.6), n = n)

# Add extremes
extreme_idx <- c(40, 85, 120, 170, 210)
x[extreme_idx] <- x[extreme_idx] + c(5, 6, 4, 7, 5)

t <- 1:n

# =========================================================
# 1. BLOCK MAXIMA PLOT
# =========================================================

block_size <- 60
n_blocks <- n / block_size

png("slides/figures/block_maxima.png",
    width = 1200,
    height = 500,
    res = 150)

par(mar = c(4,4,3,1))

plot(t, x,
     type = "l",
     lwd = 1.5,
     xlab = "",
     ylab = "Observations",
     main = "Block Maxima")

# Block colors
cols <- c(rgb(0,0,1,0.08),
          rgb(0,1,0,0.08),
          rgb(1,0.5,0,0.08),
          rgb(0.6,0,0.8,0.08))

# Shade blocks
for(i in 1:n_blocks){

  xmin <- (i-1)*block_size + 1
  xmax <- i*block_size

  rect(xmin,
       par("usr")[3],
       xmax,
       par("usr")[4],
       col = cols[i],
       border = NA)
}

# Replot series
lines(t, x, lwd = 1.5)

# Vertical separators
abline(v = seq(block_size, n-block_size, by = block_size),
       lty = 2,
       col = "gray40")

# Find maxima
block_max_t <- numeric(n_blocks)
block_max_x <- numeric(n_blocks)

for(i in 1:n_blocks){

  idx <- ((i-1)*block_size + 1):(i*block_size)

  j <- idx[which.max(x[idx])]

  block_max_t[i] <- j
  block_max_x[i] <- x[j]
}

# Plot maxima
points(block_max_t,
       block_max_x,
       pch = 19,
       cex = 1.2,
       col = c("blue","forestgreen","orange","purple"))

dev.off()

# =========================================================
# 2. POT PLOT
# =========================================================

threshold <- quantile(x, 0.90)

png("slides/figures/pot_threshold.png",
    width = 1200,
    height = 500,
    res = 150)

par(mar = c(4,4,3,1))

plot(t, x,
     type = "l",
     lwd = 1.5,
     xlab = "",
     ylab = "Observations",
     main = "Peaks Over Threshold (POT)")

# Threshold line
abline(h = threshold,
       col = "red",
       lty = 2,
       lwd = 2)

# Exceedances
exc <- which(x > threshold)

# Highlight exceedances
points(exc,
       x[exc],
       pch = 19,
       col = "red",
       cex = 1.2)

# Draw exceedance segments
segments(exc,
         threshold,
         exc,
         x[exc],
         col = "red",
         lty = 2)

# Label threshold
text(10,
     threshold + 1,
     labels = expression(u),
     col = "red",
     cex = 1.2)

dev.off()

# =========================================================
# Files produced:
# - block_maxima.png
# - pot_threshold.png
# =========================================================