block_max <- function(x, block_size) {
  # return block maxima given the block size
  sapply(split(x, ceiling(seq_along(x) / block_size)), max)
}

lag_plot <- function(x) {
  n <- length(x)
  
  # Get the name of the object passed
  varname <- deparse(substitute(x))
  
  plot(
    x[-n], x[-1],
    pch = 16,
    col = rgb(0, 0, 1, 0.1),
    main = paste("Lag plot of", varname),
    xlab = bquote(.(as.name(varname))[t-1]),
    ylab = bquote(.(as.name(varname))[t])
  )
  
  abline(0, 1, col = "red")
}
