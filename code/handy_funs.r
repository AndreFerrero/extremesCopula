block_max <- function(x, block_size) {
  # return block maxima given the block size
  sapply(split(x, ceiling(seq_along(x) / block_size)), max)
}

lag_plot <- function(x, main = NULL) {
  n <- length(x)
  
  plot(
    x[-n], x[-1],
    pch = 16,
    main = main,
  )
  
  abline(0, 1, col = "red")
}
