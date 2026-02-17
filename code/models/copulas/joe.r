# ------------------------------
# Joe Copula Object
# ------------------------------
copula_joe <- list(

  name = "joe",

  # --------------------------
  # Generator and inverse generator
  # --------------------------
  psi = function(t, theta) 1 - (1 - exp(-t))^(1/theta),   # generator
  inv_psi = function(u, theta) -log(1 - (1 - u)^theta),

)
