# Load core packages automatically
core_packages <- c("here", "tidyverse", "evd", "copula", "rmarkdown")

for (pkg in core_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

message("Core packages loaded: ", paste(core_packages, collapse = ", "))

options(radian.auto_match = FALSE)
options(radian.auto_indentation = FALSE)
options(radian.complete_while_typing = FALSE)