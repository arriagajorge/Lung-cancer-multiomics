#!/usr/bin/env Rscript
setwd("/home/jvasquez/Documents/Omics/LUAD/2_SGCCA/")
# arguments as we need
argms <- list(c)
cont <- 1

for (i in seq(0, 1, 0.01)){ # seq(0, 1, 0.01) = 0.00, 0.01, 0.02, ..., 0.99, 1.00
  argms[[cont]] = rep(i, 3)
  cont <- cont + 1
}

# iterate arguments in the list and execute 2_1FitCommand with the respective arg
for (arg in argms) {
  commands <- paste("Rscript 2_1FitAux.R", paste(arg, collapse = " "))
  system(commands)
  print(arg)
}
