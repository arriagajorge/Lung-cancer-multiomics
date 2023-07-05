setwd("/home/mdiaz/lungsquamouscells/sgcca")
subtype = "LUSC"
#pb <- txtProgressBar(min = start, max = end, style = 3)
start = 1
end = 100
for (arg in start:end) {
  #setTxtProgressBar(pb, arg) #update progress bar
  commands = paste(paste("Rscript 2_4sgcca_subsamples.R", subtype), arg)
  system(commands)
}
