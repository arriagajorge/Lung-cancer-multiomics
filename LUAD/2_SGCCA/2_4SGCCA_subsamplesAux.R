# iterate arguments in the list and execute 2_4SGCCA_subsamples.R with the respective arg
setwd("/home/jvasquez/Documents/Omics/LUAD/2_SGCCA/")
subtype = "LUAD"
#pb <- txtProgressBar(min = start, max = end, style = 3)
start = 1
end = 100
for (arg in start:end) {
  #setTxtProgressBar(pb, arg) #update progress bar
  commands = paste(paste("Rscript 2_4SGCCA_subsamples.R", subtype), arg)
  system(commands)
}
 
