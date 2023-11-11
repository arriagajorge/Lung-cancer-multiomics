setwd("/home/jorge/Documents/redes")
library(data.table)

##################
#CAPN2 network###
##################

CAPN2_network = read.csv("CAPN2.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality",
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables
CAPN2__df= CAPN2_network[,columns]
#Get mean, sd and combine in a single data frame
CAPN2__mean = colMeans(CAPN2__df)
CAPN2__sd = sapply(CAPN2__df, sd)
CAPN2__summarize = data.frame(CAPN2__mean, CAPN2__sd)

##################
#mir125b network###
##################

mir125b_network = read.csv("mir125b.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality",
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables
mir125b__df= mir125b_network[,columns]
#Get mean, sd and combine in a single data frame
mir125b__mean = colMeans(mir125b__df)
mir125b__sd = sapply(mir125b__df, sd)
mir125b__summarize = data.frame(mir125b__mean, mir125b__sd)

##################
#mir199a2 network###
##################

mir199a2_network = read.csv("mir199a2.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality",
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables
mir199a2__df= mir199a2_network[,columns]
#Get mean, sd and combine in a single data frame
mir199a2__mean = colMeans(mir199a2__df)
mir199a2__sd = sapply(mir199a2__df, sd)
mir199a2__summarize = data.frame(mir199a2__mean, mir199a2__sd)

##################
#flt3 network###
##################

flt3_network = read.csv("mir199a2flt3csv.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality",
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables
flt3__df= flt3_network[,columns]
#Get mean, sd and combine in a single data frame
flt3__mean = colMeans(flt3__df)
flt3__sd = sapply(flt3__df, sd)
flt3__summarize = data.frame(flt3__mean, flt3__sd)


# Mean LUAD
df_mean = data.frame(CAPN2__mean, mir125b__mean, mir199a2__mean, flt3__mean)
meansLUAD=rowMeans(df_mean)
# meansLUAD
# AverageShortestPathLength     BetweennessCentrality       ClosenessCentrality  NeighborhoodConnectivity 
# 1.7589286                 0.1409112                 0.5979097                 5.8247619 
# TopologicalCoefficient 
# 0.2738757 

# sd LUAD
df_sd = data.frame(CAPN2__sd, mir125b__sd, mir199a2__sd, flt3__sd)
sdLUAD=apply(df_sd,1,sd)
# AverageShortestPathLength     BetweennessCentrality       ClosenessCentrality  NeighborhoodConnectivity 
# 0.04874035                0.13867177                0.03458323                1.13526798 
# TopologicalCoefficient 
# 0.18394156
