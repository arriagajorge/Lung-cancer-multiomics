setwd("~/lungsquamouscells/MI/filteredCys/")
library(data.table)

##################
#EIF4G1 network###
##################

EIF4G1_network = read.csv("/home/mdiaz/lungsquamouscells/MI/filteredCys/EIF4G1_LUSC-Table.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality", 
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables 
EIF4G1_df= EIF4G1_network[,columns]
#Get mean, sd and combine in a single data frame 
EIF4G1_mean = colMeans(EIF4G1_df)
EIF4G1_sd = sapply(EIF4G1_df, sd)
#EIF4G1__summarize = data.frame(EIF4G1__mean, EIF4G1__sd)


##################
#TRIO network###
##################

TRIO_network = read.csv("TRIO_LUSC-Table.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality", 
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables 
TRIO_df= TRIO_network[,columns]
#Get mean, sd and combine in a single data frame 
TRIO_mean = colMeans(TRIO_df)
TRIO_sd = sapply(TRIO_df, sd)

##################
#RPS18 network###
##################

RPS18_network = read.csv("RPS18_LUSC-Table.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality", 
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables 
RPS18_df= RPS18_network[,columns]
#Get mean, sd and combine in a single data frame 
RPS18_mean = colMeans(RPS18_df)
RPS18_sd = sapply(RPS18_df, sd)

##################
#LRP1 network###
##################

LRP1_network = read.csv("LRP1_LUSC-Table.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality", 
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables 
LRP1_df= LRP1_network[,columns]
#Get mean, sd and combine in a single data frame 
LRP1_mean = colMeans(LRP1_df)
LRP1_sd = sapply(LRP1_df, sd)

##################
#RPS6 network###
##################

RPS6_network = read.csv("RPS6_LUSC-Table.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality", 
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables 
RPS6_df= RPS6_network[,columns]
#Get mean, sd and combine in a single data frame 
RPS6_mean = colMeans(RPS6_df)
RPS6_sd = sapply(RPS6_df, sd)


##################
#PFN2 network###
##################

PFN2_network = read.csv("PFN2_LUSC-Table.csv")
#Define your variables of interest
columns=c("AverageShortestPathLength", "BetweennessCentrality", "ClosenessCentrality", 
          "NeighborhoodConnectivity", "TopologicalCoefficient")
#Filter according to variables 
PFN2_df= PFN2_network[,columns]
#Get mean, sd and combine in a single data frame 
PFN2_mean = colMeans(PFN2_df)
PFN2_sd = sapply(PFN2_df, sd)


df_mean = data.frame(EIF4G1_mean, LRP1_mean, PFN2_mean, RPS18_mean, RPS6_mean, TRIO_mean)
meansLUSC=rowMeans(df_mean)
# meansLUSC
# AverageShortestPathLength     BetweennessCentrality       ClosenessCentrality  NeighborhoodConnectivity    TopologicalCoefficient 
# 2.1093558                 0.1776303                 0.5165609                 3.9211216                 0.2386346 

df_sd = data.frame(EIF4G1_sd, LRP1_sd, PFN2_sd, RPS18_sd, RPS6_sd, TRIO_sd)
sdLUSC=apply(df_sd, 1, sd)
# sdLUSC
# AverageShortestPathLength     BetweennessCentrality       ClosenessCentrality  NeighborhoodConnectivity    TopologicalCoefficient 
# 0.13254688                0.08673272                0.04303972                0.81222930                0.14513460 