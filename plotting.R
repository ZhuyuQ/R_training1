### R TRAINING PRACTICE ####################################
# AUTHOR: ZHUYU
# DATE: 10/17/2019

### NOTES ########################################################


### COMMAND LINE PARAMETERS #########################################  


### PREABMBLE ########################################################
library("tidyverse");
library("ggplot2");
library("dplyr");
library("tidyr");
library("broom");
library("car");
library("boot");
library("doParallel");
library("rbenchmark");
library("latticeExtra");
library("hexbin");
library("gridExtra");
library("e1071");
library("datasets");
path = "/cloud/project/";
input = c("BoutrosLab.utilities_1.9.10.tar.gz",
          "BoutrosLab.dist.overload_1.0.2.tar.gz",
          "BoutrosLab.statistics.general_2.1.3.tar.gz",
          "BoutrosLab.statistics.survival_0.4.20.tar.gz",
          "BoutrosLab.prognosticsignature.general_1.3.13.tar.gz",
          "BoutrosLab.plotting.general_5.9.8.tar.gz",
          "BoutrosLab.plotting.survival_3.0.10.tar.gz"
);
for (i in 1:length(input)){
  
  library(unlist(strsplit(input[i], "_"))[1],character.only=TRUE);
}


### Question 1 ##################################################################################



### Question 2 #######################################################################################
# 2a. Create a simple scatterplot using BoutrosLab.plotting.general
cars = datasets::cars;
summary(cars);

create.scatterplot(
  formula = dist ~ speed,
  data = cars,
  main = "Scatter plot of distance and speed",
  xlab.label = "Distance",
  ylab.label = "Speed",
  yat = seq(0,150,30),
  xaxis.cex = 0.6,
  yaxis.cex = 0.6,
  xlab.cex = 1,
  ylab.cex = 1,
  main.cex = 1,
  pch  =  1,  
  col  =  "black",
  type  =  c  (  "p"  ,  "g"  ,  "r"  ),
)

# 2b. Create a heatmap displaying data found in the 'Loblolly' dataset 
Loblolly = datasets::Loblolly







### Question 2 ######################################################################################



