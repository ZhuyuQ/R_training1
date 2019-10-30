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



### Function 1 ##############################################################
##Input variables:
#number of row, number of columns, column data type, and column names
#Output variables:
#empty dataset 
#Description:
#Function that create empty dataframe 

emptydf = function(numrow, numcol, type, name){
  df = data.frame(matrix(NA, nrow=numrow, ncol=numcol));
  for (i in 1:numcol){
    print(type[i])
    if('numeric' == type[i]) {df[,i] = as.numeric(df[,i])
    colnames(df)[i] = name[i]};
    if('character' == type[i]) {df[,i] = as.character(df[,i])
    colnames(df)[i] = name[i]};
    if('logical' == type[i]) {df[,i] = as.logical(df[,i])
    colnames(df)[i] = name[i]};
    if('factor' == type[i]) {df[,i] = as.factor(df[,i])
    colnames(df)[i] = name[i]};
  }
  return(df);
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
loblolly = datasets::Loblolly;
loblolly = reshape(data = loblolly, idvar = "Seed",
                   v.names = "height",
                   timevar = "age",
                   direction = "wide");
loblolly$Seed = as.numeric(loblolly$Seed);
loblolly.sort = loblolly[order(loblolly$Seed),];
loblolly.matrix = as.matrix(loblolly.sort[,-1]);
row.names(loblolly.matrix) = paste("seed", 1:14, sep = ".");
colnames(loblolly.matrix) = c("3yrs", "5yrs", "10yrs", "15yrs", "20yrs", "25yrs");
create.heatmap(
  x = loblolly.matrix,
  # format the colour key
  colourkey.cex = 1,
  colourkey.labels.at = seq(0, 70, 5),
  # set labels to NA -- results in default labels
  xaxis.lab = NA,
  yaxis.lab = NA,
  xaxis.cex = 0.8,
  yaxis.cex = 0.8,

  # set font style (default is bold, 1 is roman)
  xaxis.fontface = 1,
  yaxis.fontface = 1,
  
  # specify clustering method
  # if no clustering is desired, set this to "none"
  clustering.method = "complete",
  # select distance measure
  rows.distance.method = "euclidean",
  cols.distance.method = "manhattan"
);


# 2c. Take a look at the 'ChickWeight' dataset in the R Datasets package
chickweight = datasets::ChickWeight;

# chicken weight by time, spaghetti plot

ggplot(data = chickweight, 
       aes(x=Time, y = weight, group = Chick))+
  geom_line()+ facet_grid(.~ Diet, scales='free')+
  ggtitle("Chicken Weight By Time Under Different Diet, Spaghetti Plot");

# From the spaghetti plot, chicken tend to have a bigger increasement in weight under Diet 3.

  
chick.residual = chickweight %>%
  group_by(Time,Diet) %>%
  mutate(meanweight = mean(weight)) %>%
  ungroup() %>%
  mutate(residual = weight - meanweight) %>%
  group_by(Chick,Diet) %>%
  mutate(median.residual = median(residual)) %>%
  ungroup();

chick.stat = chick.residual %>%
  group_by(Diet) %>%
  mutate(min = min(median.residual),
         q1 = quantile(median.residual,c(.25)),
         q2 = quantile(median.residual,c(.5)),
         q3 = quantile(median.residual,c(.75)),
        max = max(median.residual)) %>%
  ungroup()%>%
  select(Diet,min, q1, q2, q3, max) %>%
  unique();

chick.id = chick.residual %>%
  filter(median.residual %in% as.matrix(chick.stat[,-1]));

ggplot(chick.id, aes(x = Time, y = weight, group = Chick)) + 
  geom_line()+facet_grid(.~ Diet, scales='free') +
  ggtitle("Weight By Time Under Different Diet, Spaghetti Plot",
          subtitle = "min,25%,50%,75%,max of the weight residuals");

# Chicken tend to have higher vairation in weight as tim egoes by under Diet2. 
# Least variation in weight was shown in chicken under diet4.


# chicken weight difference by time, box plot
chickweight.dif = emptydf(nlevels(as.factor(chickweight$Chick)),
                          nlevels(as.factor(chickweight$Time))+2 , 
                      c('character', 'character', rep('numeric',nlevels(as.factor(chickweight$Time)) )), 
                        c('chick', 'diet', paste("dif", 0:nlevels(as.factor(chickweight$Time)) ,sep = ".")));

chickweight$Time = as.character(chickweight$Time);
chickweight.wide = chickweight %>%
  group_by(Chick, Diet) %>%
  spread(Time,weight, fill=NA, sep = ".")

chickweight.wide  = chickweight.wide[, c(1,2,3,9,12,13,14,4,5,6,7,8,10,11)]

# Calculate weight difference between each measurement
for(i in 1: (nlevels(as.factor(chickweight$Time))-1) ){
  chickweight.dif [,i+3]= chickweight.wide[,i+3]-chickweight.wide[,i+2]
}

chickweight.dif[,1:3] = chickweight.wide[,1:3]

chickweight.dif.long = chickweight.dif %>%
  gather(interval,dif,dif.0:dif.11 )


create.violinplot(formula = dif.0 ~ diet, 
               data = chickweight.dif,
               ylimits =  c(35,50)                        
)


ggplot(data = chickweight.dif.long,
       aes(x = diet, y =dif, fill = diet)) + 
  geom_boxplot()+
  facet_grid(.~ interval)

# From the boxplot, we can find out that there is not a much difference in weight between 4 group at time0.
# Diet 1 tends to have the least weight on chicken, while Diet 3&4 have the most changes.
# Diet 3 tends to have the most changes in the last few weeks. 
# In the last few weeks the weight changes tend to have higher vairance in all diet groups.



### Question 3 ######################################################################################
seq.control = read.table("/cloud/project/Q3_SeqControl_data", header=T);


