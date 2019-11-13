### R TRAINING PRACTICE ####################################
# AUTHOR: ZHUYU
# DATE: 10/7/2019

### NOTES ########################################################


### COMMAND LINE PARAMETERS #########################################  


### PREABMBLE ########################################################
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(car)
library(boot)
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

### FUNCTION_1 ######################################################
#Input variables:
#filename filename to be read
#Output variables:
#output.data vector of sample names
#Description:
#Function that reads ina file and outputs names listed in the file
function_1 = function(filename){
  
# read filename file
  input.data = read.table(
    file = filename,
    header = T,
    sep = ' '
  );
  
# return vector of sample names  
output.data = input.date$name;

return(output.data);
}

## Question 1###################################################

# 1. read file
ahr.test = read.table("/cloud/project/AHR-test-file.txt", header=T);
summary(ahr.test);
#convert data to long format 
ahr.test$number = paste("mouse", 1:9, sep = "");
ahr.test.long=gather(ahr.test, Treatment, test, Control:Treated, 
                     factor_key=TRUE);
#plot induvidual change curve
ggplot(data=ahr.test.long, aes(x=Treatment, y=test, group=number)) + geom_line()

# 2. Perform a t-test between control and treated
# paired t-test
t.test(ahr.test$Control, ahr.test$Treated, paired = T);
# Paired t-test gave the result that the p-value=0.1345>0.05, therefore, we failed to reject the null at the significance level of 0.05
# And we can conclude that there is no significant differnce between the Control group and the Treated group 

# 3. Perform a wilcoxon test between control and treated
# wilconxon signed-rank test, asssume an underlying continuous symmetric distribution
wilcox.test(ahr.test$Control, ahr.test$Treated, paired = T);
# Wilcoxon test gave the result that the p-value=0.1953>0.05

# 4. Calculate a fold-change between control and treated
foldchange = log2(ahr.test$Treated/ahr.test$Control);
hist(foldchange);
# 7 mouses showed an AHR increasement after being treated and 1 mouse decreased. 

### Theory #########################################################################
# 1. different type of ttest
# 1 sample/ paired/ 2 sample(equal variance? var.test(x,y))

# 2. different two sample tests
# paired t test/Wilcoxon signed-rank test(median)
# 2 sample t test/Wilcoxon rank-sum test(median)
# McNemar's Test(two sample test for binomial proportions for matched-pair data) 

# 3. Which test to use and when
# t-test follows the asumption that the sample follows normal distribution,
# when the assumption is violated, like the sample size is too small,use non-parametric



### Question2 #########################################################################

# 1.read file
tumor1 = read.table("/cloud/project/input1.txt", header=T);
tumor2 = read.table("/cloud/project/input2.txt", header=T);

# 2.combine two files
# 2.1 Sort the data and cbind
tumor1.sort = arrange(tumor1,GeneID);
tumor2.sort = arrange(tumor2,GeneID);
tumor.sort.cbind = cbind(tumor1.sort, tumor2.sort)
tumor.sort.cbind = tumor.sort.cbind[,-5]
#2.2 Merge
tumor.merge = merge(tumor1, tumor2, all=T)
#2.3 Verify if they get the same results
table(tumor.sort.cbind == tumor.merge, useNA = 'ifany')
# after removing the repeated gene ID colums in the sort&cbind way, the two ways gave the same results.

# 3. Perform a t-test comparing the first three tumours to the last nine tumours for *each* gene using a for-loo
# create empty dataframe to put p values

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
   
tumor.ttest = emptydf(500, 2, c('character','numeric'), c('GeneID', 'pvalue'));

tumor.ttest$GeneID=tumor.merge$GeneID;

# perform ttest
for (i in 1:nrow(tumor.merge)){
  tumor.ttest[i,2] = tidy(t.test(tumor.merge[i,2:4],tumor.merge[i,5:13]))$p.value;
}  

# 4. histogram of the p values
ggplot(tumor.ttest,aes(pvalue)) + geom_histogram(binwidth = 0.025);

# 5. It's not rotated

# 6. plot pvalues in log space
tumor.ttest$log.pvalue = log(tumor.ttest$pvalue);
ggplot(tumor.ttest, aes(log.pvalue)) + geom_histogram();

# 7. Since log(0.05)=-3, from the plot we can find out that most p values from the ttest
# are greater than 0.05. Therefore, among the 500 genes, there don't exsit a significant 
# mRNA levels differnece between Tumor type A and type B.
nrow(tumor.ttest[tumor.ttest$pvalue<=0.05,]);#33


### Question3 ########################################################################

# 1.Wilcoxon test
# create empty dataframe
tumor.wilcoxon = emptydf(500, 2, c('character','numeric'), c('GeneID', 'pvalue'));

tumor.wilcoxon$GeneID = tumor.merge$GeneID;
# perform wilcoxon rank sum test
for (i in 1:nrow(tumor.merge)){
  tumor.wilcoxon[i,2] = tidy(wilcox.test(unlist(tumor.merge[i,2:4]), unlist(tumor.merge[i,5:13]), 
                          alternative = "two.sided"))$p.value;
}
# 14 warnings show there are tied data

# plot histogram
ggplot(tumor.wilcoxon,aes(pvalue)) + geom_histogram(binwidth = 0.05);

#log transform the p values and plot the histogram
tumor.wilcoxon$log.pvalue = log(tumor.wilcoxon$pvalue);
ggplot(tumor.wilcoxon, aes(log.pvalue)) + geom_histogram();

nrow(tumor.wilcoxon[tumor.wilcoxon$pvalue<=0.05,]);#24


# 3. fold change
tumor.fc=tumor.merge %>%
  mutate(mean.a = (Patient1+Patient2+Patient3)/3,
            mean.b = rowSums(.[5:13])/9,
         log.mean.a = log(mean.a),
         log.mean.b = log(mean.b),
         logfc=log.mean.b - log.mean.a) %>%
  select(GeneID, logfc)

ggplot(tumor.fc, aes(logfc)) + geom_histogram(binwidth = 0.01)

# Both the t test and wilcoxon log p values histgrams shows a log normal trend
# and gave similar results about the p value
# The fold change histgram shows a normal trend with mean around 0. shows that most genes 
# do not have a big difference in mRNA level between type A and B.


# 4. Use apply
# apply(X, MARGIN, FUN)
# -x: an array or matrix
# -MARGIN:  take a value or range between 1 and 2 to define where to apply the function:
# -MARGIN=1`: the manipulation is performed on rows
# -MARGIN=2`: the manipulation is performed on columns
# -MARGIN=c(1,2)` the manipulation is performed on rows and columns
# -FUN: tells which function to apply. Built functions like mean, median, sum, min, max and even user-defined functions can be applied>

### Function_2 ###################################################################################
# Input variables:
# filename test type
# output variables:
# p values
# description:
# function that read in the dataset and test type, then conduct the test and output p values
sample.test = function(df,  test.type){
 
    if('ttest' == test.type ) {result = tidy(t.test(df[2:4], df[5:13]))$p.value};
    if('wilcoxon' == test.type ) {result = tidy(wilcox.test(df[2:4], df[5:13], alternative = "two.sided"))$p.value};
    if('fc' == test.type ) {result = log2(mean(df[5:13])/mean(df[2:4]))};
  BoutrosLab.statistics.general
  return(result);
}

result.ttest = apply(as.matrix(sapply(tumor.merge, as.numeric)), 1, FUN =  sample.test,  test.type='ttest');
result.wilcoxon = apply(as.matrix(sapply(tumor.merge, as.numeric)), 1, FUN =  sample.test,  test.type='wilcoxon');
result.fc = apply(as.matrix(sapply(tumor.merge, as.numeric)), 1, FUN =  sample.test,  test.type='fc');

### Question 4 #############################################################################################
# Adjustment for multiple comparisons are needed to avoid spurious associations.
# MUiltiple comparison can increase type I error
# Bonferroni is conservative, less powerful: alpha* = alpha/(k  2)
Bonferroni.ttest =p.adjust(tumor.ttest$pvalue, method = "bonferroni");
hist(log(Bonferroni.ttest), breaks=50, xlim=c(-1.5,0));

Bonferroni.wilcoxon =p.adjust(tumor.wilcoxon$pvalue, method = "bonferroni");
hist(log(Bonferroni.wilcoxon),breaks = 50, xlim=c(-3,0));

# False Discover Rate controls the number of false positives copared to the total number of positives
# The methods BH (Benjamini–Hochberg, which is the same as FDR in R) and BY control the false discovery rate. 
FDR.ttest =p.adjust(tumor.ttest$pvalue, method = "BH");
hist(log(FDR.ttest), breaks = 50, xlim=c(-3,0));

FDR.wilcoxon =p.adjust(tumor.wilcoxon$pvalue, method = "BH");
hist(log(FDR.wilcoxon), breaks = 50, xlim = c(-1,0));

x=tumor.ttest$pvalue;
y=cbind(Bonferroni.ttest, FDR.ttest);


# After the adjustment, the number of significant results was shrinked to 0. And the Bonferroni correction 
# shows a more big effect on the original p value

### question 5 #############################################################################
# 1. Calculate the median of the first three columns for each gene
random.typea = apply(tumor.merge[,2:4],1,median);

# 2. Use a permutation test to estimate the expected value for each gene:
# a. Randomly select three columns from amongst all 12
set.seed(11)
random.3 = tumor.merge[, sample(2:ncol(tumor.merge),3, replace = F)];

# b. Calculate their median
random.median = apply(random.3,1,median);

# c. Determine if this value is larger or smaller than that of the first 3 columns
median.2 = cbind(random.typea,random.median);

### Function_3 ##############################################################################
# Input variables:
# filename
# output variables:
# number of genes which have larger median than the first 3 columns
# description:
# function that read in the dataset compare the two columns of medians, output 1 if the
# gene has larger median in the first 3 columns, else output 0.

compare = function(df){
  
    if(df[1] > df[2]){
      i=1
    }else if(df[1] <= df[2]){
      i=0
    }
  
  return(i)
}

median.compare = apply(median.2, 1, compare);
median.compare
table(median.compare)

# We can find out among 500 genes, 200 have a smaller median in the first 3 columns.

# d. Repeat 1000 times

### function_4 ###############################################################################
# Input variables:
# filename
# output variables:
# comparison results of two groups of median
# description:
# function that read in the dataset, generate the median of first 3 columns and random 3 columns,
# compare the two columns of medians and if the random group median is bigger, i=1, if not, i=0. 
# The final output is a list of 1 and 0

sample.median = function(df){
  # select 3 random columns from the dataframe
  df1 = df[,sample(2:ncol(df), 3, replace = T)];
  # calculate the median of the 3 values
  median.1 = apply(df1, 1, median);
  # calculate the median of first 3 column
  median.a = apply(df[,2:4],1,median);
  #put the two columns of median together
  df2 = cbind(median.a, median.1);
  # compare the two columns 
  num = apply(df2, 1, compare);
  return(num);
  
}

set.seed(16)
rep = tumor.merge$GeneID;
for (i in 1:1000){
  rep = cbind(rep,sample.median(tumor.merge))
}

# 3. Use the frequencies in 2. to estimate a p-value for each gene
freq = apply(rep[, 2:1001], 1, sum)/1000
hist(freq)

# 4.Perform a false-discovery adjustment on the p-values (?p.adjust)
FDR.freq =p.adjust(freq, method = "BH");
hist(FDR.freq)

# 5.Write your results (gene ID, observed median, expected median, p-value, adjusted p-value) to file in a tab-delimited format
output = cbind(as.character(tumor.merge$GeneID), random.typea)
output = cbind(output, apply(tumor.merge[,2:13],1,median))
output = cbind(output,freq)
output = cbind(output, FDR.freq)
colnames(output) = c("GeneID", "Oberserved Median", "Expected Median", "P-value", "Adjusted P-value");
write.table(output, "/cloud/project/result_q5.txt", append = FALSE, sep = " ", dec = ".",
            col.names = TRUE)

# 6.Plot a histogram of the (unadjusted) p-values. What does this tell you?
hist(freq, breaks = 50)
# h0: ma>mb h1: ma<mb
# among 500 genes, most genes do not have a higher leverl of mRNA in typeA tumor.

# In the case where a set of observations can be assumed to be from an independent and 
# identically distributed population, this can be implemented by constructing a number of 
# resamples with replacement, of the observed dataset (and of equal size to the observed dataset).
# Bootstrap would underestimate the skewed tails. so the distributionwould be more centered than the tru distribution
# Cross-validation evaluates how well an algorithm generalizes, whereas bootstrapping actually helps 
# the algorithm to generalize better.
# Bootstrapping is more about building ensemble models or estimating parameters.


### Question 6 ##########################################################################
# download packages
path = "/cloud/project/";

input = c("BoutrosLab.utilities_1.9.10.tar.gz",
           "BoutrosLab.dist.overload_1.0.2.tar.gz",
           "BoutrosLab.statistics.general_2.1.3.tar.gz",
           "BoutrosLab.statistics.survival_0.4.20.tar.gz",
           "BoutrosLab.prognosticsignature.general_1.3.13.tar.gz",
           "BoutrosLab.plotting.general_5.9.8.tar.gz",
           "BoutrosLab.plotting.survival_3.0.10.tar.gz"
);
library("doParallel")
library("rbenchmark")
library("latticeExtra")
library("hexbin")
library("gridExtra")
library("e1071")

for (i in 1:length(input)){
  install.packages(paste(path, input[i], sep="/"), repo=NULL, type="source", dependencies=TRUE);
  library(unlist(strsplit(input[i], "_"))[1],character.only=TRUE);
}


### Question2 #########################################################################

# 1.read file
tumor1 = read.table("/cloud/project/input1.txt", header=T);
tumor2 = read.table("/cloud/project/input2.txt", header=T);

# 2.combine two files
# 2.1 Sort the data and cbind
tumor1.sort = arrange(tumor1,GeneID);
tumor2.sort = arrange(tumor2,GeneID);
tumor.sort.cbind = cbind(tumor1.sort, tumor2.sort)
tumor.sort.cbind = tumor.sort.cbind[,-5]
#2.2 Merge
tumor.merge = merge(tumor1, tumor2, all=T)
#2.3 Verify if they get the same results
table(tumor.sort.cbind == tumor.merge, useNA = 'ifany')
# after removing the repeated gene ID colums in the sort&cbind way, the two ways gave the same results.

# 3. Perform a t-test comparing the first three tumours to the last nine tumours for *each* gene using a for-loo
# create empty dataframe to put p values

tumor.ttest2 = emptydf(500, 2, c('character','numeric'), c('GeneID', 'pvalue'));

tumor.ttest2$GeneID=tumor.merge$GeneID;

# perform ttest

for (i in 1:nrow(tumor.merge)){
  tumor.ttest2[i,2] = get.ttest.p(tumor.merge[i,],
                               group1 = c(F, T, T, T, rep(F,9)),
                               group2 = c(rep(F,4), rep(T,9)),
                               paired = FALSE,
                               var.equal = FALSE,
                               alternative = 'two.sided');
} 


# 4. histogram of the p values
create.histogram(tumor.ttest2$pvalue, breaks = 50,type = 'count')

# 5. It's not rotated

# 6. plot pvalues in log space
tumor.ttest2$log.pvalue = log(tumor.ttest2$pvalue);
create.histogram(tumor.ttest2$log.pvalue, type = 'count', breaks = 20, ylab.label = 'Counts', xlab.label = 'log(Pvalues)', xlab.cex = 1.5,
                 ylab.cex = 1.5)

# 7. Since log(0.05)=-3, from the plot we can find out that most p values from the ttest
# are greater than 0.05. Therefore, among the 500 genes, there don't exsit a significant 
# mRNA levels differnece between Tumor type A and type B.
nrow(tumor.ttest[tumor.ttest$pvalue<=0.05,]);#33


### Question3 ########################################################################

# 1.Wilcoxon test
# create empty dataframe
tumor.wilcoxon2 = emptydf(500, 2, c('character','numeric'), c('GeneID', 'pvalue'));

tumor.wilcoxon2$GeneID = tumor.merge$GeneID;
# perform wilcoxon rank sum test
for (i in 1:nrow(tumor.merge)){
  tumor.wilcoxon2[i,2] = get.utest.p(as.matrix(sapply(tumor.merge[i,], as.numeric)),  group1 = c(F, rep(TRUE, 3), rep(FALSE, 9)),
                                     group2 = c(rep(FALSE,4), rep(TRUE, 9)),  
                                     paired = FALSE,
                                     alternative = 'two.sided');
}


# 14 warnings show there are tied data

# plot histogram

create.histogram(tumor.wilcoxon2$pvalue, breaks = 20, type = 'count');
#log transform the p values and plot the histogram
tumor.wilcoxon2$log.pvalue = log(tumor.wilcoxon2$pvalue);
create.histogram(tumor.wilcoxon2$log.pvalue, breaks = 20,type = 'count', ylab.label = 'Counts', xlab.label = 'log(Pvalues)', xlab.cex = 1.5,
                 ylab.cex = 1.5);

nrow(tumor.wilcoxon[tumor.wilcoxon$pvalue<=0.05,]);#24


# 3. fold change
tumor.fc2 = emptydf(500, 2, c('character','numeric'), c('GeneID', 'fc'));

for (i in 1:nrow(tumor.merge)){
tumor.fc2[i,] = get.foldchange(as.matrix(sapply(tumor.merge[i,], as.numeric)),  
                               group1 = c(F, rep(TRUE, 3), rep(FALSE, 9)),
                           group2 = c(rep(FALSE,4), rep(TRUE, 9)), 
                           logged = F);
}


create.histogram(tumor.fc2$fc, type = 'count', breaks = 40, ylab.label = 'Counts', xlab.label = 'log(FC)', xlab.cex = 1.5,
                 ylab.cex = 1.5)


# Both the t test and wilcoxon log p values histgrams shows a log normal trend
# and gave similar results about the p value
# The fold change histgram shows a normal trend with mean around 0. shows that most genes 
# do not have a big difference in mRNA level between type A and B.


# 4. Use apply
# apply(X, MARGIN, FUN)
# -x: an array or matrix
# -MARGIN:  take a value or range between 1 and 2 to define where to apply the function:
# -MARGIN=1`: the manipulation is performed on rows
# -MARGIN=2`: the manipulation is performed on columns
# -MARGIN=c(1,2)` the manipulation is performed on rows and columns
# -FUN: tells which function to apply. Built functions like mean, median, sum, min, max and even user-defined functions can be applied>

### Function_2 ###################################################################################
# Input variables:
# filename test type
# output variables:
# p values
# description:
# function that read in the dataset and test type, then conduct the test and output p values

result.ttest2 =apply(as.matrix(sapply(tumor.merge, as.numeric)), 1, 
                   FUN = get.ttest.p, 
                   group1 = c(rep(FALSE, 1), rep(TRUE, 3), rep(FALSE, 9)),
                   group2 = c(rep(FALSE,4), rep(TRUE, 9)),  paired = FALSE, var.equal = FALSE, 
                   alternative = 'two.sided')

result.wilcoxon2 = apply(as.matrix(sapply(tumor.merge, as.numeric)), 1, 
                        FUN = get.utest.p,
                        group1 = c(rep(FALSE, 1), rep(TRUE, 3), rep(FALSE, 9)),
                        group2 = c(rep(FALSE,4), rep(TRUE, 9)),  
                        paired = FALSE,  
                        alternative = 'two.sided'
                        );
result.fc2 = apply(as.matrix(sapply(tumor.merge, as.numeric)), 1, 
                  FUN =  get.foldchange,  
                  group1 = c(rep(FALSE, 1), rep(TRUE, 3), rep(FALSE, 9)),
                  group2 = c(rep(FALSE,4), rep(TRUE, 9)),  
                  logged = F);

### Question 4 #############################################################################################
# Adjustment for multiple comparisons are needed to avoid spurious associations.
# MUiltiple comparison can increase type I error
# Bonferroni is conservative, less powerful: alpha* = alpha/(k  2)
Bonferroni.ttest =p.adjust(tumor.ttest$pvalue, method = "bonferroni");
create.histogram(log(Bonferroni.ttest), type = 'count',breaks=50, ylab.label = 'Counts', xlab.label = 'log(Pvalues)', xlab.cex = 1.5,
                 ylab.cex = 1.5);

Bonferroni.wilcoxon =p.adjust(tumor.wilcoxon$pvalue, method = "bonferroni");
create.histogram(log(Bonferroni.wilcoxon), type = 'count',ylab.label = 'Counts', xlab.label = 'log(Pvalues)', xlab.cex = 1.5,
                 ylab.cex = 1.5);
# False Discover Rate controls the number of false positives copared to the total number of positives
# The methods BH (Benjamini–Hochberg, which is the same as FDR in R) and BY control the false discovery rate. 
FDR.ttest =p.adjust(tumor.ttest$pvalue, method = "BH");
create.histogram(log(FDR.ttest), breaks=50,type = 'count',xlim=c(-3,0), ylab.label = 'Counts', xlab.label = 'log(Pvalues)', xlab.cex = 1.5,
                 ylab.cex = 1.5);

FDR.wilcoxon =p.adjust(tumor.wilcoxon$pvalue, method = "BH");
create.histogram(log(FDR.wilcoxon), type = 'count', breaks=50);

x=tumor.ttest$pvalue;
y=cbind(Bonferroni.ttest, FDR.ttest);


# After the adjustment, the number of significant results was shrinked to 0. And the Bonferroni correction 
# shows a more big effect on the original p value

### question 5 #############################################################################
# 1. Calculate the median of the first three columns for each gene
random.typea = apply(tumor.merge[,2:4],1,median);

# 2. Use a permutation test to estimate the expected value for each gene:
# a. Randomly select three columns from amongst all 12
set.seed(11)
random.3 = tumor.merge[, sample(2:ncol(tumor.merge),3, replace = F)];

# b. Calculate their median
random.median = apply(random.3,1,median);

# c. Determine if this value is larger or smaller than that of the first 3 columns
median.2 = cbind(random.typea,random.median);

### Function_3 ##############################################################################
# Input variables:
# filename
# output variables:
# number of genes which have larger median than the first 3 columns
# description:
# function that read in the dataset compare the two columns of medians, output 1 if the
# gene has larger median in the first 3 columns, else output 0.
compare = function(df){
  
  if(df[1] > df[2]){
    i=1
  }else if(df[1] <= df[2]){
    i=0
  }
  
  return(i)
}

median.compare = apply(median.2, 1, compare);
median.compare
table(median.compare)
# We can find out among 500 genes, 200 have a smaller median in the first 3 columns.

# d. Repeat 1000 times


### function_4 ###############################################################################
# Input variables:
# filename
# output variables:
# comparison results of two groups of median
# description:
# function that read in the dataset, generate the median of first 3 columns and random 3 columns,
# compare the two columns of medians and if the random group median is bigger, i=1, if not, i=0. 
# The final output is a list of 1 and 0

sample.median = function(df){
  # select 3 random columns from the dataframe
  df1 = df[,sample(2:ncol(df), 3, replace = T)];
  # calculate the median of the 3 values
  median.1 = apply(df1, 1, median);
  # calculate the median of first 3 column
  median.a = apply(df[,2:4],1,median);
  #put the two columns of median together
  df2 = cbind(median.a, median.1);
  # compare the two columns 
  num = apply(df2, 1, compare);
  return(num);
  
}

set.seed(16)
rep = tumor.merge$GeneID;
for (i in 1:1000){
  rep = cbind(rep,sample.median(tumor.merge))
}

# 3. Use the frequencies in 2. to estimate a p-value for each gene
freq = apply(rep[, 2:1001], 1, sum)/1000

create.histogram(freq, type = 'count',ylab.label = 'Counts', xlab.label = 'P values', xlab.cex = 1.5,
                 ylab.cex = 1.5)


# 4.Perform a false-discovery adjustment on the p-values (?p.adjust)
FDR.freq =p.adjust(freq, method = "BH");
create.histogram(FDR.freq, type = 'count',ylab.label = 'Counts', xlab.label = 'BH P values', xlab.cex = 1.5,
                 ylab.cex = 1.5)
# 5.Write your results (gene ID, observed median, expected median, p-value, adjusted p-value) to file in a tab-delimited format
output = cbind(as.character(tumor.merge$GeneID), random.typea)
output = cbind(output, apply(tumor.merge[,2:13],1,median))
output = cbind(output,freq)
output = cbind(output, FDR.freq)
colnames(output) = c("GeneID", "Oberserved Median", "Expected Median", "P-value", "Adjusted P-value");
write.table(output, "/cloud/project/result_q5.txt", append = FALSE, sep = " ", dec = ".",
            col.names = TRUE)

# 6.Plot a histogram of the (unadjusted) p-values. What does this tell you?
hist(freq, breaks = 50)
# h0: ma>mb h1: ma<mb
# among 500 genes, most genes do not have a higher leverl of mRNA in typeA tumor.

