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

#2.2 Merge
tumor.merge = merge(tumor1, tumor2, all=T)





