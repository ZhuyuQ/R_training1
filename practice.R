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
    if(type[i] == 'numeric') {df[,i] = as.numeric(df[,i])
                               colnames(df)[i] = name[i]};
    if(type[i] == 'character') {df[,i] = as.character(df[,i])
                                colnames(df)[i] = name[i]};
    if(type[i] == 'logical') {df[,i] = as.logical(df[,i])
                            colnames(df)[i] = name[i]};
    if(type[i] == 'factor') {df[,i] = as.factor(df[,i])
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
tumor.fc = tumor.merge %>%
  as_tibble()%>%
  mutate(mean.a = (Patient1+Patient2+Patient3)/3,
            mean.b = rowSums(.[5:13])/9,
         log.mean.a = log2(mean.a),
         log.mean.b = log2(mean.b),
         logfc=log.mean.b - log.mean.a) %>%
  select(GeneID, logfc)

ggplot(tumor.fc, aes(logfc)) + geom_histogram(binwidth = 0.01)


# Both the t test and wilcoxon log p values histgrams shows a log normal trend
#and gave similar results about the p value
# the fold change histgram shows a more normal trend.


