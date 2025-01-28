# Companion scripts for "Landscape of microRNA and target expression variation and covariation in single mouse embryonic stem cells"

This repository contains custom scripts related to the publication "Landscape of microRNA and target expression variation and covariation in single mouse embryonic stem cells". It also features some pre-processed data tables containing aggregated information about the abundance and variation of miRNAs and their mRNA targets in mESCs, as well as example data for the functions presented.

## Expression matched controls:

When performing enrichment analyses or comparing gene sets, it is important to consider potential biases in the data. For instance, it is quite likely that a gene set of interest is biases towards highly expressed genes in the model system, since lowly expressed genes may simply not be detected in the essay or have high levels of technical noise, overshadowing the effects one aims to observe. It is therefore useful to compare a set of interest with a random control set that shares certain features, for instance the same gene expression distribution.

The `sample_constraint_to_one` fuction samples a random set of genes that are statistically not different for a given metric (e.g., gene expression) from a given set of genes of interest.
```
# 'all' is a list of all genes and the measurement to balance for, list names are gene names (or gene IDs)
# 'subset' is a list of gene names (or IDs)
# 'size' is the number of genes that are to be samples
# 'N' is the number of bins

sample_constraint_to_one <- function(all, subset, size=NA, N=5)
{
  if (is.na(size)){size = length(subset)}
  bin_boundaries = quantile(all, seq(0,N)/N, na.rm=T); bin_boundaries[1] <- bin_boundaries[1]-1
  all_noNA = all[!is.na(all)]
  hits_in_bin_list = rep(NA, N); gene_probabilities = rep(0, length(all)); names(gene_probabilities) = names(all)
  for (i in seq(1,N))
  {
    bin_size = sum((all_noNA>bin_boundaries[i])&(all_noNA<=bin_boundaries[i+1]))
    hits_in_bin = min(sum((all_noNA[subset]>bin_boundaries[i])&(all_noNA[subset]<=bin_boundaries[i+1]))+(1/N), bin_size)
    hits_in_bin_list[i] <- hits_in_bin
    gene_probabilities[names(all_noNA[(all_noNA>bin_boundaries[i])&(all_noNA<=bin_boundaries[i+1])])] = hits_in_bin/bin_size
  }
  gene_probabilities[gene_probabilities==0] <- min(gene_probabilities[gene_probabilities!=0])
  check = F; cutoff = 0.1
  while ((!check))
  {
    test = base::sample(names(gene_probabilities), replace=F, prob=gene_probabilities, size=size)
    if (t.test(all[subset], all[test])$p.value > cutoff){check = T} else{cutoff = cutoff - 0.0001}
  }
  print(paste("succeeded with p-value of ", t.test(all[subset], all[test])$p.value, sep=""))
  return(test)
}
```
The `sample_constraint_to_many` function samples a random set of genes that are statistically not different in multiple given metrics (e.g., gene expression, GC-content and transcript length) from a given set of genes of interest.
```
# 'all' is a matrix of all genes and the measurements to balance for, column names are gene names (or gene IDs)
# 'subset' is a list of gene names (or IDs)
# 'size' is the number of genes that are to be samples
# 'N' is the number of bins

sample_constraint_to_many <- function(all, subset, size, N=5)
{
  probability_matrix = matrix(NA, nrow=nrow(all), ncol=ncol(all))
  for (j in seq(1, nrow(all)))
  {
    bin_boundaries = quantile(all[j,], seq(0,N)/N, na.rm=T); bin_boundaries[1] <- bin_boundaries[1]-1
    all_j = all[j,]; all_j_noNA = all_j[!is.na(all_j)]; subset_noNA = intersect(subset, names(all_j_noNA))
    hits_in_bin_list = rep(NA, N); gene_probabilities = rep(0, length(all_j)); names(gene_probabilities) = names(all_j)
    for (i in seq(1,N))
    {
      bin_size = sum((all_j_noNA>bin_boundaries[i])&(all_j_noNA<=bin_boundaries[i+1]))
      hits_in_bin = min(sum((all_j_noNA[subset_noNA]>bin_boundaries[i])&(all_j_noNA[subset_noNA]<=bin_boundaries[i+1]))+(1/N), bin_size)
      hits_in_bin_list[i] <- hits_in_bin
      gene_probabilities[names(all_j_noNA[(all_j_noNA>bin_boundaries[i])&(all_j_noNA<=bin_boundaries[i+1])])] = hits_in_bin/bin_size
    }
    gene_probabilities[gene_probabilities==0] <- min(gene_probabilities[gene_probabilities!=0])
    probability_matrix[j,] = gene_probabilities
  }
  total_probabilities = rep(1, ncol(all)); names(total_probabilities) = colnames(all)
  for (j in seq(1, nrow(all))) {total_probabilities = total_probabilities * probability_matrix[j,]}
  check = F; cutoff = 0.1
  while ((!check))
  {
    test = base::sample(names(total_probabilities), replace=F, prob=total_probabilities, size=size)
    check = T; for (j in seq(1, nrow(all))) {check = check & (t.test(all[j,subset], all[j,test])$p.value > cutoff)}
    if (!check) {cutoff = cutoff - 0.00001}
  }
  for (i in (1:nrow(all)))
  {print(paste(rownames(input_matrix)[i], ": succeeded with p-value of ", t.test(all[i,subset], all[i,test])$p.value, sep=""))}
  return(test)
}
```
Let's load some data and try them! So we set your working directory and load the example data.
```
setwd("~/Desktop/R_code/")
load(file='gene_expression.RData')
load(file='input_matrix.RData')
```
In this example, we use all genes that start with "Rp" as our gene set of interest. We know that these genes have a different expression distribution than or complete gene set.
```
gene_set_of_interest = names(gene_expression)[substr(names(gene_expression),1,2)=='Rp']
```
We now apply our function to sample 500 genes randomly but in a way that the expression distribution of our sample is statistically not different from the expression distribution of our gene set of interest.
```
sampled_control = sample_constraint_to_one(all=gene_expression, subset=gene_set_of_interest, size=500, N=20)
```
Then we visualize the expression distribution of the whole data set, of our gene set of interest as well as of our sampled gene set.
We can see that our gene set has higher gene expression than our complete data set and that our sampled gene set closely resembles the expression distribution of the gene set of interest.
```
vioplot(log2(gene_expression+1),
        log2(gene_expression[gene_set_of_interest]+1),
        log2(gene_expression[sampled_control]+1),
        col='#fee090',
        names=c('all data', 'subset', 'matched control'),
        ylab='gene expression in log2(RPM+1)\n')
```
<p float="center">
  <img width="500" height="500" src="https://github.com/MarcelTarbier/single-cell_miRNA_variation_and_function/blob/main/images/GitHub_MT_example_1.png">
</p>

Finally, we can check the overlap or our random samples and the gene set of interest. Of our 124 genes of interest only a small fraction is also part of our random sample.
```
length(gene_set_of_interest)
length(sampled_control)
length(intersect(gene_set_of_interest, sampled_control))
```

## miRNA target selection:
WORK IN PROGRESS
