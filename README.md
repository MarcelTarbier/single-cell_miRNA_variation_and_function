# Companion scripts for ""Landscape of microRNA and target expression variation and covariation in single mouse embryonic stem cells""

This repository contains custom scripts related to the publication "Landscape of microRNA and target expression variation and covariation in single mouse embryonic stem cells". It also features some pre-processed data tables containing aggregated information about the abundance and variation of miRNAs and their mRNA targets in mESCs, as well as example data for the functions presented.

## Expression matched controls:

When performing enrichment analyses or comparing gene sets, it is important to consider potential biases in the data. For instance, it is quite likely that a gene set of interest is biases towards highly expressed genes in the model system, since lowly expressed genes may simply not be detected in the essay or have high levels of technical noise, overshadowing the effects one aims to observe. It is therefore useful to compare a set of interest with a random control set that shares certain features, for instance the same gene expression distribution.

```
# This function samples a random set of genes that are statistically not different for a given metric (e.g., gene expression) from a given set of genes of interest.
#
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
