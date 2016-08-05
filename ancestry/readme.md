## Synopsis

This directory contains data files and scripts to estimate genetic admixture for a set of samples. 

## Instructions

###Picking AIMs
pick_aims.py is a script which reads a file called alfred_allele_freqs and picks a set of SNPs for predicting ancestry admixture. It is currently tuned to produce a set of SNPs for global ancestry prediction. To change the balance of SNPs to focus on particular regions, you can edit the pairwise counts in num_per_pair.

This script produces a file called aim_allele_freqs which is used in the next step.

###Predicting Ancestry
predict_ancesty.py is a script which predicts ancestry for a set of samples. To run it on 1000 Genomes Project data for global ancestry prediction run it as follows:

```	
python predict_ancestry.py -t <num_threads> -i aim_allele_freq -g 1kg_genos/geno_data -o 1kg_ancestry_output
```

## Requirements

This code requires Python 2.7, numpy, and scipy