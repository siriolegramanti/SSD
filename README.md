# SSD

This repository is associated with the article **A new class of nonparametric tests for second-order stochastic dominance based on the Lorenz P-P plot** ([arXiv](https://arxiv.org/abs/2308.00317)), by Tommaso Lando and Sirio Legramanti, and contains the code to reproduce the simulation studies in the paper.

Namely, the repository contains:
- [`ssd_source.R`](): the source code to perform the tests proposed in the paper and their considered competitor (see paper);
- [`ssd_experiments.R`](): the main file to reproduce the results in Table 6a of the paper. Note that Table 6b is obtained by switching samples a and b, while Tables 7 and 8 are obtained by changing the data-generating distribution parameters. The other tables in the paper can be obtained by changing the data-generating process and, in case of dependent samples, by setting "dependence=T" in 'ssd_exp()'.
