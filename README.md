# CausalLearning442
EECS442 Causal Learning Final Project

Attempted to find clinical features that cause metastasis in melanoma using the clinical data from the TCGA dataset.

Method used to generate causal directed graph was from PC algorithm in R package pcalg. 

The IDA algorithm (again, implemented in pcalg) was used to identify features with high causality.

Finally, counterfactuals and regression were used to solidify these features.
