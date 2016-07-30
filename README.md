# EECS442 Causal Learning from Data Final Project
Project Title: Causal Inference to Ascertain Causes of Metastasis in Melanoma

Spring 2015
Professor Andy Podgurski

Goal of this project was to attempt to find clinical features that cause metastasis in melanoma using the clinical data from the TCGA dataset. 

Method used to generate causal directed graph was from PC algorithm in R package pcalg. The IDA algorithm (again, implemented in pcalg) was used to identify features in the data with high causality. Finally, counterfactuals and regression were used to solidify these features.

mergingDataSources.R starts with the TCGA clinical data and builds a data frame that incorporates three different sources - New Tumor Event, Follow Up, and Patient - from the clinical dataset. The example shown is for the KIRC (kidney) cancer data; however the process was the same for the data used.

generateGraph.R starts with the output from mergingDataSources and generates a causal graph using the PCalg package on this data. It outputs:
  graph.RData - the full graph generated from the data
  graph2.RData - a smaller subgraph generated from the causal graph that represented some interesting edges
  skcmSub.RData - the data from the subgraph
  
buildModels starts with the output from the previous script and generates the causal models used in the analysis. Only simple marginal structural models and instrumental variable analysis are used.

CETest was a first test to determine causal effects.

The write-up for this analysis is available upon request.
