# tRank:
Tensor-based PageRank (tRank) algorithm is a method for prioritizing disease genes of complex disease, which is by integrating multi-type single cell RNA-seq data and multilayer gene regulatory networks.
# Decription:
In tRank, an n-mode tensor is introduced as a core to represent the multilayer network, then the n-mode multiplication operation between tensor and vector is adapted to implement the PageRank algorithm.
# Steps and preparation:
##(1) Single cell RNA-seq data preprocess
##(2) Multi-layer network construction
##(3) Combining omics data and network data
##(4) Implement tensor-based PageRank algorithm
##(5) Ranking nodes based on their PageRank values
# Case studies description:
We validate our method in type 2 diabetes disease and Alzheimer's disease dataset and we select 3 cell type data each dataset. We download gene regulatory network from Regetwork and use differential mutual information to combine the single cell data and network data. All the data preprocess is done in R with Scanpy and SAVER package. Then, we input them into tRank algoritm in matlab which implemented by ttv.m in MATLAB Tensor Toolbox.
