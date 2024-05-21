# **sciNMF: Single-Cell Integration by  Non-Negative Matrix Factorization**

sciNMF is an R package designed for exploring the heterogeneity of cellular transcriptional states across individuals using single-cell RNA-seq data. This package is developed by HuangLab at the School of Life Sciences, Xiamen University, Xiamen, Fujian, China.

## **Overview**
![alt text](overview.png)
There are 3 main steps to identify cell states from scRNA-seq data:  

**Step 1:** Perform multiple ranks Non-Negative Matrix Factorization (NMF) on single-cell gene expression matrices for each individual.

**Step 2:** Filter out low-quality programs by considering their Interquartile Range (IQR) and median usages. Identify robust programs based on their intra- and inter-sample reproducibility.

**Step 3:** Cluster the robust programs based on their overlapping gene numbers. Generate meta-programs from the clustering results, utilizing the top genes with the highest average weight to represent these meta-programs.


## **Installation**
To install the sciNMF package, you can use the following commands: 
```
# Install devtools if not already installed 
if(!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools") 
} 
devtools::install_github('Tang-RH/sciNMF')
```

## **Tutorial**
For detailed instructions and examples, please refer to our comprehensive tutorial.[https://github.com/Tang-RH/sciNMF/tree/master/Tutorial](https://github.com/Tang-RH/sciNMF/tree/master/Tutorial)

## **Citation**
If you find sciNMF useful for your research, please consider citing our publication.
