# sciNMF

sciNMF (Single-Cell Individual Non-Negative Matrix Factorization) is an R package designed for exploring the heterogeneity of cellular transcriptional states across individuals using single-cell RNA-seq data.  
The package is developed by HuangLab at Schoolof Life Sciences, XiamenUniversity, Xiamen, Fujian, China.  

There are 3 main steps to identify cell states:  

**Step1**  Performs multiple ranks NMF on single-cell gene expression matrix for each individual.  

**Step2** Filter out low-quality programs by considering their Interquartile Range (IQR) and median usages. Subsequently, identify robust programs based on their intra- and intersample reproducibility.

**Step3** Cluster the robust programs based on their overlapping gene numbers. Subsequently, generate meta-programs from the clustering results, utilizing the top genes with the highest average weight to represent these meta-programs.

The tutorial can be found at:
    [https://github.com/Tang-RH/sciNMF/tree/master/Tutorial](https://github.com/Tang-RH/sciNMF/tree/master/Tutorial)  
