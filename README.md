# scvi-multiome
This repo was created to process single-cell multiome (RNA + ATAC) data using Scanpy and Seurat. It represents an end-to-end workflow that will ingest counts and use scvi and poissonvi to produce a joint representation of the data using the muon implementation of weighted-nearest neighbors (WNN). Fragment counts for the ATAC data are generated using the Signac package in R. It is recommended that users correct for ambient RNA using tools such as SoupX or CellBender prior to running the pipeline. Workflow is managed using Snakemake and Anaconda environments are described using yaml files. Links to pacakges and tools used are below:


### scVI

scVI (single-cell Variational Inference) is a scalable probabilistic framework designed to analyze single-cell RNA sequencing (scRNA-seq) data. Built on variational autoencoders, scVI offers powerful tools for batch effect correction, dimensionality reduction, and data integration across multiple samples or experiments. It enables efficient processing of large scRNA-seq datasets, supports flexible modeling of gene expression, and provides biologically interpretable latent spaces for downstream analysis such as clustering, differential expression, and cell type annotation.

### poissonVI

PoissonVI is an extension of the scVI framework specifically designed for single-cell ATAC sequencing (scATAC-seq) data. Unlike scVI, which models gene expression counts using a negative binomial distribution, PoissonVI leverages the Poisson distribution to more accurately model sparse, discrete peak accessibility counts characteristic of scATAC-seq data. This allows PoissonVI to provide more precise analysis of chromatin accessibility and better capture the underlying structure of epigenetic data for tasks such as clustering, dimensionality reduction, and identifying regulatory elements.

### WNN

WNN (Weighted Nearest Neighbors), implemented in the muon framework, is a powerful method for integrating multiple single-cell modalities such as RNA and ATAC sequencing. WNN leverages information from each modality by calculating separate nearest neighbor graphs for each data type and then weighting these graphs based on the relative contribution of each modality to a cell's overall identity. This approach allows for robust multimodal integration, enabling researchers to capture complementary information from both gene expression and chromatin accessibility, or other paired modalities, leading to improved cell type annotation, clustering, and interpretation of complex cellular states in single-cell multiome datasets.
