# Lung Cancer Multi-omics

**To replicate the analysis for both LUAD and LUSC, you need to run the folders in the following order:**

1. PrepoData.
2. SGCCA.
3. Functional_enrich.
4. Files_needed.
5. MI.

In each folder you will find the instructions to execute the respective scripts.


**Methods used in the multiomics analysis**

#PrepoData

1. Bioinformatics Data Query and Download:
GDCquery - Queries Genomic Data Commons (GDC) for transcripts, miRNA expression and DNA Methylation data.

2. Bioinformatics Database Interaction:
TxDb.Hsapiens.UCSC.hg38.knownGene - Loads the TxDb annotation database.
biomaRt::useMart - Connects to the Ensembl database.

3. Ensembl Gene ID Retrieval:
biomaRt::getBM - Retrieves Ensembl transcript IDs for a list of gene symbols.

4. PCA Analysis:
NOISeq::dat - Creates a NOISeq data object for PCA analysis.
explo.plot - Generates exploratory plots for PCA analysis.

5. Exploratory Data Analysis (EDA):
NOISeq::readData - Reads data for use with the NOISeq package.
NOISeq::dat - Prepares data for analysis using NOISeq.

6. Normalization:
betweenLaneNormalization - Normalizes sequencing data between lanes.
NOISeq::filtered.data - Filters low-count features from the data.
tmm - Calculates the TMM normalization factor.
uqua - Calculates upper-quartile normalization.

7. Differential Expression Analysis:
DESeqDataSetFromMatrix - Constructs a DESeqDataSet for differential expression analysis.
estimateSizeFactors - Estimates size factors for normalization.
counts - Retrieves count data from a DESeqDataSet.

8. Batch Effect Correction:
Batch effect correction using NOISeq.

9. Genomic Region Operations:
makeGRangesFromDataFrame - Converts data frames to GRanges objects for genomic region manipulation.


#SGCCA

1. SGCCA (Sparse Generalized Canonical Correlation Analysis):
SGCCA is a method for finding linear combinations of variables (canonical variables) that have maximum correlation across multiple data sets while promoting sparsity.
The script uses the mixOmics library, particularly the wrapper.sgcca function, to perform SGCCA.

2. Wrapper Function (wrapper.sgcca): The script uses a wrapper function that encapsulates the SGCCA implementation and handles the specifics of the analysis, such as penalty setting, scaling, and the number of components.

3. Cross-Validation:
k-fold cross-validation (k=5) is used to assess the performance of the SGCCA model.
The script divides the samples into training and testing sets for each tumor subtype and performs SGCCA on each fold.

4. Bias-Variance Trade-off:
The bias-variance trade-off is associated with the choice of k-fold cross-validation, emphasizing the common choices of k=5 or k=10 for balanced performance.

5. Random Sampling:
Random sampling is used to select a subset of samples (20% of the total) for each fold in the cross-validation.

6. Handling Infinite Values:
Infinite values in the calculated slopes are replaced with NA.

7. Suggested Penalty Selection:
Penalty values corresponding to the maximum slopes are selected as suggested penalties for each omic.

8. AVE (Average Variance Explained): AVE is calculated for the X-blocks in both the initial and final SGCCA results. AVE measures the proportion of variance explained by the canonical components.

9. Penalty Parameter Setting: The penalty parameter is set to control the amount of sparsity in the SGCCA. It determines the trade-off between fitting the data well and having sparse loadings.

#FUNCTIONAL ENRICHMENT

1. Gene Expression Analysis:
Differential gene expression information is obtained from a file named "DE.genes.tsv".

2. Literature Search:
A literature search is performed using the entrez_search function to find relevant literature in the PubMed database based on gene and subtype information.

3. Enrichment Analysis: The script conducts enrichment analysis using two methods:
Gene Ontology (GO) Enrichment Analysis: The enrichGO function from the clusterProfiler library is used to perform GO enrichment analysis on the gene sets.
KEGG Pathway Enrichment Analysis: The enrichKEGG function from the clusterProfiler library is used to perform KEGG pathway enrichment analysis on the gene sets.

4. Exclusive Functions Analysis:
Exclusive functions for both GO Biological Processes (BP) and KEGG pathways are extracted and written to files ("KEGG.exclusive").
The script involves the manipulation of classes and descriptions, merging with external data, and formatting data frames.

5. Gene Ontology (GO) Analysis:
GO analysis is performed using the GO.db and GSEABase libraries. The script retrieves a slim subset of GO terms ("goslim_agr.obo") and counts the number of terms in the Biological Process (BP) ontology for each subtype.

6. Over-Representation Analysis:
A function bias is defined to perform over-representation analysis and identify categories that are significantly over-represented in each subtype. Fisher's exact test is used for statistical testing.

7. Cytoscape Setup:
The script sets the Cytoscape REST API endpoint to "http://localhost:1234"

8. Finding "CROSSLINKED" Functions:
The script calculates the Jaccard index for the components of enriched functions across different biological subtypes.
Hierarchical clustering (hclust) is performed on the Jaccard index to create trees.
Groups of functions that are enriched exactly in the same components are identified.

9. Network Visualization with Cytoscape:
Phylogenetic trees are created from the hierarchical clustering results.
The trees are converted into igraph objects using the ape package.
The igraph object is then converted into a Cytoscape network using the createNetworkFromIgraph function.

10. Identifying Overlapping Functions:
The script identifies functions that map the same genes for a pair of datasets. It calculates the Jaccard index for each pair of datasets and functions.

11. Jaccard Index Calculation:
The jacc function calculates the Jaccard index given two sets.

12. Intersect Functions:
The intersect_functions function computes the Jaccard index for each pair of datasets and functions, creating a matrix of Jaccard indices.

#MI

1. ARACNE multicore:
ARACNE multicore is primarily used for reconstructing gene regulatory networks. It infers relationships between genes based on their expression patterns across different conditions or samples.
ARACNE employs an information-theoretic approach, specifically mutual information, to measure the dependence between pairs of genes. It then uses statistical tests (Data Processing Inequality) to filter out indirect associations, focusing on direct interactions.

2. Kolmogorov-Smirnov test (KS):
The KS test is a non-parametric test that compares the cumulative distributions of two datasets to assess whether they are likely to come from the same underlying distribution.
The KS test is used in the script to compare the distribution of Mutual Information (MI) values for different types of interactions.

3. Cytoscape Integration:
Integrates with Cytoscape using RCy3 for creating and visualizing networks.

#Cox model

1. The Cox proportional hazards model is employed to quantify the relationship between gene expressions and the risk of death, while the risk proportionality test evaluates the assumption that the hazards are proportional over time.

