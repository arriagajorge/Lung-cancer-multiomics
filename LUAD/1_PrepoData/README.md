# Information pre-processing.

To replicate the analysis run in the following order.

1. Execute *1_1getData.R*. The main function of this file is to load the data. The file will have as main output the file $texttt{subtypeLUAD.tsv}$.

2. Execute *1_2Prepo-mRNA.R*. The main function of this file is to process the mRNA data. The file will have as its main output the file $texttt{RNAseqnormalized.tsv}$. 
**Warning** Verify that the resulting matrix ($\texttt{RNAseqnormalized}$) has no rows of zeros.

3. Run *1_3prepo-miRNA.R*. The main function of this file is to process the miRNA data. The file will have as main output the file $texttt{miRNAseqNormi.tsv}$.

4. Run *1_4prepoMethy.R*. The main function of this file is to process the Methylation data. The file will have as main output the file $texttt{MethyM.tsv}$.

5. Execute *1_5CONCAT_FINAL.R*. The main function of this file is to concatenate the previously obtained information and separate it by subtypes. The main outputs are $$$texttt{normal.MTRX}$, $$texttt{prox.-inflam.MTRX}$, $$texttt{prox.-prolif..MTRX}$, $$texttt{TRU.MTRX}$.

6. Execute *1_6mfanormi.R*, this will give as output the file $$texttt{LUAD.eigenNormi}$.
To replicate the analysis run in the following order.

**Warning** The libraries $texttt{NOISeq}$ and $texttt{data.table}$ have certain conflicts with the functions $texttt{dat}, \texttt{ReadData}$ so it is recommended to verify that before each of these functions there is either $texttt{NOISeq::}$ or $texttt{data.table::}$ before executing the respective function.
