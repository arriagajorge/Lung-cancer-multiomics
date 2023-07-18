# Preprocesamiento de la información.

Para replicar el análists ejecutar en el siguiente orden.

1. Ejecutar *1_1getData.R*. La principal función de este archivo es cargar la información. El archivo tendra como output principal el archivo $\texttt{subtypeLUAD.tsv}$

2. Ejecutar *1_2Prepo-mRNA.R*. La principal función de este archivo es procesar los datos mRNA. El archivo tendra como output pricipal el archivo $\texttt{RNAseqnormalized.tsv}$. 
**Warning** Verificar que la matriz resultante ($\texttt{RNAseqnormalized}$) no tenga filas de ceros.

3. Ejecutar *1_3prepo-miRNA.R*. La principal función de este archivo es procesar los datos de miRNA. El archivo tendrá como output principal el archivo $\texttt{miRNAseqNormi.tsv}$.

4. Ejecutar *1_4prepoMethy.R*. La principal función de este archivo es procesar los datos de Methylation. El archivo tendrá como output principal el archivo $\texttt{MethyM.tsv}$.

5. Ejecutar *1_5CONCAT_FINAL.R*. La principal función de este archivo es concatenar la información anteriormente obtenida y separala por subtipos. Los outputs principales son $\texttt{normal.MTRX}$, $\texttt{prox.-inflam.MTRX}$, $\texttt{prox.-prolif..MTRX}$, $\texttt{TRU.MTRX}$.

6. Ejecutar *1_6mfanormi.R*, esto dara como output el archivo $\texttt{LUAD.eigenNormi}$

**Warning** Las librerias $\texttt{NOISeq}$ y $\texttt{data.table}$ tienen ciertos conflictos con las funciones $\texttt{dat}, \texttt{ReadData}$ por lo cual se recomienda verificar que antes de cada una de estas funciones se encuentre ya sea $\texttt{NOISeq::}$ o $\texttt{data.table::}$ antes de ejecutar la respectiva función.
