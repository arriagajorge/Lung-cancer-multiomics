# SGCCA.

To replicate the analysis run in the following order.

1. Execute the *2_1FitComm.R* file. The file adjusts penalty values for each omics. The main outputs of this file are **0.01.0.01.0.01.tsv ... 0.99.0.99.99.0.99.tsv**. 
2. Execute the *2_2ChoosePen.R* file. Select the penalty values with higher AVE and lower nfeatures.
3. Execute file *2_3SGCCA.R*. Run SGCCA, using the adjusted penalty values. The main output of this file is **LUAD.selected**.
4. Run the file *2_4SGCCA_subsamplesAux.R*. SGCCA with 100 subsets of data to test the stability of the features. The main outpus are **LUAD.1.selected ... LUAD.100.selected**.
5. Run the *2_5JoinSubsamples.R* file. Count the selected features in the subsamples to check stability. The main output of this file is **LUAD.subsampled**.


