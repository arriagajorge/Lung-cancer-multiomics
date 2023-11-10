# Functional enrichment.

To replicate the analysis run in the following order.

1. Run the *3_1Selected_features.R* file. Select the most stable feautures. The main output is **LUAD.stable**.
2. Execute file *3_2functions_overrepre.R*. Functional enrichment (overrepresentation). The main output are the files **BP.enrichment, KEGG.enrichment**.
3. Execute the *3_3function_over.R* file. Check if the biological process categories are overrepresented.
4. Run file *3_4*. Group the functions that are over-represented in the same SGCCA components. The main output are the files **Groups_per_component.tsv** and **groups.RDS**.
