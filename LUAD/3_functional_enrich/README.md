# Functional enrichment.

Para replicar el análists ejecutar en el siguiente orden.

1. Ejecutar el archivo *3_1Selected_features.R*. Selecciona las feautures más estables. El output principal es **LUAD.stable**
2. Ejecutar el archivo *3_2functions_overrepre.R*. Enriquecimiento funcional (sobrerrepresentación). El output principal son los archivos **BP.enrichment, KEGG.enrichment**.
3. Ejecutar el archivo *3_3function_over.R*. Comprueba si las categorías de procesos biológicos están sobrerrepresentadas.
4. Ejecutar el archivo *3_4*. Group the functions that are over-represented in the same SGCCA components. El output principal son los archivos **Groups_per_component.tsv**  y **groups.RDS**.
