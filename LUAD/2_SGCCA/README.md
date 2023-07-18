# SGCCA.

Para replicar el análists ejecutar en el siguiente orden.

1. Ejecutar el archivo *2_1FitComm.R*. El archivo ajusta valores de penalización por cada ómica.
2. Ejecutar el archivo *2_2ChoosePen.R*. Selecciona los valores de penalización con mayor AVE y menor nfeatures.
3. Ejecutar el archivo *2_3SGCCA.R*. Ejecuta SGCCA, utilizando los valores de penalización ajustados. El output principal de este archivo es **LUAD.selected**.
4. Ejecutar el archivo *2_4SGCCA_subsamplesAux.R*. SGCCA con 100 subconjuntos de datos para comprobar la estabilidad de las features.
5. Ejecutar el archivo *2_5JoinSubsamples.R*. Recuento de las características seleccionadas en los subconjuntos para comprobar la estabilidad.El output principal de este archivo es **LUAD.subsampled**.
