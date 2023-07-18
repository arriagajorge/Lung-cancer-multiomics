# MI Networks.

Para replicar el análists ejecutar en el siguiente orden.

1. Ejecutar el archivo *4_1get_submatrixAux.R*. Obteneción de una matriz con todas las características en los componentes SGCCA donde
la función está sobrerrepresentada. Se obtienen archivos de la forma *funcion.LUAD.mtrx*.
2. Correr ARACNE con la matriz obtenida. Ocupamos el repositorio https://github.com/CSB-IG/ARACNE-multicore.
   1. Copiar las matrices obtenidas al directorio launch.
   2. Tambien copiar los archivos **rename.sh** y **mult-tsvs.sh** al directorio launch.
   3. Ejecutar el archivo **mult-tsvs.sh** con el comando:
        ```
        bash mult-tsvs.sh
        ```
    Se obtienen archivos de la forma *funcion.LUAD.sort*.

3. Ejecutar el comando:
    ```
    bash 4_3MI_filter.sh
    ```
    Mantiene sólo las aristas con MI superior a un umbral específico para el tipo de arista. Se obtienen archivos de la forma *funcion.LUAD.filtered.alt*.
  
4. Ejecutar el comando:
    ```
    bash 4_4plot_graph.sh
    ```
    Gráfica las interacciones devueltas por *4_3MIfilter.sh*. Se obtienen archivos de la forma *funcion.LUAD.filtered.cys*.
5. Ejecutar el archivo *4_5Networks.R*. Obteniene información sobre los nodos y aristas gráficados.
