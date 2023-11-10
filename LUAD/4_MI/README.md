# MI Networks.

To replicate the analysis run in the following order.

1. Execute the *4_1get_submatrixAux.R* file. Obtain a matrix with all the features in the SGCCA components where the function is overrepresented.
the function is overrepresented. Files of the form *function.LUAD.mtrx* are obtained.
2. Run ARACNE with the obtained matrix. We use the repository https://github.com/CSB-IG/ARACNE-multicore.
   1. Copy the obtained matrices to the launch directory.
   2. Also copy the files **rename.sh** and **mult-tsvs.sh** to the launch directory.
   3. Execute the **mult-tsvs.sh** file with the command:
        ```
        bash mult-tsvs.sh
        ```
    You get files of the form *function.LUAD.sort*.

3. Execute the command:
    ```
    bash 4_3MI_filter.sh
    ```
    Keeps only edges with MI greater than a specific threshold for the edge type. You get files of the form *function.LUAD.filtered.alt*.
  
4. Execute the command:
    ```
    bash 4_4plot_graph.sh
    ```
    Graph the interactions returned by *4_3MIfilter.sh*. You get files of the form *function.LUAD.filtered.cys*.
5. Execute the *4_5Networks.R* file. Obtain information about the graphed nodes and edges.

