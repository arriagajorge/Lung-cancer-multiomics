for i in $(ls *filtered.alt)
do
	#nombre_archivo=$(basename "$i" .filtered.alt)
        Rscript 4_4plot_graph.R $i
	#ext_cys="$nombre_archivo.cys"
	#Rscript 4_5.R $ext_cys
done

