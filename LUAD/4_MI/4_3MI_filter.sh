for i in $(ls *.sort)
do
        Rscript 4_3MI_filter.R $i
done

