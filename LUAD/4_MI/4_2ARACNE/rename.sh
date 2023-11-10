#!/bin/bash

# Obtener el nombre del archivo como argumento
nombre=$1

# Verificar si se proporcionó un nombre de archivo
if [ -z "$nombre" ]; then
  echo "No se ha proporcionado un nombre de archivo."
  exit 1
fi

# Guardar el nombre original del archivo
nombre_original="$nombre"

#echo $nombre_original

# Cambiar el nombre del archivo a GO_temp
mv "$nombre" "temp.mtrx"
#echo "El archivo ha sido renombrado a temp.mtrx"

bash run.sh temp.mtrx &> "${nombre_original%.mtrx}_salida.txt" &
echo "ARACNE has been run"

# esperar a que termine ARACNE
wait

#rename el .sif
mv temp.sif "${nombre_original%.mtrx}.sif"

# rename el .sort
mv temp.sort "${nombre_original%.mtrx}.sort"

# rename the matrix temp
mv temp.mtrx "$nombre_original"

wait
