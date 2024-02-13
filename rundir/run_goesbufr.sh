#!/bin/bash

# Run modulo para crear prepbufr/agregar observaciones
#------------------------------------------------------
namelist=namelist
filelist=flist.txt


anio=2018

for d in {317..323} 
do

diaj=$d

for h in {0..23} 
do

hora=`printf %2.2i $h`			# Hora de asimilación
fecha=$(date -d "`date +$anio`-01-01 +$(( ${diaj} - 1 ))days" +%Y%m%d)

echo '*** Starting ***'
echo $fecha ' hora: ' $hora

obs_dir='/home/paola.corrales/datosmunin3/DATA/GOES16/ABI-L1b-RadF/'${anio}'/'${diaj}'/'${hora}
mask_dir='/home/paola.corrales/datosmunin3/DATA/GOES16/ABI-L2-ACMF/OR_ABI-L2-ACMF-M3_G16_s'${anio}${diaj}${hora}

files_rad=($(ls -d ${obs_dir}/*))
file_mask=($(ls -d ${mask_dir}*))

echo '*** copying files ***'

for file in ${files_rad[*]}
do
  cp ${file} /home/paola.corrales/goes_16/rundir/
done

cp ${file_mask} /home/paola.corrales/goes_16/rundir/

echo '*** writing filillist ***'

rm ${filelist}
touch ./$filelist

files_all=($(ls *.nc))

for file in ${files_all[*]}
do
  echo ${file} >> ${filelist}
done

echo '*** the magic is happening ***'

./GOES16_nc2bufr.exe

echo '*** cleaning ***'

rm *.nc
mv *.bufr /home/paola.corrales/datosmunin3/DATA_EXP/GOESOBS/abig16.${fecha}.t${hora}z.bufr_d

done

done
