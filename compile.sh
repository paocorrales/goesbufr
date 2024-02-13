#!/bin/bash

export GSI_COREDIR=/home/paola.corrales/datosmunin3/comGSIv3.7_EnKFv1.3/build

export NETCDF=/usr/local/

 gfortran -o GOES16_nc2bufr.exe -ffree-line-length-none -fbacktrace -fcheck=mem GOES_GRB_nc2bufr.f90 \
    -L${NETCDF}/lib -lnetcdf -lnetcdff -lm -I${NETCDF}/include  \
    -L${GSI_COREDIR}/lib -lbufr_v10.2.5

