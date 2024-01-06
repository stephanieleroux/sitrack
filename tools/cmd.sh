#!/bin/bash

DATADIR="/Users/leroux/DATA/SASIP/DATA_SItrack/sitrack_demo"

#python generate_idealized_seeding.py -d '1996-12-15_00:00:00' \
#                              -m ${DATADIR}/NANUK4/mesh_mask_NANUK4_L31_4.2_1stLev.nc \
#                              -i ${DATADIR}/NANUK4/NANUK4-BBM23U06_1h_19961215_19970420_icemod.nc \
#                              -k 0 \
#                              -S 5 \
#                              -f  ${DATADIR}/NANUK4/mask_SeedInit_TrackIce_NANUK4.nc \
#                              -N NANUK4

#python generate_sidfex_seeding.py -d '1996-12-15_00:00:00' \
#-f  ${DATADIR}/NANUK4/mask_SeedInit_TrackIce_NANUK4.nc \
#-m ${DATADIR}/NANUK4/mesh_mask_NANUK4_L31_4.2_1stLev.nc \
# --lsidfex 1  -k 0 \
#             -S 5 \
#             -N NANUK4


python generate_sidfex_seeding.py -d '1996-12-15_00:00:00' --lsidfex 1  -k 0 -S 5 
