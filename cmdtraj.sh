#!/bin/bash

DATADIR="/Users/leroux/DATA/SASIP/DATA_SItrack/sitrack_demo"

python si3_part_tracker.py -i ${DATADIR}/NANUK4/NANUK4-BBM23U06_1h_19961215_19970420_icemod.nc \
                    -m ${DATADIR}/NANUK4/mesh_mask_NANUK4_L31_4.2_1stLev.nc \
                    -s ./tools/nc/sitrack_seeding_sidfex_19961215_00_HSS5.nc -e 1996-12-26 -N NANUK4 -F 



