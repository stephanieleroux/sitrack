#!/bin/bash

YEAR="1997"
#YEAR="2016"

SITRACK_DIR="${HOME}/DEV/sitrack"

NEMO_CONF="NANUK4"

cxtra=""

export DATE1="${YEAR}0101"

# /data/gcm_setup/NANUK4/BBM00$ NANUK4_ICE-BBM00_1h_19970101_19970331_icemod.nc4
NEMO_EXP="BBM00"  ; SBDIR="TEST" ; export DATE2="${YEAR}0331"
#NEMO_EXP="BBM2300" ; SBDIR="00000001-00002976"; cxtra="crndg2nc1"
#NEMO_EXP="BBM2300" ; SBDIR="00000001-00002976" ; cxtra="nc1"
#NEMO_EXP="BBM2300" ; SBDIR="00000001-00002976" ; cxtra="nc0p5"
#NEMO_EXP="BBM2300" ; SBDIR="00000001-00002976" ; cxtra="ref23"
#export DATE2="${YEAR}0131"

export iHSS=4 ; RESKM=12

NJPAR=4 ; # number of jobs we can launch in //

FREQ_AN_DAYS=3 ; # frequency in days of the deformation analysis...

xtra_sfx=""

host=`hostname | cut -d '.' -f2`
case ${host} in
    "merlat")
        #/MEDIA/data/NANUK4/BBM00/NANUK4_ICE-BBM00_1h_19970101_19970331_icemod_LIGHT480.nc4
        export DATA_DIR="/MEDIA/data"
        #export iHSS=6 ; RESKM=73
        #export iHSS=4 ; RESKM=49
        export iHSS=10 ; RESKM=120
        Ndays=6
        xtra_sfx="_LIGHT480"
        #
        ;;
    "mcp-oceannext-01")
        export DATA_DIR="/data/gcm_setup"
        #
        iHSS=8
        NJPAR=4       
        #
        ;;
    "frazilo")
        #export DATA_DIR="/data"
        export DATA_DIR="/home/data/laurent/tmp"
        Ndays=31
        #
        NJPAR=30
        #NJPAR=1
        #
        #SI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${YEAR}0101_${YEAR}0205_icemod.nc4" ; # 1 month !!!
        #
        ;;
    *)
        echo "Unsupported host: ${host} !"
        exit
esac


FSI3IN="${NEMO_CONF}_ICE-${NEMO_EXP}_1h_${DATE1}_${DATE2}_icemod${xtra_sfx}.nc4"

export FSI3IN="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}_ICE-${NEMO_EXP}-S/${SBDIR}${cxtra}/${FSI3IN}"

export FNMM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L31-I/mesh_mask_${NEMO_CONF}_L31_4.2_1stLev.nc"

export FFSM="${DATA_DIR}/${NEMO_CONF}/${NEMO_CONF}.L31-I/mask_RGPS_${NEMO_CONF}.nc"

mkdir -p ./figs ./npz

YYYY=`echo ${DATE1} | cut -c1-4`
MM=`echo ${DATE1} | cut -c5-6`
DD=`echo ${DATE1} | cut -c7-8`
export NDATE1="${YYYY}-${MM}-${DD}"
export LDATE1="${NDATE1}_00:00:00"

