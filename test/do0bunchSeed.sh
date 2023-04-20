#!/bin/bash

. ./conf.bash

EXE="${SITRACK_DIR}/tools/generate_idealized_seeding.py"


#LIST_DATES="19970101 19970104 19970107 19970110"

LIST_DATES="19970107"

for date in ${LIST_DATES}; do

    YYYY=`echo ${date} | cut -c1-4`
    MM=`echo ${date} | cut -c5-6`
    DD=`echo ${date} | cut -c7-8`
    export NDATE0="${YYYY}-${MM}-${DD}"
    export LDATE0="${NDATE0}_00:00:00"


    echo $NDATE0
    echo

    fout="./nc/sitrack_seeding_nemoTsi3_${NDATE0}_HSS${iHSS}.nc"


    CMD="${EXE} -d ${LDATE0} -m ${FNMM} -S ${iHSS} -i ${FSI3IN} -k 0 -f ${FFSM} -C ${ICOARSEN}"
    echo
    echo " *** About to launch:"; echo "     ${CMD}"; echo
    ${CMD}

done
