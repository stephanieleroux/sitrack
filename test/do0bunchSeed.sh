#!/bin/bash

. ./conf.bash

EXE="${SITRACK_DIR}/tools/generate_idealized_seeding.py"

EXE2="${SITRACK_DIR}/tools/listofdates.py"

LIST_DATES=`${EXE2} ${DATE1} ${DATE2} 3`

echo
echo "LIST_DATES = ${LIST_DATES}"
echo


mkdir -p ./logs

ijob=0


for date in ${LIST_DATES}; do

    YYYY=`echo ${date} | cut -c1-4`
    MM=`echo ${date} | cut -c5-6`
    DD=`echo ${date} | cut -c7-8`
    export NDATE0="${YYYY}-${MM}-${DD}"
    export LDATE0="${NDATE0}_00:00:00"


    echo $NDATE0
    echo

    for ICOARSEN in ${LCOARSEN}; do

        str="${NDATE0}_crsn{ICOARSEN}km"

        flog="seeding_${str}"

        fout="./nc/sitrack_seeding_nemoTsi3_${NDATE0}_${ICOARSEN}km.nc"

        if [ ! -f ${fout} ]; then

            CMD="${EXE} -d ${LDATE0} -m ${FNMM} -i ${FSI3IN} -k 0 -f ${FFSM} -C ${ICOARSEN}"
            echo
            echo " *** About to launch:"; echo "     ${CMD}"; echo
            ijob=$((ijob+1))
            ${CMD} 1>"./logs/out_${flog}.out" 2>"./logs/err_${flog}.err" &
            echo
            if [ $((ijob%NJPAR)) -eq 0 ]; then
                echo "Waiting! (ijob = ${ijob})...."
                wait; echo; echo
            fi

        fi
    done
done
