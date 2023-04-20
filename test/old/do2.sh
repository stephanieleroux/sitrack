#!/bin/bash

. ./conf.bash

EXE="python3 -u ${SITRACK_DIR}/generate_quad_mesh.py"

fin=`\ls ./nc/NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_tracking_nemoTsi3_${DATE1}h00_${YEAR}????h??.nc`


if [ `echo ${fin} | wc -w ` -ne 1 ]; then
    echo "ERROR: problem with available NC files!!!"
    exit
fi


Nb3ds=$((Ndays/FREQ_AN_DAYS-1))

istep=$((FREQ_AN_DAYS*24))

echo "Nb3ds=$Nb3ds, istep=$istep"

ijob=0
mkdir -p logs

jk=1
ik2=0
while [ ${jk} -lt ${Nb3ds} ]; do
    ik1=${ik2}
    ik2=$((ik1+istep))

    flog="quadgen_${NEMO_CONF}-${NEMO_EXP}_${RESKM}_${ik1}_${ik2}"
    
    CMD="${EXE} ${fin} ${ik1},${ik2} ${RESKM}"
    echo "    ==> will launch:"; echo "     ${CMD}"; echo
    ${CMD} 1>"./logs/out_${flog}.out" 2>"./logs/err_${flog}.err" &
    sleep 1
    echo

    ijob=$((ijob+1))

    if [ $((ijob%NJPAR)) -eq 0 ]; then
        echo "Waiting! (ijob = ${ijob})...."
        wait
        echo; echo
    fi

    jk=$((jk+1))
done

wait
