#!/bin/bash

# Deformation only
. ./conf.bash

EXE="${SITRACK_DIR}/deformation.py"

mkdir -p logs

ijob=0

# Populating the batches available:
listQ=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_nemoTsi3_${YEAR}????t0_${YEAR}????_*km.npz`

echo "${listQ}"

echo; echo
echo " *** ${RESKM} km ***"
echo

list=`\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_nemoTsi3_${YEAR}????t0_${YEAR}????_${RESKM}km.npz`
nbf=`echo ${list} | wc -w`

echo " *** Number of files = ${nbf}"

list_date_ref=""
for ff in ${list}; do
    date_ref=`echo ${ff} | cut -d_ -f6`
    list_date_ref+=" ${date_ref}"
done
list_date_ref=$(echo ${list_date_ref} | tr ' ' '\n' | sort -nu) ; # unique and sorted !

echo; echo " *** List of reference dates:"; echo "${list_date_ref}"; echo

for dr in ${list_date_ref}; do
    echo
    lst=( `\ls npz/Q-mesh_NEMO-SI3_${NEMO_CONF}_${NEMO_EXP}_nemoTsi3_${dr}_${YEAR}????_${RESKM}km.npz` )
    nf=`echo ${lst[*]} | wc -w` ; #echo " => ${nf} files "

    if [ ${nf} -eq 2 ]; then

        fQ1=${lst[0]}
        fQ2=${lst[1]}
        echo " ==> will use:"
        echo " * ${fQ1}"
        echo " * ${fQ2}"

        flog="quadgen_${NEMO_CONF}-${NEMO_EXP}_${RESKM}_${dr}"

        ijob=$((ijob+1))

        CMD="${EXE} ${fQ1} ${fQ2} 0"
        echo "  ==> ${CMD}"; echo
        ${CMD} 1>logs/out_${flog}.out 2>logs/err_${flog}.err &
        echo; echo

        if [ $((ijob%NJPAR)) -eq 0 ]; then
            echo "Waiting! (ijob = ${ijob})...."
            wait
            echo; echo
        fi

    fi

done

wait
