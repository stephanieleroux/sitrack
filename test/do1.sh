#!/bin/bash

. ./conf.bash

EXE1="${SITRACK_DIR}/tools/generate_idealized_seeding.py"
EXE2="python3 -u ${SITRACK_DIR}/si3_part_tracker.py"

echo $NDATE1
echo

#ICOARSEN

fout="./nc/sitrack_seeding_nemoTsi3_${NDATE1}_HSS${iHSS}.nc"

if [ ! -f ${fout} ]; then
    #${EXE1} -d ${LDATE1} -m ${FNMM} -S ${iHSS} -i ${FSI3IN} -k 0 -f ${FFSM}  ; #-C ${ICOARSEN}
    ${EXE1} -d ${LDATE1} -m ${FNMM}  -i ${FSI3IN} -k 0 -f ${FFSM}  -C 10
fi

exit


# Actually that the ice tracker that should look inside the nc file to get date 1 and 2:
CMD="${EXE2} -i ${FSI3IN} -m ${FNMM} -s ${fout}" ; # with nc file for init seed...
echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
