#!/bin/bash

. ./conf.bash

EXE1="${SITRACK_DIR}/tools/generate_idealized_seeding.py"
EXE2="python3 -u ${SITRACK_DIR}/si3_part_tracker.py"

echo $NDATE1
echo

fout="./nc/sitrack_seeding_nemoTsi3_${NDATE1}_HSS${iHSS}.nc"

if [ ! -f ${fout} ]; then
    #${EXE1} ${LDATE1} ${FNMM},${iHSS}
    ${EXE1} ${LDATE1} ${FNMM},${iHSS} ${FSI3IN} 0
fi


# Actually that the ice tracker that should look inside the nc file to get date 1 and 2:
CMD="${EXE2} -i ${FSI3IN} -m ${FNMM} -s ${fout}" ; # with nc file for init seed...
echo
echo " *** About to launch:"; echo "     ${CMD}"; echo
${CMD}
