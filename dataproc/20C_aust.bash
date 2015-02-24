#!/bin/bash

for FILENAME in prmsl.????.nc
do
    #cdo monavg $FILENAME monthly_$FILENAME
    ncks -v prmsl -d lat,-60.,0. -d lon,90.,200. $FILENAME aus_$FILENAME
done
