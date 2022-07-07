#!/usr/bin/env sh

cd ..

DATE_FROM=22050601
DATE_TO=22050608
FILENAME=relu_yeswane_sev2.0_22021104

Rscript fit-omi.R relu all \
    ${DATE_FROM} \
    ${DATE_TO} \
    yeswane \
    central \
    2.0 1.0 seaslate 0.1 ${FILENAME} 2.5 0.551 0.5 0.5 1.5 0.3
