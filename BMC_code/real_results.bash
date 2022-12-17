#!/bin/bash

cd /home/zzh294/BMC_Codes/BMC_code

./myMCMC 2501 5000 > flux_7501_10000.txt &

./myMCMC 5001 8000 > flux_10001_13000.txt &

wait
