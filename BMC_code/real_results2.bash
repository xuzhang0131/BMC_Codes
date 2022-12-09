#!/bin/bash

cd /home/zzh294/BMC_Codes/BMC_code

./myMCMC 2001 5000 &

./myMCMC 5001 8000 &

./myMCMC 8001 11000 &

./myMCMC 11001 13500 &

wait
