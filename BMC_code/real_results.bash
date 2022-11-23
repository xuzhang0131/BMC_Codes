#!/bin/bash

cd /home/zzh294/BMC_Codes/BMC_code

./myMCMC 101 500 &

./myMCMC 501 1000 &

./myMCMC 1001 1500 &

./myMCMC 1501 2000 &

wait
