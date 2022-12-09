#!/bin/bash

cd /home/zzh294/BMC_Codes/noncancer

./myMCMC 1 2000 &

./myMCMC 2001 4000 &

./myMCMC 4001 6000 &

./myMCMC 6001 8400 &

wait
