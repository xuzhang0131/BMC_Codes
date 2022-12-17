#!/bin/bash

cd /home/zzh294/BMC_Codes/noncancer

./myMCMC 2501 5000 &

./myMCMC 5001 8000 &


wait
