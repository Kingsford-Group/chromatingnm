#!/bin/bash
salmon_index=/mnt/disk17/RNAtoTADs/salmon_human_index

time salmon quant -p 16 -i ${salmon_index} -r <(gunzip -c $1) -l IU -o $2
