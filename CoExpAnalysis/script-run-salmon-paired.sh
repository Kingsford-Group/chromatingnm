#!/bin/bash
#salmon_exec="/mnt/scratch13/salmon-0.6.0/bin/salmon"
salmon_index="/mnt/disk17/RNAtoTADs/salmon_human_index"
#"/mnt/scratch13/salmon-0.6.0/bin/cdna_index"

#${salmon_exec} quant -p 16 -i ${salmon_index} -l ISR -1 <(gunzip -c $1) -2 <(gunzip -c $2) -o $3
#${salmon_exec} quant -p 16 -i ${salmon_index} -l ISR -1 <(gunzip -c $1) -2 <(gunzip -c $2) -o $3 --useVBOpt --numBootstraps 30

time salmon quant -p 16 -i ${salmon_index} -l IU -1 <(gunzip -c $1) -2 <(gunzip -c $2) -o $3
