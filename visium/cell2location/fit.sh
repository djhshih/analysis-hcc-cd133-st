#!/bin/bash
# Fit models

set -euo pipefail
IFS=$'\n\t'

samples=(DTA18 DTA24)

#mkdir -p log

for x in ${samples[@]}; do
	./fit-model.py $x
done

