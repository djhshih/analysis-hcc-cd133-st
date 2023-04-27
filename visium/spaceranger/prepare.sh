#!/bin/bash
# Set up symlink

# scanpy expects tissue_positions_list.csv but
# we have tissue_positions.csv instead
#
set -euo pipefail
IFS=$'\n\t'

for s in */; do
	echo $s
	(
		cd $s/outs/spatial
		if [[ ! -f tissue_positions_list.csv ]]; then
			ln -s tissue_positions.csv tissue_positions_list.csv
		fi
	)
done

