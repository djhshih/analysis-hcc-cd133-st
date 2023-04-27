#!/bin/bash
# Set up symlink

# scanpy expects tissue_positions_list.csv (which has header row)
# but we have tissue_positions.csv (which has not have header row) instead
#
set -euo pipefail
IFS=$'\n\t'

for s in */; do
	echo $s
	(
		cd $s/outs/spatial
		rm -f tissue_positions_list.csv
		sed 1d tissue_positions.csv > tissue_positions_list.csv
	)
done

