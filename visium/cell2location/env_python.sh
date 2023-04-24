#!/bin/bash
# 

conda create -y -n st python=3.9
conda activate st

conda install -y -c conda-forge \
	numpy matplotlib pandas \
	scanpy python-igraph leidenalg scvi-tools stlearn mudata squidpy \

pip install git+https://github.com/BayraktarLab/cell2location.git#egg=cell2location[tutorials]
pip install tangram-sc

