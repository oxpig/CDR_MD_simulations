#!/usr/bin/env bash

### script to install environments for simulation


## prep environment
echo "[INFO] Setting up prep_env"
# For prep environment install openmm 8.2.0 (in yaml) from pip followed by openmm 7.7.0
# This is a hack for a current issue and will be fixed in future releases
conda env create -f prep_env.yaml
conda activate prep_env
conda install openmm=7.7.0
conda deactivate


## calvados environment
echo "[INFO] Setting up calvados environment"
conda create -n calvados python=3.10 -y
conda activate calvados
cd CALVADOS3_Fv
pip install .
cd ..
conda deactivate


## cg2all environment
echo "[INFO] Setting up cg2all environment"
# This is an example with cudatoolkit=11.3.
# Set a proper cudatoolkit version that is compatible with your CUDA driver and DGL library.
# dgl>=1.1 occasionally raises some errors, so please use dgl<=1.0.
conda create --name cg2all pip cudatoolkit=12.8 dgl=1.0 -c dglteam/label/cu128 python=3.10 -y
conda activate cg2all
pip install git+http://github.com/huhlim/cg2all
conda deactivate

echo "Done"