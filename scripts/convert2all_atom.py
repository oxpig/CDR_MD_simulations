import pandas as pd
import subprocess
import argparse
import os

argparser = argparse.ArgumentParser()
argparser.add_argument('--simulation_dir', type=str)
argparser.add_argument('--pdb_id', type=int)
argparser.add_argument('--h_chain', type=str)
argparser.add_argument('--l_chain', type=str)

def main(path, pdb_id, h_chain, l_chain):

    print(f'Processing {path}: {pdb_id}_{h_chain}{l_chain}')

    print('Converting topology')
    os.makedirs(path, exist_ok=True)

    command = f'python CALVADOS3_Fv/convert_top.py --top {path}/top.pdb --out {path}/updated_top.pdb'
    subprocess.run(command, shell=True)

    print('Converting trajectory')
    command = f'convert_cg2all --pdb {path}/updated_top.pdb --cg ResidueBasedModel --dcd {path}/{pdb_id}_Fv_{h_chain}{l_chain}.dcd --out {path}/{pdb_id}_Fv_{h_chain}{l_chain}_all_atom.xtc'
    subprocess.run(command, shell=True)
    command = f'convert_cg2all --pdb {path}/updated_top.pdb --cg ResidueBasedModel --out {path}/{pdb_id}_Fv_{h_chain}{l_chain}_top.pdb'
    subprocess.run(command, shell=True)

    print('Done')


if __name__ == '__main__':
    args = argparser.parse_args()
    main(args.simulation_dir, args.pdb_id, args.h_chain, args.l_chain)
