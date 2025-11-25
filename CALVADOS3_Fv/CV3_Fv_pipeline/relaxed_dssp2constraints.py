'''Script to extract constraints from a DSSP file.'''
import yaml
import argparse
import os
from dssp_parsing import (
    OneLineListDumper,
    get_h_bonded_residues,
    load_map,
    sequential_loop_boundaries,
    fw_boundaries_str,
    get_ss_resi
) 


parser = argparse.ArgumentParser(description='Convert an AB IMGT PDB file into Calvados PDB format')

parser.add_argument(
    '--dssp-file',
    type=str,
    help='path to the input DSSP file',
    required=True
)

parser.add_argument(
    "--constraints-dir",
    type=str,
    help="path to the dir to place the constraints file",
    required=True
)

parser.add_argument(
    "--map-file",
    type=str,
    help="path to file with imgt --> sequential mapping",
    required=True
)

parser.add_argument(
    "--chains",
    type=str,
    help="Ligh and heavy chain labels, L need to be first",
    required=True
)

parser.add_argument(
    "--energy-threshold",
    type=float,
    help="Threshold to use for H-bonding energy",
    required=False,
    default=-0.75
)


def resi_to_constraints_file(h_bonded_resi, *, file_name, chains, mAb_name, imgt2sequential):

    constraints = []
    for chain in chains:
        if chain != '-': # ignore loop for - chain in nanobodies
            mp = imgt2sequential[chain]
            chain_start = min(mp.values())
            constraints.append([chain_start, mp[fw_boundaries_str[0][0]]])
        
            for i in range(3):
                # h bonded resi in cdr
                constraints += [[int(resi), int(resi)] for resi in h_bonded_resi[chain][i]]
                # next fw
                constraints.append([mp[fw_boundaries_str[i][1]], mp[fw_boundaries_str[i+1][0]]] if i < 2 else 
                                [mp[fw_boundaries_str[i][1]], max(mp.values())])

    with open(file_name, 'w') as yaml_file:
        yaml.dump({mAb_name: [constraints]}, yaml_file, Dumper=OneLineListDumper, default_flow_style=None)

def filter_dictionary(data_dict, valid_values):
    # Convert valid_values to strings for comparisonrelaxed_dssp2constraints
    valid_values_str = set(str(v) + ' ' for v in valid_values)

    # Filter the dictionary
    filtered_dict = {
        key: [
            [val for val in sublist if val in valid_values_str]  # Filter sublist
            for sublist in value
        ]
        for key, value in data_dict.items()
    }

    # Remove empty sublists
    filtered_dict = {key: [sublist for sublist in value if sublist] for key, value in filtered_dict.items()}

    return filtered_dict
def filter_dictionary_with_empty_sublists(data_dict, valid_values):
    # Convert valid_values to strings for comparison
    valid_values_str = set(str(v) + ' ' for v in valid_values)

    # Filter the dictionary, retaining empty sublists
    filtered_dict = {
        key: [
            [val for val in sublist if val in valid_values_str] or []  # Retain empty sublist
            for sublist in value
        ]
        for key, value in data_dict.items()
    }

    return filtered_dict

def filter_third_sublists(data_dict, valid_values):
    # Convert valid_values to strings for comparison
    valid_values_str = set(str(v) + ' ' for v in valid_values)

    # Filter only the third sublist in each key
    filtered_dict = {}
    for key, value in data_dict.items():
        # Check if the third sublist exists
        if len(value) >= 3:
            # Apply filtering only to the third sublist
            value[2] = [val for val in value[2] if val in valid_values_str] or []
        # Retain the structure of the dictionary
        filtered_dict[key] = value

    return filtered_dict

def main(dssp_file, constraints_dir, chains, energy_threshold, map_file):
    abname = dssp_file.split('/')[-1].split('_')[0] + '_Fv_' +chains
    constraints_file = os.path.join(constraints_dir, f'{abname}_constraints.yaml')

    imgt2sequential = load_map(map_file)
    
    seq_loop_boundaries = sequential_loop_boundaries(imgt2sequential, chains)

    constrained_HB_resi={}
    for chain in chains:
        if chain != '-': # ignore loop for - chain in nanobodies
            constrained_HB_resi[chain] = get_h_bonded_residues(
                dssp_file, seq_loop_boundaries[chain], chains=chain, energy_threshold=energy_threshold
                )[chain]
        else:
            constrained_HB_resi[chain] = []

    ss_res= get_ss_resi(dssp_file)
    constrained_resi=filter_dictionary_with_empty_sublists(constrained_HB_resi,ss_res)

    resi_to_constraints_file(
        constrained_resi,
        file_name=constraints_file,
        chains=chains,
        mAb_name=abname,
        imgt2sequential=imgt2sequential
    )

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.dssp_file, args.constraints_dir, args.chains, args.energy_threshold, args.map_file)
