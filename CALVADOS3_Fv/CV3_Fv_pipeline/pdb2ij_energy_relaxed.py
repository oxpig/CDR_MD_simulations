import os
from dssp_parsing import *
import argparse
import numpy as np
import MDAnalysis as mda

parser = argparse.ArgumentParser(description='Get indices of hydrogen bond donors and acceptors and energy')

parser.add_argument(
    '--pdb-file',
    type=str,
    help='path to the input PDB file',
    required=True
)

parser.add_argument(
    '--energy-HB-file',
    type=str,
    help='path to the input HB energy file',
    required=True
)

parser.add_argument(
    "--out-dir",
    type=str,
    help="path to the dir to place output",
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
    "--dist-threshold",
    type=float,
    help="Minimum threshold to use for H-bonding energy",
    required=False,
    default=5
)

parser.add_argument(
    '--flag-value',
    type=float,
    help='value to assign to the contacts',
    default=1
)   

def sidechain_contacts(u, selection="protein", cutoff=4.5, frame=-1, filter_resids=None):
    """
    Calculate sidechain contacts between residues in a protein using the center of mass of sidechains,
    except for Glycine where the CA position is used.
    
    Parameters:
    - u: MDAnalysis Universe object
    - selection: Atom selection string (default: "protein")
    - cutoff: Distance threshold in Angstroms (default: 4.5 Å)
    - frame: Frame to analyze (-1 for last frame)
    - filter_resids: List of residue IDs to filter results (default: None, returns all contacts)
    
    Returns:
    - contacts: List of residue ID pairs in contact
    """
    # Compute the center of mass for each sidechain (or CA for Glycine)
    residues = []
    positions = []
    
    for res in u.select_atoms(selection).residues:
        if res.resname == "GLY":
            pos = res.atoms.select_atoms("name CA").positions[0]
        else:
            pos = res.atoms.select_atoms("not backbone and not name H*").center_of_mass()
        residues.append(res.resid)
        positions.append(pos)
    
    positions = np.array(positions)
    
    # Compute pairwise distances
    u.trajectory[frame]  # Move to the requested frame
    dist_matrix = np.linalg.norm(positions[:, np.newaxis] - positions, axis=-1)
    # Identify contacts based on cutoff
    contact_indices = np.where((dist_matrix < cutoff) & (dist_matrix > 0))  # Exclude self-contacts
    contacts = set()
    for i, j in zip(*contact_indices):
        if abs(residues[i] - residues[j]) >= 2:  # Skip contacts where |j - i| < 2
            contact_pair = tuple(sorted((residues[i], residues[j])))  # Ensure (i, j) with i < j
            if filter_resids is None or contact_pair[0] in filter_resids or contact_pair[1] in filter_resids:
                contacts.add(contact_pair)
                # print(contact_pair)
    return sorted(contacts)

def remove_duplicate_contacts(df):
    """
    Remove duplicate (i, j) pairs in the DataFrame, keeping the row with the highest 'parameters' value.
    
    Parameters:
    - df: Pandas DataFrame with columns ["i", "j", "parameters"]
    
    Returns:
    - df_filtered: Filtered DataFrame with duplicates removed.
    """
    df_filtered = df.sort_values(by="parameter", ascending=True).drop_duplicates(subset=["i", "j"], keep="first")
    df_filtered = df_filtered.sort_values(by=["i", "j"]).reset_index(drop=True)
    return df_filtered

def main(pdb_file, map_file, chains, output_dir, dist_threshold, energy_HB_file, flag_value):
    abname = pdb_file.split('/')[-1].split('_')[0] + '_Fv_' +chains
    out_file = os.path.join(output_dir, f'{abname}_ij_energy_contact.csv')

    imgt2sequential = load_map(map_file)
    seq_loop_boundaries = sequential_loop_boundaries(imgt2sequential, chains)

    u = mda.Universe(pdb_file)

    loop_resids = []
    for chain in chains:
        if chain != '-':
            for loop_idx in range(3):
                loop_boundary = seq_loop_boundaries[chain][loop_idx]
                loop_resids = loop_resids + list(range(loop_boundary[0], loop_boundary[1] + 1))
        else:
            continue

    contacts = sidechain_contacts(u, selection="protein", filter_resids=loop_resids,cutoff=dist_threshold)

    parameters = [flag_value for i in range(len(contacts))]

    assert len(contacts) == len(parameters), "Length of contacts and parameters do not match"

    df_in=pd.read_csv(energy_HB_file)

    # Create DataFrame
    df_contacts = pd.DataFrame(contacts, columns=["i", "j"])
    df_contacts["parameter"] = parameters  # Ensure correct length
    
    # check that the dfs have the same columns
    assert df_in.columns.equals(df_contacts.columns), "Columns in input and output DataFrames do not match"

    # Save DataFrame
    df_out = pd.concat([df_in, df_contacts], axis=0)

    df_out=remove_duplicate_contacts(df_out)
    df_out.to_csv(out_file, index=False)
    print(f"Output saved to {out_file}")

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.pdb_file, args.map_file, args.chains, args.out_dir, args.dist_threshold, args.energy_HB_file, args.flag_value)
