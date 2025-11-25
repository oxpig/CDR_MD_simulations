import os
from dssp_parsing import *
import argparse

parser = argparse.ArgumentParser(description='Get indices of hydrogen bond donors and acceptors and energy')

parser.add_argument(
    '--dssp-file',
    type=str,
    help='path to the input DSSP file',
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
    "--energy-threshold",
    type=float,
    help="Minimum threshold to use for H-bonding energy",
    required=False,
    default=-0.41
)


def add_dssp_line_idx(protein_df):
    first_chain = protein_df.iloc[0].chain#
    n_resi_first_chain = protein_df.chain.value_counts()[first_chain].item()


    second_chain = protein_df.iloc[-1].chain
    if first_chain != second_chain: # ignore second chain for nanobodies
        n_resi_second_chain = protein_df.chain.value_counts()[second_chain].item()
    else:
        n_resi_second_chain = 0

    dssp_line_idx = np.concatenate([
        np.arange(0,n_resi_first_chain,1),
        np.arange(n_resi_first_chain, n_resi_first_chain+n_resi_second_chain,1)+1
    ])

    protein_df['dssp_line_idx'] = dssp_line_idx

    return protein_df

def check_j(j, j_relative, energy):
    correct = False
    if j['N_h_bond_1_energy'] == energy and j['N_h_bond_1'] == -j_relative:
        correct = True
    elif j['O_h_bond_1_energy'] == energy and j['O_h_bond_1'] == -j_relative:
        correct = True
    elif j['N_h_bond_2_energy'] == energy and j['N_h_bond_2'] == -j_relative:
        correct = True
    elif j['O_h_bond_2_energy'] == energy and j['O_h_bond_2'] == -j_relative:
        correct = True

    assert correct, "H-bond acceptor doesn't match"

def update_j_numbering(protein_df, i_jrel_energy, chain_label):
    ij_energy = []
    for i, j_relative, energy in i_jrel_energy:
        j_index = protein_df.query(
            f'chain == "{chain_label}" and residue == {i}'
            ).dssp_line_idx.item() + j_relative
        j = protein_df[protein_df.dssp_line_idx == j_index].squeeze()
        try:
            check_j(j, j_relative, energy)
        except AssertionError:
            # print(f"Skipping {i} {j_relative} {energy}")
            continue
        except ValueError:
            # print(f"Skipping {i} {j_relative} {energy}, no infomation found for j")
            continue
        ij_energy.append((i, int(j['residue']), energy))
    
    return ij_energy

def h_bond_pairings_and_energy(protein_df, chain_label, loop_boundary, energy_threshold):
    loop_region = protein_df[(protein_df['chain'] == chain_label) &
                             (protein_df['residue'] >= loop_boundary[0]) &
                             (protein_df['residue'] <= loop_boundary[1])]
    
    i_jrel_energy = h_bonds_and_energy(loop_region, energy_threshold)
    ij_energy = update_j_numbering(protein_df, i_jrel_energy, chain_label)

    return ij_energy

def main(dssp_file, map_file, chains, output_dir, energy_threshold):
    abname = dssp_file.split('/')[-1].split('_')[0] + '_Fv_' +chains
    out_file = os.path.join(output_dir, f'{abname}_ij_energy.csv')

    imgt2sequential = load_map(map_file)
    seq_loop_boundaries = sequential_loop_boundaries(imgt2sequential, chains)

    with open(dssp_file) as f:
        dssp = f.readlines()

    # extract region around CDR3 loop to inspect for loops in DSSP file
    header = dssp_header()
    protein_lines = dssp[dssp.index(header)+1:]
    protein_df = dssp_format_to_df(protein_lines)
    protein_df = add_dssp_line_idx(protein_df)

    ij_energy = []
    for chain in chains:
        if chain != '-': # ignore loop for - chain in nanobodies
            for loop_idx in range(3):
                loop_boundary = seq_loop_boundaries[chain][loop_idx]
                ij_energy += h_bond_pairings_and_energy(protein_df, chain, loop_boundary, energy_threshold)
        else:
            continue

    pd.DataFrame(ij_energy, columns=['i', 'j', 'parameter']).to_csv(out_file, index=False)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.dssp_file, args.map_file, args.chains, args.out_dir, args.energy_threshold)
