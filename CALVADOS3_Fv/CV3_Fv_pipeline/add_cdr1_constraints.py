import sys
import yaml
from dssp_parsing import OneLineListDumper
from dssp_parsing import dssp_header, dssp_format_to_df, load_map
import argparse

parser = argparse.ArgumentParser(description='Add CDR1 constraints to a constraints file')

parser.add_argument(
    '--dssp-file',
    type=str,
    help='path to the input DSSP file',
    required=True
)

parser.add_argument(
    "--map-file",
    type=str,
    help="path to file with imgt --> sequential mapping",
    required=True
)

parser.add_argument(
    "--constraints-file",
    type=str,
    help="path to the constraints file",
    required=True
)


parser.add_argument(
    "--pdb",
    type=str,
    help="PDB ID",
    required=True
)

parser.add_argument(
    "--chain",
    type=str,
    help="Chain labels",
    required=True
)



MaxASA = { # Tien et al. 2013
    'A': 129.0,
    'R': 274.0,
    'N': 195.0,
    'D': 193.0,
    'C': 167.0,
    'E': 223.0,
    'Q': 225.0,
    'G': 104.0,
    'H': 224.0,
    'I': 197.0,
    'L': 201.0,
    'K': 236.0,
    'M': 224.0,
    'F': 240.0,
    'P': 159.0,
    'S': 155.0,
    'T': 172.0,
    'W': 285.0,
    'Y': 263.0,
    'V': 174.0
}

def get_pinned_down_residue(DSSP_file, map_file, hchain, lchain, CDR, threshold=0.075, return_min=False):
    
    header = dssp_header()
    with open(DSSP_file) as f:
        dssp = f.readlines()
    protein_lines = dssp[dssp.index(header)+1:]
    protein_df = dssp_format_to_df(protein_lines)
    protein_df['RSA'] = protein_df.apply(lambda row: row.solvent_accessibility / MaxASA[row.aa], axis=1)

    mp = load_map(map_file)

    if CDR == 'H1':
        inv_map = {v: k for k, v in mp[hchain].items()}
        cdr_start, cdr_end = mp[hchain]['27.0'], mp[hchain]['38.0']
    elif CDR == 'L1':
        cdr_start, cdr_end = mp[lchain]['27.0'], mp[lchain]['38.0']
        inv_map = {v: k for k, v in mp[lchain].items()}

    CDR1 = protein_df[(protein_df.residue >= cdr_start) & (protein_df.residue <= cdr_end)]
    pinned = CDR1[CDR1.RSA <= threshold]
    pinned = pinned[~pinned.residue.apply(lambda x: inv_map[x] in ['27.0', '36.0', '37.0', '38.0'])]
    
    if not return_min: # return all pinned residues
        if pinned.shape[0] == 0:
            return None, None, None
        return [x for x in pinned.residue], [x for x in pinned.aa], [x for x in pinned.RSA]
    else: # return only the lowest RSA
        if pinned.shape[0] == 0:
            return None, None, None, None
        i = pinned.RSA.values.argmin()
        return pinned.iloc[i].residue.item(), pinned.iloc[i].aa, pinned.iloc[i].RSA.item()

def post_process(resi, aas, rsa):
    selection_order = ['F', 'L', 'I', 'V']

    if resi == None:
        return None, None, None
    elif len(resi) == 1:
        return resi[0], aas[0], rsa[0]
    else: # select hydrophobic aa if present
        selected_aa, selected_resi, selected_RSA = [], [], []
        for j, aa in enumerate(aas):
            if aa in selection_order:
                selected_aa.append(aa)
                selected_resi.append(resi[j])
                selected_RSA.append(rsa[j])
        if len(selected_RSA) > 0:
            i = selected_RSA.index(min(selected_RSA))
            return selected_resi[i], selected_aa[i], selected_RSA[i]
        else:

            i = rsa.index(min(rsa)) # select lowest RSA
            return resi[i], aas[i], rsa[i]
        
def main(dssp_file, map_file, constraints_file, pdb, chain):
    h_chain = chain[0]
    l_chain = chain[1]

    # CDRH1
    resi_H1, aa_H1, rsa_H1 = get_pinned_down_residue(dssp_file, map_file, h_chain, l_chain, 'H1', threshold=0.075)
    if resi_H1 == None:
        resi_H1, aa_H1, rsa_H1 = get_pinned_down_residue(dssp_file, map_file, h_chain, l_chain, 'H1', threshold=0.2)
    resi_H1, aa_H1, rsa_H1 = post_process(resi_H1, aa_H1, rsa_H1)

    # CDRL1
    if l_chain != '-': # ignore loop for - chain in nanobodies
        resi_L1, aa_L1, rsa_L1 = get_pinned_down_residue(dssp_file, map_file, h_chain, l_chain, 'L1', threshold=0.075)
        if resi_L1 == None:
            resi_L1, aa_L1, rsa_L1 = get_pinned_down_residue(dssp_file, map_file, h_chain, l_chain, 'L1', threshold=0.2)
        resi_L1, aa_L1, rsa_L1 = post_process(resi_L1, aa_L1, rsa_L1)
    else:
        resi_L1, aa_L1, rsa_L1 = None, None, None
    

    # Update constraints file
    if resi_H1 != None or resi_L1 != None:
        with open(constraints_file) as f:
            constraints = yaml.safe_load(f)
        
        if resi_H1 != None:
            constraints[f'{pdb}_Fv_{h_chain}{l_chain}'][0].append([resi_H1, resi_H1])
        if resi_L1 != None:
            constraints[f'{pdb}_Fv_{h_chain}{l_chain}'][0].append([resi_L1, resi_L1])

        with open(constraints_file, 'w') as yaml_file:
                yaml.dump(constraints, yaml_file, Dumper=OneLineListDumper, default_flow_style=None)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.dssp_file, args.map_file, args.constraints_file, args.pdb, args.chain)