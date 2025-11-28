import argparse
import os

from Bio.PDB import PDBParser, PDBIO, Select
# from Bio.PDB.Polypeptide import three_to_one, one_to_three
from Bio.PDB.Polypeptide import one_to_index, index_to_three
import MDAnalysis as mda

def create_parser():
    parser = argparse.ArgumentParser(description='Convert topology file to pdb file')
    parser.add_argument('--top', type=str, help='Topology file')
    parser.add_argument('--out', type=str, help='Updated topology file')
    return parser

class CAOnlySelect(Select):
    def accept_atom(self, atom):
        # Accept only CA atoms
        return atom.name == 'CA'
import MDAnalysis as mda

def update_atom_names(input_pdb, output_pdb):
    # Load the PDB file
    u = mda.Universe(input_pdb)

    # Update atom names to 'CA' for all atoms
    for atom in u.atoms:
        atom.name = 'CA'

    # Save the updated PDB file
    with mda.Writer(output_pdb, multiframe=False) as W:
        W.write(u.atoms)

def convert_residue_and_atom_names(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb)

    # Dictionary to convert one-letter codes to three-letter codes
    one_to_three_dict = {c: index_to_three(one_to_index(c)) for c in "ACDEFGHIKLMNPQRSTVWY"}
    print(one_to_three_dict)
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname
                if len(resname) == 1 and resname in one_to_three_dict:
                    three_letter_code = one_to_three_dict[resname]
                    residue.resname = three_letter_code

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

def main(args):
    
    params={}
    for key in vars(args):
        params[key]=vars(args)[key]

    outdir = os.path.dirname(params['out'])
    tmp_file = os.path.join(outdir, 'tmp.pdb')
    convert_residue_and_atom_names(params['top'], tmp_file)
    update_atom_names(tmp_file, params['out'])

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args)
