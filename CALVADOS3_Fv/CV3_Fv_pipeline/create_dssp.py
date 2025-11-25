'''Script that produce a DSSP file a PDB file.'''

import subprocess
import os
import argparse


def reformat_pdbs_for_dssp_processing(in_file, out_file):
    with open(in_file) as f:
        lines = f.readlines()

    header = ['CRYST1\n']
    lines_select = [line for line in lines if line.startswith("HETATM") or
                    line.startswith("ATOM")]
    lines_new = header + lines_select

    with open(out_file, 'w') as f:
        for line in lines_new:
            f.write(line)


def run_dssp(in_file, out_file):
    command = f"mkdssp -i {in_file} -o {out_file}"
    subprocess.Popen(command,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     shell=True)


def main(infile, outdir):

    pdb_file_processed = os.path.join(outdir, os.path.basename(infile.split('.')[0] + 'struc4dssp.pdb'))
    dssp_file =  os.path.join(outdir, os.path.basename(infile)).replace('.pdb', '.dssp')

    reformat_pdbs_for_dssp_processing(infile, pdb_file_processed)
    run_dssp(pdb_file_processed, dssp_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Script that produces a DSSP file from a PDB file.')
    parser.add_argument('--infile', help='Input PDB file')
    parser.add_argument('--outdir', help='Output dir for processed pdb and DSSP')
    args = parser.parse_args()

    main(args.infile, args.outdir)
