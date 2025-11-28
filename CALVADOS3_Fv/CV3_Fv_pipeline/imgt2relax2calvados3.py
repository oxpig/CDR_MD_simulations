import argparse
import io,re,os,shutil,sys
import json
import numpy as np
from Bio import PDB

import subprocess

from itertools import chain as iter_chain
from itertools import tee

from pdbtools import pdb_selchain,pdb_selres,pdb_tofasta,pdb_keepcoord,pdb_reres,pdb_merge,pdb_chain,pdb_reatom,pdb_delhetatm,pdb_tidy,pdb_delelem

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RELAX_path = os.path.abspath(f'{SCRIPT_DIR}/openmm_relax.py')

def create_parser():
    def parse_tuple(value):
        try:
            # Split the input string by comma
            parts = value.split(',')
            # Ensure there are exactly two elements
            if len(parts) != 2:
                raise ValueError("Expected two elements separated by a comma.")
            return tuple(parts)

        except Exception as e:
            raise argparse.ArgumentTypeError(f"Invalid input for a tuple: {value}. Error: {str(e)}")

    parser = argparse.ArgumentParser(description='Convert an AB IMGT PDB file into Calvados PDB format')

    parser.add_argument(
        '--job-id',
        type=str,
        help='Job ID to use for the output file',
        required=True
    )

    parser.add_argument(
        "--pdb-file",
        type=str,
        help="AB IMGT PDB file to convert",
        required=True
    )

    parser.add_argument(
        '--chains',
        type=parse_tuple,
        help='Chains to use for the output file. Example: "H,L". The chain order will be used also in concatenation',
        required=True
    )

    parser.add_argument(
        '--output-dir',
        type=str,
        help='Output directory for the output files',
        default='./'
    )

    parser.add_argument(
        '--max-idx',
        type=int,
        help='Maximum residue `IMGT index to keep in each chain (same for each chain)',
        default=130
    )

    parser.add_argument(
        '--length-loop',
        type=int,
        help='Length of the loop to connect the two chains (loop is GGGGS, so only multiple of 5 are allowed)',
        default=25
    )

    # parser.add_argument(
    #   '--keep-modeller-files',
    #   action='store_true',
    #   help='Keep the modeller files in the output directory',
    #   default=False
    # )

    parser.add_argument(
        '--relax',
        action='store_true',
        help='Relax the structure using openmm Amber (AF2 like)',
    )

    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite the output file if it already exists',
        default=False
    )

    return parser


class PDBediting():
    def  __init__(self):
        pass

    def __call__(self,params):

        # create pdb first chain
        chain0_pdb,chain0_dict,chain0_sequence = self.extract_chain(params['pdb_file'],params['chains'][0],params['max_idx'],1)

        # create pdb second chain
        ## evaluate the starting residue number for chain 2 including the first chain ending point and the connecting loop length
        starting_chain2 = len(chain0_sequence) + 1

        chain1_pdb,chain1_dict,chain1_sequence = self.extract_chain(params['pdb_file'],params['chains'][1],params['max_idx'],starting_chain2)

        # merge the two chains
        merged_pdb = self.merge_chains(chain0_pdb,chain1_pdb)

        if params['relax']:
            ## refine the merged pdb
            merged_pdb = self.refine_merged(merged_pdb,atom_start=1)

            ## write the merged pdb to file
            combined_pdb_loc = os.path.join(params['output_dir'],params['job_id']+f"_Fv_{params['chains'][0]}{params['chains'][1]}.pdb")
            with open(combined_pdb_loc, "w") as f:
                f.writelines(merged_pdb)

            ## run the relaxation script
            relaxed_pdb_loc = self.run_relaxation(combined_pdb_loc,params)

            ## adjust the relaxed  pdb
            merged_pdb = self.refine_relaxed(relaxed_pdb_loc,params['chains'],re2_start=starting_chain2)

        ## refine the merged pdb

        merged_pdb = self.refine_merged(merged_pdb,atom_start=1)

        ## write the merged pdb to file
        combined_pdb_loc = os.path.join(params['output_dir'],params['job_id']+f"_Fv_{params['chains'][0]}{params['chains'][1]}.pdb")
        with open(combined_pdb_loc, "w") as f:
            f.writelines(merged_pdb)

        ## remove the temporary relaxed pdb
        if params['relax']:
            os.remove(relaxed_pdb_loc)

        ## create a dictionary of the chain dictionaries and save as json:
        chain_dict = {params['chains'][0]:chain0_dict,params['chains'][1]:chain1_dict}

        chain_dict_loc = os.path.join(params['output_dir'],params['job_id']+f"_dict_Fv_{params['chains'][0]}{params['chains'][1]}.txt")

        with open(chain_dict_loc, 'w') as f:
            json.dump(chain_dict, f)

        return os.path.abspath(combined_pdb_loc),os.path.abspath(chain_dict_loc)

    def run_relaxation(self,pdb_location,params):
        command = [
                    "python",
                    RELAX_path,
                    "--job-id", params['job_id']+f"_Fv_{params['chains'][0]}{params['chains'][1]}",
                    "--pdb-file", pdb_location,
                    "--output-dir", params['output_dir']]

        try:
        # Run the script and wait for it to finish
            result = subprocess.run(command, check=True, text=True, capture_output=True)
            print("Script output:", result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error while running script: {e.stderr}")

        return os.path.join(params['output_dir'],params['job_id']+f"_Fv_{params['chains'][0]}{params['chains'][1]}_relaxed.pdb")

    def extract_resrange(self, fh, resrange):
    # Read file handle to extract residue numbers
    # Because sys.stdin is not seekable we store the
    # lines again in an iterator.

    # Options can be single numbers or ranges.

        def remove_keys_larger_than_last_set_value(d, s):
        # Convert set to a sorted list and get the last value
            last_value_of_set = sorted(s)[-1] if s else None

            if last_value_of_set is None:
            # If the set is empty, return the original dictionary
                return d

            # Create a new dictionary with keys not larger than the last value of the set
            filtered_dict = {k: v for k, v in d.items() if k <= last_value_of_set}

            return filtered_dict

        def _validate_opt_numeric(value):
            """Returns a valid numerical option or dies trying"""
            try:
                num = int(value)
            except ValueError:
                emsg = "ERROR!! Not a valid number: '{}'\n"
                sys.stderr.write(emsg.format(value))
                sys.exit(1)
            else:
                if (-999 <= num < 10000):
                    return num
                else:
                    emsg = "ERROR!! Residue numbers must be between -999 and 9999: '{}'\n"
                    sys.stderr.write(emsg.format(value))
                    sys.exit(1)

        def _validate_opt_range(value, resid_list):
            """Returns a numerical range or dies trying"""

            # Validate formatting
            if not (1 <= value.count(':') <= 2):
                emsg = "ERROR!! Residue range must be in 'a:z:s' where a and z are "
                emsg += 'optional (default to first residue and last respectively), and'
                emsg += 's is an optional step value (to return every s-th residue).\n'
                sys.stderr.write(emsg)
                sys.exit(1)

            start, end, step = None, None, 1
            slices = [_validate_opt_numeric(num)
                      if num.strip() else None for num in value.split(':')]

            if len(slices) == 3:
                start, end, step = slices
            elif len(slices) == 2:
                start, end = slices
            elif len(slices) == 1:
                if value.startswith(':'):
                    end = slices[0]
                else:
                    start = slices[0]

            # Upper/Lower limits, resid max 4 char
            if start is None:
                start = -1000
            if end is None:
                end = 10000

            # extra validation for step
            if step is None:
                step = 1
            else:
                if step < 1:
                    emsg = "ERROR!! Step value must be a positive number: '{}'\n"
                    sys.stderr.write(emsg.format(step))
                    sys.exit(1)

            # validate proper order in range
            if start > end:
                emsg = 'ERROR!! Start ({}) cannot be larger than end ({})\n'
                sys.stderr.write(emsg.format(start, end))
                sys.exit(1)

            # Build range
            bounded_resid = [r for r in resid_list if start <= r <= end]
            return bounded_resid[::step]


        buffer = iter([])
        resid_list = []
        multi_resid = {}
        records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
        prev_res = None
        for line in fh:

            if line.startswith(records):
                res_id = line[21:26]  # include chain ID
                if res_id != prev_res:
                    prev_res = res_id
                    resid_list.append(int(line[22:26]))
                    multi_resid[int(line[22:26])]=1
                    extra_res=[]
                else:
                    if str(line[26]).isalpha() and str(line[26])  not in extra_res:
                        multi_resid[int(line[22:26])]+=1
                        extra_res.append(str(line[26]))
            buffer = iter_chain(buffer, [line])

        fh = buffer
        residue_range = set()  # stores all the residues to write.
        for entry in resrange:
            if ':' in entry:
                resrange = _validate_opt_range(entry, resid_list)
                residue_range.update(resrange)
            else:
                singleres = _validate_opt_numeric(entry)
                residue_range.add(singleres)

        # Filter the multi_resid dictionary to only include residues in the residue_range
        multi_resid = remove_keys_larger_than_last_set_value(multi_resid, residue_range)

        return (fh,residue_range, multi_resid)

    def dictionary_imgt_calvados_pos(self,multi_resid_imgt,calvados_pos):
        out_dict={}
        calvado_idx=0
        calvados_pos=list(calvados_pos)
        for key in multi_resid_imgt:
            if multi_resid_imgt[key]>10:
                print('This software does not support more than 10 residues with the same index, AB skipped')
                sys.exit(1)

            for i in range(multi_resid_imgt[key]):
                if key == 112:
                    rev=multi_resid_imgt[key]-1-i
                    out_dict[float(key)+float(rev)/10.]=calvados_pos[calvado_idx]
                else:
                    out_dict[float(key)+float(i)/10.]=calvados_pos[calvado_idx]

                calvado_idx+=1
        return out_dict

    def filtered_fasta_sequence(self,fasta_gen):
            ### this takes in input a generator and return the sequences included (generator is a fasta here)
        return re.sub('\n','',"".join(line.strip() for line in fasta_gen if not line.startswith('>')))

    def extract_chain(self,pdb_loc,chain,max_res,re_start=1,extras=True):

        ## select chain
        #chain_pdb = run_selchain(open(pdb_loc,'r'), {chain})
        chain_pdb = pdb_selchain.run(open(pdb_loc,'r'), {chain})

        ## remove HETATM
        chain_pdb = pdb_delhetatm.run(chain_pdb)

        ## select residues
        inpdb,resrange_imgt,multi_resid_imgt  = self.extract_resrange(chain_pdb, [f'1:{str(max_res)}'])
        chain_pdb = pdb_selres.run(inpdb, resrange_imgt)

        ## renumber residues
        chain_pdb = pdb_reres.run(chain_pdb,re_start)
        if extras:
            # polish pdb, keeping only coordinates
            chain_pdb,extra1,extra2 = tee(pdb_keepcoord.run(chain_pdb),3)

            ## extract residue range
            _,resrange_calvados,_ = self.extract_resrange(extra1, [f'1:9999'])
            chain_dict=self.dictionary_imgt_calvados_pos(multi_resid_imgt,resrange_calvados)
            #extract fasta from pdb
            sequence_fasta = self.filtered_fasta_sequence(pdb_tofasta.run(extra2,None))

            return (chain_pdb, chain_dict,sequence_fasta)
        else:
            chain_pdb = pdb_keepcoord.run(chain_pdb)
            return (chain_pdb, None, None)

    def merge_chains(self,chain1, chain2):
        """Merge two generators into a single generator."""
        # Yield from the first generator
        yield from chain1
        # Yield from the second generator
        yield from chain2

    def refine_relaxed(self,relaxed_pdb,chains,re2_start=1):

        #openmm renames any chain to A,B so we have to switch back to the original chain names after extraction
        max_res=999
        chain1,_,_=self.extract_chain(relaxed_pdb,'A',max_res,re_start=1,extras=False)

        chain1 = pdb_chain.run(chain1,chains[0])

        chain2,_,_=self.extract_chain(relaxed_pdb,'B',max_res,re_start=re2_start,extras=False)

        chain2 = pdb_chain.run(chain2,chains[1])

        merged_pdb = self.merge_chains(chain1,chain2)

        return merged_pdb

    def refine_merged(self,merged_pdb,atom_start=1):

        # remove hydrogens if presents
        merged_pdb = pdb_delelem.run(merged_pdb,'H')

        # remve any extra lines
        merged_pdb = pdb_tidy.run(merged_pdb)

        # renumber atoms
        merged_pdb = pdb_reatom.run(merged_pdb,atom_start)

        return merged_pdb

def main(args):

    params={}
    for key in vars(args):
        params[key]=vars(args)[key]

    if os.path.abspath(params['output_dir']) == os.path.abspath(os.getcwd()):
        ## error: out directory cannot be execution one because of modeller output
        raise ValueError("Output directory cannot be the current directory. Please provide a different output directory.")
    else:
        os.makedirs(params['output_dir'], exist_ok = True)

    ### check if the output directory exist, if yes do not run unless overwrite is set to True
    if os.path.exists(os.path.join(params['output_dir'],params['job_id']+f"_Fv_{params['chains'][0]}{params['chains'][1]}.pdb")):
        if params['overwrite']:
            print('Output file already exists, overwriting')
        else:
            print('Output file already exists, skipping')
            sys.exit(0)

    pdbediting = PDBediting()
    pdb_merged_loc,dict_imgt_location=pdbediting(params)

    print('Done, file are stored in ',params['output_dir'])

if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
