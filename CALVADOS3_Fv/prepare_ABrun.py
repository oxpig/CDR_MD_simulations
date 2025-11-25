import os
from calvados.cfg import Config, Job, Components
import subprocess
import argparse

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare script for Calvados2COM AB simulation')
    parser.add_argument('--temp', type=float, help='Temperature in K',default=293)
    parser.add_argument('--ionic', type=float, help='Ionic strength in molar',default=0.15)
    parser.add_argument('--pH', type=float, help='pH',default=7.0)
    parser.add_argument('--k-harmonic', type=float, help='Restraint force constant',default=700)
    parser.add_argument('--fresidues-file', type=str, help='Residue definitions file',default='/projects/prism/people/bqm193/software/CALVADOS/residues_C3.csv')
    parser.add_argument('--wfreq', type=int, help='DCD writing frequency',default=5000)
    parser.add_argument('--steps', type=int, help='Number of simulation steps (total: wfreq*steps)',default=1010)
    parser.add_argument('--platform', type=str, help='Platform (CPU or CUDA)',default='CPU')
    parser.add_argument('--threads', type=int, help='Number of threads (if CPU)',default=4)
    parser.add_argument('--verbose', type=bool, help='Verbose',default=True)
    parser.add_argument('--custom-restraints-file', default='', type=str, help='Restraints for custom pairs - columns: i,j,parameter')
    
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument('--input-pdb', type=str, help='Input PDB for the simulation',required=True)
    requiredNamed.add_argument('--out-dir', type=str, help='Output directory',required=True)
    requiredNamed.add_argument('--fdomains-file', type=str, help='Restrainted folded regions definition file',required=True)
    
    return parser

def extract_file_info(file_path):
    # Get the absolute path of the folder containing the file
    folder_path = os.path.abspath(os.path.dirname(file_path))
    
    # Get the file name with extension
    file_name_with_ext = os.path.basename(file_path)
    
    # Get the file name without extension
    file_name_without_ext = os.path.splitext(file_name_with_ext)[0]

    return folder_path, file_name_without_ext, file_name_with_ext

def main(args):

    params={}
    for key in vars(args):
        params[key] = vars(args)[key]

    pdb_folder, pdb_name, pdb_fullname = extract_file_info(params['input_pdb'])
    print(f'pdb_folder: {pdb_folder}',f'pdb_name: {pdb_name}',f'pdb_fullname: {pdb_fullname}')
    cwd = os.getcwd()
    
    config = Config(
      # GENERAL
      sysname = pdb_name,
      box = [25., 25., 150.], # nm
      temp = params['temp'], # K
      ionic = params['ionic'], # molar
      pH = params['pH'],
      topol = 'slab',
    
      # RUNTIME SETTINGS
      wfreq = params['wfreq'], # dcd writing frequency, 1 = 10fs
      steps = params['wfreq']*params['steps'], #5000*1010, # number of simulation steps
      runtime = 0, # overwrites 'steps' keyword if > 0
      platform = params['platform'], # 'CPU' or 'CUDA'
      threads = params['threads'], # number of threads
      restart = None,
      verbose = params['verbose'], # verbose output
    )
    
    # PATH
    path = os.path.join(cwd,params['out_dir'])
    subprocess.run(f'mkdir -p {path}',shell=True)
    
    config.write(path,name='config.yaml')
    
    components = Components(
      # Defaults
      molecule_type = 'protein',
      nmol = 1, # number of molecules
      restraint = True, # apply restraints
      charge_terminus = 'both', # charge N or C or both
      
      # RESTRAINTS
      restraint_type = 'harmonic', # harmonic or go
      use_com = True, # apply on centers of mass instead of CA
      colabfold = 1, # PAE format (EBI AF=0, Colabfold=1&2)
      k_harmonic = params['k_harmonic'], # force harmonic for restraint
      
      #k_custom=[000,000,000], # TO IMPLEMENT ORIGINAL CV3: force harmonic for custom restraints order: (i) lower, (ii) higher than condition ands (iii) extra restraints for a positions
      #k_custom=[700,500,250], # force harmonic for custom restraints order: (i) lower, (ii) higher than condition ands (iii) extra restraints for a positions
      k_custom=[350,175,175],
      #k_custom=[175,90,90],
      #k_custom=[700,700,700], # TO TRY FULL HB:  force harmonic for custom restraints order: (i) lower, (ii) higher than condition ands (iii) extra restraints for a positions
      # INPUT
      fresidues = params['fresidues_file'], # residue definitions
      fdomains = os.path.abspath(params['fdomains_file']), # domain definitions (harmonic restraints)
      pdb_folder = pdb_folder, # directory for pdb and PAE files
      custom_restraints_file= os.path.abspath(params['custom_restraints_file']), # restraints for custom pairs: i,j,parameter
      )
    components.add(name=pdb_name)
    
    components.write(path,name='components.yaml')

if __name__ == '__main__':

    parser = create_parser()
    args = parser.parse_args()
    main(args)

