import os
import io
import time
import numpy as np
from typing import Any, Dict, Tuple, Optional, Collection, Sequence
from pathlib import Path
import argparse
import openmm
import cleanup_amber
from openmm import  unit
from openmm import app as openmm_app
from openmm.app.internal.pdbstructure import PdbStructure

ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms

def create_parser():

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
		help="PDB file to relax",
		required=True
	)
    parser.add_argument(
        "--output-dir",
        type=str,
        help="Output directory",
        default="./"
    )

    return parser

class AmberRelaxation(object):
  """Amber relaxation."""

  def __init__(self,
               max_iterations: int,
               tolerance: float,
               stiffness: float,
               max_outer_iterations: int,
               use_gpu: bool):
    """Initialize Amber Relaxer.

    Args:
      max_iterations: Maximum number of L-BFGS iterations. 0 means no max.
      tolerance: kcal/mol, the energy tolerance of L-BFGS.
      stiffness: kcal/mol A**2, spring constant of heavy atom restraining
        potential.
      max_outer_iterations: Maximum number of violation-informed relax
       iterations. A value of 1 will run the non-iterative procedure used in
       CASP14. Use 20 so that >95% of the bad cases are relaxed. Relax finishes
       as soon as there are no violations, hence in most cases this causes no
       slowdown. In the worst case we do 20 outer iterations.
      use_gpu: Whether to run on GPU.
    """
    self._max_iterations = max_iterations
    self._tolerance = tolerance
    self._stiffness = stiffness
    self._max_outer_iterations = max_outer_iterations
    self._use_gpu = use_gpu

  def process(self, *,
              pdb_path: str = None,
              ) -> Tuple[str, Dict[str, Any], Sequence[float]]:
    """Runs Amber relax on a prediction, adds hydrogens, returns PDB string."""

    pdb_string = clean_protein(pdb_path)

    out = _run_one_iteration(
        pdb_string=pdb_string,
        max_iterations=self._max_iterations,
        tolerance=self._tolerance, stiffness=self._stiffness,
        restraint_set="non_hydrogen",
        max_attempts=100,
        use_gpu=self._use_gpu)

    min_pos = out['pos']
    start_pos = out['posinit']

    min_pdb = out['min_pdb']
    #min_pdb = utils.overwrite_b_factors(min_pdb, prot.b_factors)
    #utils.assert_equal_nonterminal_atom_types(
    #    protein.from_pdb_string(min_pdb).atom_mask,
    #           prot.atom_mask)
    #violations = out['structural_violations'][
    #    'total_per_residue_violations_mask'].tolist()
    return min_pdb#, debug_data, violations

def _get_pdb_string(topology: openmm_app.Topology, positions: unit.Quantity):
  """Returns a pdb string provided OpenMM topology and positions."""
  with io.StringIO() as f:
    openmm_app.PDBFile.writeFile(topology, positions, f)
    return f.getvalue()

def _check_cleaned_atoms(pdb_cleaned_string: str, pdb_ref_string: str):
  """Checks that no atom positions have been altered by cleaning."""
  cleaned = openmm_app.PDBFile(io.StringIO(pdb_cleaned_string))
  reference = openmm_app.PDBFile(io.StringIO(pdb_ref_string))

  cl_xyz = np.array(cleaned.getPositions().value_in_unit(LENGTH))
  ref_xyz = np.array(reference.getPositions().value_in_unit(LENGTH))

  for ref_res, cl_res in zip(reference.topology.residues(),
                             cleaned.topology.residues()):
    assert ref_res.name == cl_res.name
    for rat in ref_res.atoms():
      for cat in cl_res.atoms():
        if cat.name == rat.name:
          if not np.array_equal(cl_xyz[cat.index], ref_xyz[rat.index]):
            raise ValueError(f"Coordinates of cleaned atom {cat} do not match "
                             f"coordinates of reference atom {rat}.")

def load_pdb_to_memory(file_path: str):
    """
    Loads a PDB file into memory using StringIO.
    
    Args:
        file_path (str): Path to the PDB file.
    
    Returns:
        StringIO: In-memory file-like object containing the PDB content.
    """
    with open(file_path, 'r') as pdb_file:
        pdb_string = pdb_file.read()
    
    # Create an in-memory file-like object
  #  pdb_in_memory = io.StringIO(pdb_content)
    return pdb_string

def clean_protein(
    pdb_file,
    checks: bool = True):
  """Adds missing atoms to Protein instance.

  Args:
    prot: A `protein.Protein` instance.
    checks: A `bool` specifying whether to add additional checks to the cleaning
      process.

  Returns:
    pdb_string: A string of the cleaned protein.
  """
  prot_pdb_string=load_pdb_to_memory(pdb_file)
  pdb_mem = io.StringIO(prot_pdb_string)
  # Clean pdb.
  alterations_info = {}
  fixed_pdb = cleanup_amber.fix_pdb(pdb_mem, alterations_info)
  fixed_pdb_file = io.StringIO(fixed_pdb)
  pdb_structure = PdbStructure(fixed_pdb_file)
  cleanup_amber.clean_structure(pdb_structure, alterations_info)

  # Write pdb file of cleaned structure.
  as_file = openmm_app.PDBFile(pdb_structure)
  pdb_string = _get_pdb_string(as_file.getTopology(), as_file.getPositions())
  if checks:
    _check_cleaned_atoms(pdb_string, prot_pdb_string)
  return pdb_string


def will_restrain(atom: openmm_app.Atom, rset: str) -> bool:
  """Returns True if the atom will be restrained by the given restraint set."""

  if rset == "non_hydrogen":
    return atom.element.name != "hydrogen"
  elif rset == "c_alpha":
    return atom.name == "CA"

def _add_restraints(
    system: openmm.System,
    reference_pdb: openmm_app.PDBFile,
    stiffness: unit.Unit,
    rset: str):

  """Adds a harmonic potential that restrains the system to a structure."""
  assert rset in ["non_hydrogen", "c_alpha"]

  force = openmm.CustomExternalForce(
      "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
  force.addGlobalParameter("k", stiffness)
  for p in ["x0", "y0", "z0"]:
    force.addPerParticleParameter(p)

  for i, atom in enumerate(reference_pdb.topology.atoms()):
    if will_restrain(atom, rset):
      force.addParticle(i, reference_pdb.positions[i])
  system.addForce(force)

def _openmm_minimize(
    pdb_string: str,
    max_iterations: int,
    tolerance: unit.Unit,
    stiffness: unit.Unit,
    restraint_set: str,
    use_gpu: bool):
  """Minimize energy via openmm."""
  
  pdb_file = io.StringIO(pdb_string)
  pdb = openmm_app.PDBFile(pdb_file)

  force_field = openmm_app.ForceField("amber99sb.xml")
  constraints = openmm_app.HBonds
  system = force_field.createSystem(
      pdb.topology, constraints=constraints)
  if stiffness > 0 * ENERGY / (LENGTH**2):
    _add_restraints(system, pdb, stiffness, restraint_set)

  integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
  platform = openmm.Platform.getPlatformByName("CUDA" if use_gpu else "CPU")
  simulation = openmm_app.Simulation(
      pdb.topology, system, integrator, platform)
  simulation.context.setPositions(pdb.positions)

  ret = {}
  state = simulation.context.getState(getEnergy=True, getPositions=True)
  ret["einit"] = state.getPotentialEnergy().value_in_unit(ENERGY)
  ret["posinit"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)

  simulation.minimizeEnergy(maxIterations=max_iterations,
                            tolerance=tolerance)

  state = simulation.context.getState(getEnergy=True, getPositions=True)
  ret["efinal"] = state.getPotentialEnergy().value_in_unit(ENERGY)
  ret["pos"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
  ret["min_pdb"] = _get_pdb_string(simulation.topology, state.getPositions())

  return ret

def _run_one_iteration(
    *,
    pdb_string: str,
    max_iterations: int,
    tolerance: float,
    stiffness: float,
    restraint_set: str,
    max_attempts: int,
    use_gpu: bool):
  """Runs the minimization pipeline.

  Args:
    pdb_path: A pdb path.
    max_iterations: An `int` specifying the maximum number of L-BFGS iterations.
    A value of 0 specifies no limit.
    tolerance: kcal/mol, the energy tolerance of L-BFGS.
    stiffness: kcal/mol A**2, spring constant of heavy atom restraining
      potential.
    restraint_set: The set of atoms to restrain.
    max_attempts: The maximum number of minimization attempts.
    use_gpu: Whether to run on GPU.
    exclude_residues: An optional list of zero-indexed residues to exclude from
        restraints.

  Returns:
    A `dict` of minimization info.
  """

  # Assign physical dimensions.
  tolerance = tolerance * ENERGY
  stiffness = stiffness * ENERGY / (LENGTH**2)

  start = time.time()
  minimized = False
  attempts = 0
  while not minimized and attempts < max_attempts:
    attempts += 1
    #try:
    #  logging.info("Minimizing protein, attempt %d of %d.",
    #               attempts, max_attempts)
    ret = _openmm_minimize(
        pdb_string=pdb_string, max_iterations=max_iterations,
        tolerance=tolerance, stiffness=stiffness,
        restraint_set=restraint_set,
        use_gpu=use_gpu)
    minimized = True
    #except Exception as e:  # pylint: disable=broad-except
     # logging.info(e)
  #if not minimized:
   # raise ValueError(f"Minimization failed after {max_attempts} attempts.")
  ret["opt_time"] = time.time() - start
  ret["min_attempts"] = attempts
  return ret

def patch_openmm():
    from openmm import app
    from openmm.unit import nanometers, sqrt

    # applied https://raw.githubusercontent.com/deepmind/alphafold/main/docker/openmm.patch
    # to OpenMM 7.7.1 (see PR https://github.com/openmm/openmm/pull/3203)
    # patch is licensed under CC-0
    # OpenMM is licensed under MIT and LGPL
    # fmt: off
    def createDisulfideBonds(self, positions):
        def isCyx(res):
            names = [atom.name for atom in res._atoms]
            return 'SG' in names and 'HG' not in names
        # This function is used to prevent multiple di-sulfide bonds from being
        # assigned to a given atom.
        def isDisulfideBonded(atom):
            for b in self._bonds:
                if (atom in b and b[0].name == 'SG' and
                    b[1].name == 'SG'):
                    return True

            return False

        cyx = [res for res in self.residues() if res.name == 'CYS' and isCyx(res)]
        atomNames = [[atom.name for atom in res._atoms] for res in cyx]
        for i in range(len(cyx)):
            sg1 = cyx[i]._atoms[atomNames[i].index('SG')]
            pos1 = positions[sg1.index]
    # applied https://raw.githubusercontent.com/deepmind/alphafold/main/docker/openmm.patch
    # to OpenMM 7.7.1 (see PR https://github.com/openmm/openmm/pull/3203)
    # patch is licensed under CC-0
    # OpenMM is licensed under MIT and LGPL
    # fmt: off
    def createDisulfideBonds(self, positions):
        def isCyx(res):
            names = [atom.name for atom in res._atoms]
            return 'SG' in names and 'HG' not in names
        # This function is used to prevent multiple di-sulfide bonds from being
        # assigned to a given atom.
        def isDisulfideBonded(atom):
            for b in self._bonds:
                if (atom in b and b[0].name == 'SG' and
                    b[1].name == 'SG'):
                    return True

            return False

        cyx = [res for res in self.residues() if res.name == 'CYS' and isCyx(res)]
        atomNames = [[atom.name for atom in res._atoms] for res in cyx]
        for i in range(len(cyx)):
            sg1 = cyx[i]._atoms[atomNames[i].index('SG')]
            pos1 = positions[sg1.index]
            candidate_distance, candidate_atom = 0.3*nanometers, None
            for j in range(i):
                sg2 = cyx[j]._atoms[atomNames[j].index('SG')]
                pos2 = positions[sg2.index]
                delta = [x-y for (x,y) in zip(pos1, pos2)]
                distance = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2])
                if distance < candidate_distance and not isDisulfideBonded(sg2):
                    candidate_distance = distance
                    candidate_atom = sg2
            # Assign bond to closest pair.
            if candidate_atom:
                self.addBond(sg1, candidate_atom)
    # fmt: on
    app.Topology.createDisulfideBonds = createDisulfideBonds

def relax(pdb_filename=None, use_gpu=False,max_iterations=0,stiffness=1.):#sc start #Emil added sampled_angles pass
    if "relax" not in dir():#sc
        patch_openmm()
    
    amber_relaxer = AmberRelaxation(
        max_iterations=max_iterations,
        tolerance=2.39,
        stiffness=stiffness, #10.0 originally
        max_outer_iterations=3,
        use_gpu=use_gpu)

    relaxed_pdb_lines = amber_relaxer.process(pdb_path=pdb_filename)
    return relaxed_pdb_lines

def main(args):

    params={}
    for key in vars(args):
        params[key]=vars(args)[key]
    try:   
        pdb_file = params['pdb_file']
        job_id = params['job_id']
        output_dir = os.path.abspath(params['output_dir'])

        relaxed_pdb_lines = relax(pdb_filename=pdb_file, use_gpu=False, max_iterations=0, stiffness=1.0)

        Path(os.path.join(output_dir,job_id+'_relaxed.pdb')).write_text(relaxed_pdb_lines)
        print('Relaxation completed successfully')
    except Exception as e:
        print(e)

if __name__ == "__main__":
	parser = create_parser()
	args = parser.parse_args()
	main(args)

