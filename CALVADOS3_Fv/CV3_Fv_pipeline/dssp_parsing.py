import pandas as pd
import numpy as np
import ast
import yaml


chain2num = {'H': 0, 'L': 1}

# range of loop residues
loop_boundaries=[
    [27, 38],
    [56, 65],
    [105, 117]
]

loop_boundaries_str=[
    ['27.0', '38.0'],
    ['56.0', '65.0'],
    ['105.0', '117.0']
]

fw_boundaries_str=[
    ['26.0', '39.0'],
    ['55.0', '66.0'],
    ['104.0', '118.0']
]

def dssp_header():
    # return '  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA\n'  # noqa
    return '  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN\n'


def get_info(line):
    '''Formats a line in the dssp file into a list with required
    information'''
    try:
        info = [int(line[5:10].strip()),  # residue number
                line[10],  # insertion code
                line[11],  # chain
                line[13],  # amino acid
                line[16],  # summary of secondary structure
                line[18],  # 3-turn
                line[19],  # 4-turn
                line[20],  # 5-turn
                line[21],  # bend
                line[22],  # chirality
                line[23],  # beta bridge label 1
                line[24],  # beta bridge label 2
                line[33],  # beta sheet label
                int(line[35:38].strip()),  # Solvent accessibility
                int(line[41:45].strip()),  # N-H h bond 1
                float(line[46:50].strip()),  # N-H h bond 1 energy
                int(line[52:56].strip()),  # O h bond 1
                float(line[57:61].strip()),  # O h bond 1 energy
                int(line[63:67].strip()),  # N-H h bond 2
                float(line[68:72].strip()),  # N-H h bond 2 energy
                int(line[75:78].strip()),  # O h bond 2
                float(line[79:83].strip()),  # O h bond 2 energy
                ]
        return info
    except:
        return 'invalid line'


def dssp_format_to_df(dssp_formatted):
    '''Converts a list of formatted dssp lines into a pandas dataframe'''
    col_names = ['residue',
                 'insertion_code',
                 'chain',
                 'aa',
                 'structure_summary',
                 '3-turn',
                 '4-turn',
                 '5-turn',
                 'bend',
                 'chirality',
                 'beta_bridge_label1',
                 'beta_bridge_label2',
                 'beta_sheet_label',
                 'solvent_accessibility',
                 'N_h_bond_1',
                 'N_h_bond_1_energy',
                 'O_h_bond_1',
                 'O_h_bond_1_energy',
                 'N_h_bond_2',
                 'N_h_bond_2_energy',
                 'O_h_bond_2',
                 'O_h_bond_2_energy',
                 ]
    return pd.DataFrame([x for x in list(map(get_info, dssp_formatted)) if x != 'invalid line'],
                        columns=col_names)

def get_ss_resi(dssp_file):
    header = dssp_header()
    with open(dssp_file) as f:
        dssp = f.readlines()

    protein_lines = dssp[dssp.index(header)+1:]
    protein_df = dssp_format_to_df(protein_lines)

    return protein_df.query('structure_summary in ["E", "B", "S"]').residue.values

def get_antiprallel_b_strand_labels(protein_df, same_column=False):
    '''Returns a list of all antiparallel beta strand labels in a protein'''
    # structure summary 'E' indicates residue in beta strand
    beta_strands = protein_df[protein_df.structure_summary.isin(['E', 'B'])]

    if not same_column:
        # anti parallel strands have upper case labels
        labels1 = beta_strands.beta_bridge_label1[beta_strands.beta_bridge_label1.apply(lambda x: x.isupper())].values
        labels2 = beta_strands.beta_bridge_label2[beta_strands.beta_bridge_label2.apply(lambda x: x.isupper())].values

        # return the beta sheet labels that occur in both columns
        return list(np.intersect1d(labels1, labels2))

    elif same_column:
        # anti parallel strands have upper case labels
        labels1 = beta_strands.beta_bridge_label1[beta_strands.beta_bridge_label1.apply(lambda x: x.isupper())].values
        labels2 = beta_strands.beta_bridge_label2[beta_strands.beta_bridge_label2.apply(lambda x: x.isupper())].values

        return list(set(np.concatenate([labels1, labels2])))

    else:
        raise ValueError('same_column must be a boolean')


def get_loop_idx(strands_plus_loop, b_strand_label):
    '''Extracts the longest consecutive loop region between the two
    antiparallel b-strands
    '''
    strand_p_loop_copy = strands_plus_loop.copy()
    # single column indicating if residue is in a strand
    strand_p_loop_copy['is_strand'] = ~((strand_p_loop_copy.beta_bridge_label1 != b_strand_label) &
                          (strand_p_loop_copy.beta_bridge_label2 != b_strand_label))

    # Create a boolean mask for residues not in the strand
    mask = strand_p_loop_copy.is_strand.eq(False)

    # Find the start and end index of each group of consecutive False's
    groups = mask.ne(mask.shift()).cumsum()
    ranges = strand_p_loop_copy.is_strand[mask].groupby(groups).apply(lambda x: (x.index.min(), x.index.max()))

    # Calculate the size of each range
    ranges_size = ranges.apply(lambda x: x[1] - x[0])

    # return the maximum range
    return ranges[ranges_size.idxmax()]


def find_loop_region(protein_df, b_strand_label, same_column=False):
    '''Find loop residues in protein_df dataframe'''
    # strand before and after loop have labels in different colums
    # extract the region in between first and last residues of labels in different columns
    if not same_column:
        # find the first and last residues of the b-strands
        strand_id_1 = protein_df[protein_df.beta_bridge_label1 == b_strand_label].index[[0, -1]]
        strand_id_2 = protein_df[protein_df.beta_bridge_label2 == b_strand_label].index[[0, -1]]

        # extract the loop region
        if strand_id_1[1] < strand_id_2[0]:
            loop = protein_df.loc[strand_id_1[1]+1:strand_id_2[0]-1]
        elif strand_id_2[1] < strand_id_1[0]:
            loop = protein_df.loc[strand_id_2[1]+1:strand_id_1[0]-1]
        else:
            raise ValueError('Strands overlap')

    # loops with strand labels in the same columns are allowed
    # find the longest consecutive region between strand residues
    elif same_column:
        # find the first and last residues of the b-strand
        strand_idx = protein_df[(protein_df.beta_bridge_label1 == b_strand_label) |
                                (protein_df.beta_bridge_label2 == b_strand_label)].index[[0, -1]]

        # get the b-strand and loop residues
        strands_plus_loop = protein_df.loc[strand_idx[0]:strand_idx[1]]

        loop_idx = get_loop_idx(strands_plus_loop, b_strand_label)
        loop = strands_plus_loop.loc[loop_idx[0]:loop_idx[1]]

    else:
        raise ValueError('same_column must be a boolean')

    return loop


def is_loop_valid(loop, allow_secondary_structure):
    '''Check if an extracted loop is valid'''
    # No loop residues
    if loop.empty:
        return False

    # Check that the extracted region is not exclusively beta sheet
    if (sum([x in ['E', 'H', 'G', 'I'] for x in loop.structure_summary.values]) ==
            len(loop.structure_summary.values)):
        return False

    # Check for maximum number of secondary structure residues
    # E: b-strand, H, G I: types of helixes
    if sum([x in ['E', 'H', 'G', 'I'] for x in loop.structure_summary.values]) > allow_secondary_structure:
        return False

    # passed all tests
    return True


def get_valid_loops(protein_df, b_strand_labels, pdb, allow_secondary_structure, same_column=False):
    '''Returns a list of valid loops between antiparallel b-strands in a protein

    A valid loop is defined as a loop that does not contain any secondary structure
    or if allow_secondary_structure is True, a loop that contains max x residues of
    secondary structure elements between the end of the first b-strand and the start
    of the second b-strand'''
    valid_loops = []
    non_valid_count = 0
    for b_strand_label in b_strand_labels:
        try:
            loop = find_loop_region(protein_df, b_strand_label, same_column=same_column)
        except: # if anything breaks the function, ignore
            continue

        if not is_loop_valid(loop, allow_secondary_structure):
            non_valid_count += 1
            continue

        valid_loops.append([pdb,
                            loop.chain.values[0],
                            list(loop.insertion_code.values),
                            list(loop.residue.values),
                            list(loop.aa.values),
                            len(loop)])

    return valid_loops


def process_loop_region(grouped_df, pdb, allow_secondary_structure, same_column=False):
    '''Extract CDR loop from df with approximate region of loop'''
    grouped_df = grouped_df.reset_index(drop=True)
    b_strand_labels = get_antiprallel_b_strand_labels(grouped_df, same_column=same_column)

    return get_valid_loops(grouped_df, b_strand_labels, pdb, allow_secondary_structure, same_column=same_column)


def get_loop_idx_by_summary(strands_plus_loop):
    '''Extracts the longest consecutive loop region between the strands'''
    strand_p_loop_copy = strands_plus_loop.copy()
    # single column indicating if residue is in a strand
    strand_p_loop_copy['is_strand'] = (strand_p_loop_copy.structure_summary.isin(['E', 'B']))

    # Create a boolean mask for residues not in the strand
    mask = strand_p_loop_copy.is_strand.eq(False)

    # Find the start and end index of each group of consecutive False's
    groups = mask.ne(mask.shift()).cumsum()
    ranges = strand_p_loop_copy.is_strand[mask].groupby(groups).apply(lambda x: (x.index.min(), x.index.max()))

    # Calculate the size of each range
    ranges_size = ranges.apply(lambda x: x[1] - x[0])

    # return the maximum range
    return ranges[ranges_size.idxmax()]


def process_loop_region_by_summary(grouped_df, pdb, allow_secondary_structure):
    '''Extract CDR loop from df with approximate region of loop. Use summary of
    secondary structure only. I.e. extract any loop region between two consecutive
    b-strands, they do not have to be bonded to each other'''
    protein_df = grouped_df.reset_index(drop=True)

    try:
        # find the first and last residues of the b-strand
        strand_idx = protein_df[protein_df.structure_summary.isin(['E', 'B'])].index[[0, -1]]

        # get the b-strand and loop residues
        strands_plus_loop = protein_df.loc[strand_idx[0]:strand_idx[1]]

        loop_idx = get_loop_idx_by_summary(strands_plus_loop)
        loop = strands_plus_loop.loc[loop_idx[0]:loop_idx[1]]

        if not is_loop_valid(loop, allow_secondary_structure):
            non_valid_count += 1
            return 'no loops'

        valid_loops = [pdb,
                       loop.chain.values[0],
                       list(loop.insertion_code.values),
                       list(loop.residue.values),
                       list(loop.aa.values),
                       len(loop)]

        return valid_loops

    except:
        return 'no loops'


def extract_CDR_from_results(results):
    '''Find the loop which corresponds to a CDR'''
    if len(results) > 1:  # multiple loops are found
        results_df = pd.DataFrame(results, columns=['pdb','chain', 'insertion_code','resi', 'aa', 'len'])
        results_df['loop_contains_resi'] = results_df.resi.apply(lambda x: (108 in np.array(x) or  # CDR3
                                                                            58 in np.array(x) or  # CDR2
                                                                            35 in np.array(x)))  # CDR1

        # select if only one conatins residue 107
        if results_df[results_df['loop_contains_resi']].shape[0] == 1:
            result = results[np.argwhere(results_df.loop_contains_resi.values)[0][0]]
        else:  # select the longer
            result = results[results_df['len'].idxmax()]
    elif not results:  # no loops found
        return 'no loops'
    else:  # only one loop found
        result = results[0]

    return result


def get_dssp_CDR(row,
                 chain_name,
                 loop_boundaries=[103, 120],
                 same_column=False,
                 structure_summary_only=False,
                 return_loop_df=False):
    '''Extract CDR3 from DSSP file.

    Function to apply to CDR_seqs dataframe to extract CDRH3 or CDRL3
    from DSSP secondary structure annotation file.

    Args:
    row: pd.Series, row in CDR_seqs dataframe
    chain_name: str, 'H' or 'L'
    loop_boundaries: list, [start, end] of region to parse for loop,
                     select depending on which CDR to extract
    same_column: bool, if True, extract loops between strands occuring in
                 the same column in the DSSP file
    structure_summary_only: bool, if True, extract loop based on summary of
                            secondary structure only, for H1
    return_loop_df: bool, if True, return the dataframe of the loop region
    '''
    header = dssp_header()
    with open(row.DSSP_file) as f:
        dssp = f.readlines()

    # extract region around CDR3 loop to inspect for loops in DSSP file
    protein_lines = dssp[dssp.index(header)+1:]
    protein_df = dssp_format_to_df(protein_lines)
    loop_region = protein_df[(protein_df['chain'] == row.chains[chain2num[chain_name]]) &
                             (protein_df['residue'] >= loop_boundaries[0]) &
                             (protein_df['residue'] <= loop_boundaries[1])]

    # return df of loop region
    if return_loop_df:
        return loop_region

    # extract loop
    if not structure_summary_only:
        results = process_loop_region(loop_region, row.pdb, 20, same_column=same_column)
        result = extract_CDR_from_results(results)
    elif structure_summary_only:
        result = process_loop_region_by_summary(loop_region, row.pdb, 20)
    if result == 'no loops':
        return 'no loops'
    else:
        return result[3], result[2], result[4], result[5]


def homogenise_DSSP_CDRs(df, CDR):
    '''Make sure CDRs with identical Fv sequences have the same DSSP CDRs.
    If multiple different defintions are found select the longest'''
    groups = []
    for i, group in enumerate(df.groupby('Fv')):
        seq, grouped_df = group

        if grouped_df[f'DSSP_{CDR}_seq'].nunique() > 1:
            grouped_df[f'DSSP_{CDR}_length'] = pd.to_numeric(grouped_df[f'DSSP_{CDR}_length'], errors='coerce')
            max_index = grouped_df[f'DSSP_{CDR}_length'].idxmax()
            select = grouped_df.loc[max_index]
            grouped_df[f'DSSP_{CDR}_seq'] = select[f'DSSP_{CDR}_seq']
            grouped_df[f'DSSP_{CDR}_insertion'] = [select[f'DSSP_{CDR}_insertion']] * len(grouped_df)
            grouped_df[f'DSSP_{CDR}_resi'] = [select[f'DSSP_{CDR}_resi']] * len(grouped_df)
            grouped_df[f'DSSP_{CDR}_length'] = select[f'DSSP_{CDR}_length']
            grouped_df[f'DSSP_{CDR}_start'] = select[f'DSSP_{CDR}_start']
            grouped_df[f'DSSP_{CDR}_end'] = select[f'DSSP_{CDR}_end']
            groups.append(grouped_df)
        else:
            groups.append(grouped_df)

    return pd.concat(groups).sort_index()


def get_h_bond_df(row, CDR='CDRH3'):
    '''Extract hydrogen bond information from DSSP file for CDR loop'''
    col_sele = ['residue',
                'insertion_code',
                'chain',
                'aa',
                'N_h_bond_1', 'N_h_bond_1_energy',
                'O_h_bond_1', 'O_h_bond_1_energy',
                'N_h_bond_2', 'N_h_bond_2_energy',
                'O_h_bond_2', 'O_h_bond_2_energy'
                ]
    boundaries = {'CDRH1': [27, 38],
                  'CDRH2': [56, 65],
                  'CDRH3': [105, 117],
                  'CDRL1': [27, 38],
                  'CDRL2': [56, 65],
                  'CDRL3': [105, 117]}

    df = get_dssp_CDR(row,
                      CDR[3],
                      loop_boundaries=boundaries[CDR],
                      return_loop_df=True)
    return df.loc[:, col_sele]

def get_h_bonded_residues_CDR(protein_df, chain_label='H', loop_boundary=[105,117], energy_threshold=-0.75):
    loop_region = protein_df[(protein_df['chain'] == chain_label) &
                             (protein_df['residue'] >= loop_boundary[0]) &
                             (protein_df['residue'] <= loop_boundary[1])]
    
    constrained_resi = loop_region[(loop_region.N_h_bond_1_energy <= energy_threshold) | 
                                   (loop_region.O_h_bond_1_energy <= energy_threshold) |
                                   (loop_region.N_h_bond_2_energy <= energy_threshold) |
                                   (loop_region.O_h_bond_2_energy <= energy_threshold)]
    return list(constrained_resi.apply(lambda row: str(row['residue']) + str(row['insertion_code']), axis=1).values)

def get_h_bonded_residues(dssp_file, loop_boundaries, chains='HL', energy_threshold=-0.75):
    with open(dssp_file) as f:
        dssp = f.readlines()

    # extract region around CDR3 loop to inspect for loops in DSSP file
    header = dssp_header()
    protein_lines = dssp[dssp.index(header)+1:]
    protein_df = dssp_format_to_df(protein_lines)

    h_bonded_resi = {}
    for chain in chains:
        h_bonded_resi[chain] = []
        for boundary in loop_boundaries:
            h_bonded_resi[chain].append(
                get_h_bonded_residues_CDR(protein_df, chain_label=chain, loop_boundary=boundary, energy_threshold=energy_threshold)
            )
    
    return h_bonded_resi

def load_map(map_file):
    with open(map_file, 'r') as file:
        data = file.read()
        imgt2sequential = ast.literal_eval(data)
    return imgt2sequential

def sequential_loop_boundaries(imgt2sequential, chains):

    seq_loop_boundaries = {}
    for chain in chains:
        if chain != '-': # ignore loop for - chain
            loop_boundaries_c = [] 
            for loop_boundary in loop_boundaries_str:
                loop_boundaries_c.append(
                    [imgt2sequential[chain][loop_boundary[0]], 
                    imgt2sequential[chain][loop_boundary[1]]]
                )
            seq_loop_boundaries[chain] = loop_boundaries_c

    return seq_loop_boundaries

def h_bonds_and_energy(loop_region, energy_threshold):
    i_jrel_energy = []

    for idx, row in loop_region.iterrows():
        if row['N_h_bond_1_energy'] < energy_threshold:
            i_jrel_energy.append((row['residue'], row['N_h_bond_1'], row['N_h_bond_1_energy']))
        if row['O_h_bond_1_energy'] < energy_threshold:
            i_jrel_energy.append((row['residue'], row['O_h_bond_1'], row['O_h_bond_1_energy']))
        if row['N_h_bond_2_energy'] < energy_threshold:
            i_jrel_energy.append((row['residue'], row['N_h_bond_2'], row['N_h_bond_2_energy']))
        if row['O_h_bond_2_energy'] < energy_threshold:
            i_jrel_energy.append((row['residue'], row['O_h_bond_2'], row['O_h_bond_2_energy']))

    return i_jrel_energy

class OneLineListDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False, *args, **kwargs):
        return super(OneLineListDumper, self).increase_indent(flow=flow, indentless=indentless, *args, **kwargs)
