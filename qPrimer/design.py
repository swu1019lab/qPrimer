# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 9:36
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : design.py

import json
import configparser
import multiprocessing
import numpy as np
import pandas as pd
from Bio import SeqIO
import primer3

# Define the global and sequence arguments for Primer3
GLOBAL_ARGS = {
    'PRIMER_NUM_RETURN': 10,
    'PRIMER_PICK_LEFT_PRIMER': 1,
    'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_MIN_GC': 40,
    'PRIMER_MAX_GC': 60,
    'PRIMER_OPT_GC_PERCENT': 50,
    'PRIMER_OPT_SIZE': 22,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 28,
    'PRIMER_OPT_TM': 60,
    'PRIMER_MIN_TM': 58,
    'PRIMER_MAX_TM': 64,
    'PRIMER_PAIR_MAX_DIFF_TM': 3,
    'PRIMER_PRODUCT_SIZE_RANGE': [[80, 300]],
    'PRIMER_EXPLAIN_FLAG': 1,
    'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
    'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT': 1,
    'PRIMER_ANNEALING_TEMP': 50,
    'PRIMER_MAX_POLY_X': 5,
    'PRIMER_MAX_SELF_ANY': 4,
    'PRIMER_MAX_SELF_END': 3,
    'PRIMER_PAIR_MAX_COMPL_ANY': 4,
    'PRIMER_PAIR_MAX_COMPL_END': 3,
    'PRIMER_MAX_END_STABILITY': 9
}

SEQ_ARGS = {}


def load_config(config_file) -> tuple[dict, dict]:
    # Load the configuration file
    config = configparser.ConfigParser()
    config.read(config_file)

    global_conf = dict(config['GLOBAL'])
    global_args = {}

    seq_conf = dict(config['SEQUENCE'])
    seq_args = {}

    # Convert the values to the correct types
    for key in global_conf:
        try:
            global_args[key.upper()] = int(global_conf[key])
        except ValueError:
            try:
                global_args[key.upper()] = float(global_conf[key])
            except ValueError:
                global_args[key.upper()] = global_conf[key]

    if global_args.get('PRIMER_PRODUCT_SIZE_RANGE'):
        global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [
            [int(x) for x in size.split('-')] for size in global_args['PRIMER_PRODUCT_SIZE_RANGE'].split()
        ]

    for key in seq_conf:
        if key.upper() == 'SEQUENCE_ID':
            seq_args[key.upper()] = seq_conf[key]
        elif key.upper() == 'SEQUENCE_TEMPLATE':
            seq_args[key.upper()] = seq_conf[key]
        elif key.upper() == 'SEQUENCE_INCLUDED_REGION':
            seq_args[key.upper()] = [[int(x.split(",")[0]), int(x.split(",")[1])] for x in seq_conf[key].split()]
        elif key.upper() == 'SEQUENCE_TARGET':
            seq_args[key.upper()] = [[int(x.split(",")[0]), int(x.split(",")[1])] for x in seq_conf[key].split()]
        elif key.upper() == 'SEQUENCE_EXCLUDED_REGION':
            seq_args[key.upper()] = [[int(x.split(",")[0]), int(x.split(",")[1])] for x in seq_conf[key].split()]
        elif key.upper() == 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST':
            seq_args[key.upper()] = [[int(e) if e else -1 for e in x.split(',')] for x in seq_conf[key].split(";")]
        elif key.upper() == 'SEQUENCE_OVERLAP_JUNCTION_LIST':
            seq_args[key.upper()] = [int(x) for x in seq_conf[key].split()]
        else:
            seq_args[key.upper()] = seq_conf[key]
    return global_args, seq_args


def design_primers(sequence) -> dict:
    # Define the sequence template
    seq_args = {
        'SEQUENCE_ID': sequence.id,
        'SEQUENCE_TEMPLATE': str(sequence.seq),
    }

    seq_args.update(SEQ_ARGS)

    # check the sequence length before running Primer3
    if len(sequence.seq) < np.min(GLOBAL_ARGS['PRIMER_PRODUCT_SIZE_RANGE']):
        print(f"Sequence {sequence.id} is too short to design primers.")
        return {}

    # Run Primer3
    results = primer3.design_primers(seq_args, GLOBAL_ARGS)

    # Add the sequence information to the results
    results['SEQUENCE_ID'] = sequence.id
    results['SEQUENCE_TEMPLATE'] = str(sequence.seq)

    return results


def run(sequence_file, config_file, out_name, out_csv, processes) -> None:
    print(f"Designing primers for sequences in {sequence_file}")
    # Parse the sequence file
    sequences = SeqIO.parse(sequence_file, 'fasta')
    # Load the configuration file
    if config_file:
        _global_args, _seq_args = load_config(config_file)
        GLOBAL_ARGS.update(_global_args)
        SEQ_ARGS.update(_seq_args)

    # Create a multiprocessing pool
    with multiprocessing.Pool(processes) as pool:
        # Design primers for each sequence
        results = pool.map(design_primers, sequences)

    # Remove empty results
    primer_results = results.copy()
    for i in range(len(results)):
        if not results[i] or results[i]['PRIMER_PAIR_NUM_RETURNED'] == 0:
            # Remove empty or failed results
            del primer_results[i]
    if len(primer_results) == 0:
        print("No primers were designed!!!")
        return

    # Save the results to a json file
    with open(out_name + ".json", 'w') as file:
        json.dump(primer_results, file)

    print(f"Results saved to {out_name}.json")
    # Save main results to a csv file
    if out_csv:
        columns = ['PENALTY', 'SEQUENCE', 'COORDS', 'TM', 'GC_PERCENT', 'END_STABILITY']
        df_f_list, df_r_list = [], []
        for res in primer_results:
            df_f_list.append(
                pd.DataFrame.from_records(res['PRIMER_LEFT'], columns=columns))
            df_r_list.append(
                pd.DataFrame.from_records(res['PRIMER_RIGHT'], columns=columns).assign(ID=res['SEQUENCE_ID']))
        df_f, df_r = pd.concat(df_f_list, ignore_index=True), pd.concat(df_r_list, ignore_index=True)
        df_f.reset_index().merge(df_r.reset_index(), on='index', suffixes=('_F', '_R')).to_csv(out_name + ".csv",
                                                                                               index=False)
        print(f"Results saved to {out_name}.csv")
