# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 9:36
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : design.py

import multiprocessing
import numpy as np
from Bio import SeqIO
import primer3


def design_primers(sequence) -> dict:
    # Define the sequence template
    seq_args = {
        'SEQUENCE_ID': sequence.id,
        'SEQUENCE_TEMPLATE': str(sequence.seq),
    }

    # Define the global arguments for Primer3
    global_args = {
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

    # check the sequence length before running Primer3
    if len(sequence.seq) < np.min(global_args['PRIMER_PRODUCT_SIZE_RANGE']):
        print(f"Sequence {sequence.id} is too short to design primers.")
        return {}

    # Run Primer3
    results = primer3.design_primers(seq_args, global_args)

    # Add the sequence information to the results
    results['SEQUENCE_ID'] = sequence.id
    results['SEQUENCE_TEMPLATE'] = str(sequence.seq)

    return results


def run(sequence_file, processes) -> list:
    # Parse the sequence file
    sequences = SeqIO.parse(sequence_file, 'fasta')

    # Create a multiprocessing pool
    with multiprocessing.Pool(processes) as pool:
        # Design primers for each sequence
        results = pool.map(design_primers, sequences)

    # Return the results
    return results
