# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 9:36
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : annotate.py

import json
import multiprocessing
from collections import defaultdict
import numpy as np


def extract_coord_from_gtf(gtf_file) -> tuple:
    """
    Extract mrna/exon coordinates from a GTF file.

    :param gtf_file: a GTF file with gene annotation information
    :return: a tuple of two dictionaries with mRNA and exon coordinates
    """
    mrna_coord = defaultdict(list)
    exon_coord = defaultdict(list)

    if not gtf_file.endswith('.gtf'):
        raise ValueError('The input file must be a GTF file.')

    with open(gtf_file) as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'transcript':
                transcript_id = fields[-1].split(';')[0].split()[1].strip('"')
                # add chromosome, start, end, strand as value to the dictionary
                mrna_coord[transcript_id] = [fields[0], int(fields[3]), int(fields[4]), fields[6]]
            if fields[2] == 'exon':
                transcript_id = fields[-1].split(';')[0].split()[1].strip('"')
                # add exon start and end as value to the dictionary
                exon_coord[transcript_id].append([int(fields[3]), int(fields[4])])
    return mrna_coord, exon_coord


def extract_coord_from_snp(snp_file, mrna_coord):
    """
    Extract SNP coordinates from a BED file and keep only the SNPs that are within the gene.

    :param snp_file: a BED file with SNP information
    :param mrna_coord: a dictionary with mRNA coordinates
    :return: a dictionary with SNP coordinates
    """
    snp_coord = defaultdict(list)
    gene_snp_coord = defaultdict(list)

    if snp_file is None:
        return gene_snp_coord

    if not snp_file.endswith('.bed'):
        raise ValueError('The SNP file must be a BED file.')

    with open(snp_file) as file:
        for line in file:
            fields = line.strip().split('\t')
            # add start, end as value to the dictionary for each chromosome
            snp_coord[fields[0]].append([int(fields[1]), int(fields[2])])

    for transcript_id, values in mrna_coord.items():
        if values[0] not in snp_coord:
            continue
        _, overlaps_sets = find_overlaps([[values[1], values[2]]], snp_coord[values[0]])
        if overlaps_sets is not None:
            # only keep snp coordinates
            gene_snp_coord[transcript_id] = overlaps_sets[:, 2:].tolist()
    return gene_snp_coord


def extract_coord_from_primer(primer_result):
    """
    Extract primer coordinates from a primer result.

    :param primer_result: a primer result for a sequence
    :return: a tuple of three lists with left primer coordinates, right primer coordinates and product coordinates
    """
    # all primer left array, primer right array and product array
    left_coords, right_coords, product_coords = [], [], []
    for left, right in zip(primer_result['PRIMER_LEFT'], primer_result['PRIMER_RIGHT']):
        left_coords.append([left['COORDS'][0], sum(left['COORDS']) - 1])
        right_coords.append([right['COORDS'][0] - right['COORDS'][1] + 1, right['COORDS'][0]])
        product_coords.append([left['COORDS'][0], right['COORDS'][0]])
    # Convert lists to numpy arrays
    left_coords = np.array(left_coords)
    right_coords = np.array(right_coords)
    product_coords = np.array(product_coords)
    return left_coords, right_coords, product_coords


def find_overlaps(set1, set2):
    """Find overlaps between two sets of genomic features.

    :param set1: a 2D numpy array with start and end positions of the first set of genomic features
    :param set2: a 2D numpy array with start and end positions of the second set of genomic features
    :return: a tuple of a dictionary with the number of overlaps for each interval in set1 and
    a numpy array with the overlapping intervals
    """
    # Convert lists to numpy arrays and reshape them to 2D arrays
    set1 = np.reshape(set1, (-1, 2))
    set2 = np.reshape(set2, (-1, 2))

    # Use broadcasting to compare each interval in set1 with each interval in set2
    start1, end1 = set1[:, None, 0], set1[:, None, 1]
    start2, end2 = set2[:, 0], set2[:, 1]

    # Find overlapping intervals
    overlap_indices = np.where((start1 < end2) & (end1 > start2))
    # if overlap_indices[0] is empty, return None
    if len(overlap_indices[0]) == 0:
        return None, None

    # count the number of overlaps for each interval in set1
    overlaps_counts = dict(zip(*np.unique(overlap_indices[0], return_counts=True)))

    # Retrieve overlapping intervals from set1 and set2
    overlaps_set1 = set1[overlap_indices[0]]
    overlaps_set2 = set2[overlap_indices[1]]

    # stack the two overlapping sets
    overlaps_sets = np.column_stack((overlaps_set1, overlaps_set2))

    return overlaps_counts, overlaps_sets


def get_junctions(exons):
    """Get the junctions between exons after joining exons end to end."""
    # Sort the exons by start position
    exons = exons[np.argsort(exons[:, 0])]
    # join the exons end to end where the start of the next exon is the end of the previous exon
    junctions = np.cumsum(np.diff(exons).ravel()) + exons[0, 0]

    return np.column_stack((junctions, junctions + 1))


def annotate_exon_span(primer_result):
    """
    Annotate the exon span of all primer pair of each sequence.

    :param primer_result: a primer result for a sequence
    :return: a primer result for a sequence with exon span
    """
    left_coords, right_coords, product_coords = extract_coord_from_primer(primer_result)
    # get the mrna start and end
    mrna_start, mrna_end = primer_result['MRNA_COORD'][1], primer_result['MRNA_COORD'][2]
    exons = np.asarray(primer_result['EXON_COORD'])
    # set start position of mran as 1
    if primer_result['MRNA_COORD'][3] == '-':
        exons = mrna_end - exons[:, ::-1] + 1
        junctions = get_junctions(exons)
    else:
        exons = exons - mrna_start + 1
        junctions = get_junctions(exons)

    # Check if any primer can span any exon and add the exon span to the primer result
    # Define default exon span, 0 means no exon span, 1 means exon span
    left_exon_count, _ = find_overlaps(left_coords, junctions)
    right_exon_count, _ = find_overlaps(right_coords, junctions)
    product_exon_count, _ = find_overlaps(product_coords, junctions)
    for i in range(primer_result['PRIMER_PAIR_NUM_RETURNED']):
        primer_result['PRIMER_LEFT'][i]['EXON_SPAN'] = int(left_exon_count.get(i, 0)) if left_exon_count else 0
        primer_result['PRIMER_RIGHT'][i]['EXON_SPAN'] = int(right_exon_count.get(i, 0)) if right_exon_count else 0
        primer_result['PRIMER_PAIR'][i]['EXON_SPAN'] = int(product_exon_count.get(i, 0)) if product_exon_count else 0
    return primer_result


def annotate_snp_span(primer_result):
    """Annotate the snp span of all primer pair of each sequence.

    :param primer_result: a primer result for a sequence
    :return: a primer result for a sequence with snp span
    """
    # extract primer coordinates from a primer result
    left_coords, right_coords, product_coords = extract_coord_from_primer(primer_result)
    # get the snp coordinates
    snp_coords = primer_result['SNP_COORD']
    if len(snp_coords) == 0:
        snp_coords = [[0, 0]]
    # get the mrna start and end and set start position of mran as 1
    mrna_start, mrna_end = primer_result['MRNA_COORD'][1], primer_result['MRNA_COORD'][2]
    snp_coords = np.asarray(snp_coords) - mrna_start + 1
    # Check if any primer can span any snp and add the snp span to the primer result
    # Define default snp span, 0 means no snp span, 1 means snp span
    left_snp_count, _ = find_overlaps(left_coords, snp_coords)
    right_snp_count, _ = find_overlaps(right_coords, snp_coords)
    product_snp_count, _ = find_overlaps(product_coords, snp_coords)
    for i in range(primer_result['PRIMER_PAIR_NUM_RETURNED']):
        primer_result['PRIMER_LEFT'][i]['SNP_SPAN'] = int(left_snp_count.get(i, 0)) if left_snp_count else 0
        primer_result['PRIMER_RIGHT'][i]['SNP_SPAN'] = int(right_snp_count.get(i, 0)) if right_snp_count else 0
        primer_result['PRIMER_PAIR'][i]['SNP_SPAN'] = int(product_snp_count.get(i, 0)) if product_snp_count else 0
    return primer_result


def run(primers, gtf_file, snp_file, out_file, processes=1):
    """
    Annotate the exon span and snp span of all primer pair of each sequence.

    :param primers: a JSON file with primer results
    :param gtf_file: a GTF file with gene annotation information
    :param snp_file: a BED file with SNP information (a basic BED format: chr\tstart\tend)
    :param out_file: a JSON file to save the annotated primer results
    :param processes: the number of processes to use, default is 1
    :return: a list of primer result for each sequence with exon span and snp span
    """
    print('Annotate module is running.')
    # Load the primer results from the JSON file
    primer_results = json.load(open(primers))
    # Extract mrna and exon coordinates from the GTF file
    mrna_coord, exon_coord = extract_coord_from_gtf(gtf_file)
    # Extract snp coordinates from the SNP file and keep only the snps that are within the gene
    snp_coord = extract_coord_from_snp(snp_file, mrna_coord)
    # Create a multiprocessing pool
    with multiprocessing.Pool(processes) as pool:
        # Add mrna and exon coordinates to the primer results
        for i in range(len(primer_results)):
            if not primer_results[i]:
                del primer_results[i]
            else:
                primer_results[i]['MRNA_COORD'] = mrna_coord[primer_results[i]['SEQUENCE_ID']]
                primer_results[i]['EXON_COORD'] = exon_coord[primer_results[i]['SEQUENCE_ID']]
                primer_results[i]['SNP_COORD'] = snp_coord[primer_results[i]['SEQUENCE_ID']]
        # Annotate the exon span of all primer pair of each sequence
        primer_results = pool.map(annotate_exon_span, primer_results)
        # Annotate the snp span of all primer pair of each sequence
        primer_results = pool.map(annotate_snp_span, primer_results)

    # Save the results to a json file
    with open(out_file, 'w') as file:
        json.dump(primer_results, file)
