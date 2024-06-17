# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 9:36
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : check.py

import json
import subprocess
from multiprocessing import Pool
from collections import Counter
from io import StringIO
import os.path


def isPcr(query, database, out='bed'):
    """
    Run isPCR to check the specificity of primer pairs.

    :param query:
    :param database:
    :param out:
    :return:
    """
    isPcr_path = os.path.join(os.path.dirname(__file__), 'bin/isPcr')
    # check isPcr whether installed
    try:
        subprocess.run([isPcr_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        raise FileNotFoundError('isPcr is not installed.')
    # isPCR command
    cmd = "{} {} stdin stdout -out={}".format(isPcr_path, database, out)
    process = subprocess.Popen(cmd.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    stdout, stderr = process.communicate(input=query.getvalue().encode())
    # Parse the isPCR results
    count = Counter()
    for line in stdout.decode().split('\n'):
        if line:
            fields = line.strip().split('\t')
            count[int(fields[3])] += 1
    return count


def check_specificity(primer_result):
    # Define the isPCR query
    query = StringIO()
    # Parse all primer results from one sequence
    name = primer_result['SEQUENCE_ID']
    database = primer_result['DATABASE']
    for i in range(primer_result['PRIMER_PAIR_NUM_RETURNED']):
        forward_primer = primer_result['PRIMER_LEFT'][i]['SEQUENCE']
        reverse_primer = primer_result['PRIMER_RIGHT'][i]['SEQUENCE']
        query.write("{}\t{}\t{}\n".format(str(i), forward_primer, reverse_primer))

    # Run isPCR
    count = isPcr(query, database, out='bed')

    # Update the primer result simply by adding a new key to the dictionary
    for i in range(primer_result['PRIMER_PAIR_NUM_RETURNED']):
        primer_result['PRIMER_PAIR'][i]['SPECIFICITY'] = count[i]
        # print('Primer pair {} of gene {} has {} specific hits'.format(i, name, count[i]))

    return primer_result


def run(primers, database, out_file, processes):
    print('Check module is running.')
    # Load the primer results from the json file
    primer_results = json.load(open(primers))
    # Create a pool of processes
    with Pool(processes=processes) as pool:
        # Add the database to the primer results
        for i in range(len(primer_results)):
            if not primer_results[i]:
                # Remove empty results
                del primer_results[i]
            else:
                primer_results[i]['DATABASE'] = database
        # Check the specificity of all primer pair of each sequence
        primer_results = pool.map(check_specificity, primer_results)

    # Save the results to a json file
    with open(out_file, 'w') as file:
        json.dump(primer_results, file)
