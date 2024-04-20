# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 9:46
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : qPrimer.py

import argparse
import logging
import time
import json
from . import design, annotate, check, visualize  # import your modules


def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Run qPrimer package.')
    parser.add_argument('--design', action='store_true', help='Run design module')
    parser.add_argument('--annotate', action='store_true', help='Run annotate module')
    parser.add_argument('--check', action='store_true', help='Run check module')
    parser.add_argument('--visualize', action='store_true', help='Run visualize module')
    parser.add_argument('--seq_file', type=str, help='Path to the sequence file with fasta format')
    parser.add_argument('--gtf_file', type=str, help='Path to the GTF file')
    parser.add_argument('--snp_file', type=str, help='Path to the SNP file with bed format')
    parser.add_argument('--database', type=str, help='Path to the database file with fasta format')
    parser.add_argument('--processes', type=int, default=1, help='Number of processes to use')
    parser.add_argument('--out_file', type=str, default='qPrimer.json',
                        help='Path to the output file with json format')
    parser.add_argument('--html_file', type=str, default='qPrimer.html', help='Path to the output html file')
    parser.add_argument('--log_file', type=str, default='qPrimer.log', help='Path to the log file')

    args = parser.parse_args()

    # Configure the logging module
    logging.basicConfig(filename=args.log_file,
                        filemode='w',
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info('Started running qPrimer package')
    logging.info('Command line arguments: {}'.format(args))
    start_time = time.time()

    # Run the qPrimer package
    primer_results = None
    if args.design:
        primer_results = design.run(args.seq_file, args.processes)
        if args.annotate:
            primer_results = annotate.run(primer_results, args.gtf_file, args.snp_file, args.processes)
        if args.check:
            primer_results = check.run(primer_results, args.database, args.processes)
        if args.visualize:
            visualize.run(primer_results, args.seq_file, args.html_file)
    else:
        raise ValueError('Please specify the design module to run with --design option.')
    if primer_results:
        # Save the results to a json file
        with open(args.out_file, 'w') as file:
            json.dump(primer_results, file)

    end_time = time.time()
    logging.info('Finished running qPrimer package')
    logging.info('Results are saved to {}'.format(args.out_file))
    logging.info('Elapsed time: {} seconds'.format(end_time - start_time))


if __name__ == '__main__':
    main()
