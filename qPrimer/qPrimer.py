# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 9:46
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : qPrimer.py

import argparse
import logging
import os.path
import time
from . import design, annotate, check, visualize  # import your modules


def main():
    # Create a parser for the main command
    parser = argparse.ArgumentParser(description='Run qPrimer package.')
    parser.add_argument('--log_file', type=str, default='qPrimer.log', help='Path to the log file')
    parser.add_argument('--out_name', type=str, default='qPrimer',
                        help='Prefix name of the output file')

    subparsers = parser.add_subparsers(dest='module', help='Specify the module to run')

    # Create a parser for the design module
    design_parser = subparsers.add_parser('design', help='Run design module')
    design_parser.add_argument('--seq_file', required=True, type=str,
                               help='Path to the sequence file with fasta format')
    design_parser.add_argument('--ini_file', type=str,
                               help='Path to the Primer3 configuration file with ini format')
    design_parser.add_argument('--processes', type=int, default=1,
                               help='Number of processes to use')

    # Create a parser for the annotate module
    annotate_parser = subparsers.add_parser('annotate', help='Run annotate module')
    annotate_parser.add_argument('--primers', required=True, type=str,
                                 help='Path to the primers results with json format')
    annotate_parser.add_argument('--gtf_file', required=True, type=str,
                                 help='Path to the GTF file')
    annotate_parser.add_argument('--snp_file', required=True, type=str,
                                 help='Path to the SNP file with bed format')
    annotate_parser.add_argument('--processes', type=int, default=1,
                                 help='Number of processes to use')

    # Create a parser for the check module
    check_parser = subparsers.add_parser('check', help='Run check module')
    check_parser.add_argument('--primers', required=True, type=str,
                              help='Path to the primers results with json format')
    check_parser.add_argument('--database', required=True, type=str,
                              help='Path to the database file with fasta format')
    check_parser.add_argument('--processes', type=int, default=1,
                              help='Number of processes to use')

    # Create a parser for the visualize module
    visualize_parser = subparsers.add_parser('visualize', help='Run visualize module')
    visualize_parser.add_argument('--primers', required=True, type=str,
                                  help='Path to the primers results with json format')
    visualize_parser.add_argument('--seq_file', required=True, type=str,
                                  help='Path to the sequence file with fasta format')

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
    if os.path.exists(args.out_name + '.json'):
        args.out_name = args.out_name + '_' + str(int(time.time()))

    if args.module == 'design':
        design.run(args.seq_file, args.ini_file, args.out_name, args.processes)
    elif args.module == 'annotate':
        annotate.run(args.primers, args.gtf_file, args.snp_file, args.out_name, args.processes)
    elif args.module == 'check':
        check.run(args.primers, args.database, args.out_name, args.processes)
    elif args.module == 'visualize':
        visualize.run(args.primers, args.seq_file, args.out_name)
    else:
        parser.print_help()

    end_time = time.time()
    logging.info('Finished running qPrimer package')
    logging.info('Results are saved to {}.json or {}.html'.format(args.out_name, args.out_name))
    logging.info('Elapsed time: {} seconds'.format(end_time - start_time))


if __name__ == '__main__':
    main()
