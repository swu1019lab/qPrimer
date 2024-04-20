# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 17:04
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_check.py

import unittest
from qPrimer import check
from io import StringIO


class TestCheck(unittest.TestCase):
    def test_isPcr(self):
        # Create a test query
        query = StringIO()
        query.write("{}\t{}\t{}\n".format('test1', 'CGTTCGTCGTTGGTTACACG', 'TCTCCGCTTCACTTGGTTCC'))
        query.write("{}\t{}\t{}\n".format('test2', 'GCGTTCGTCGTTGGTTACAC', 'TCCGCTTCACTTGGTTCCAA'))

        # Run the isPcr function
        result = check.isPcr(query, 'test.fa', out='bed')
        # Check the result
        self.assertIsInstance(result, str)

    def test_check_specificity(self):
        # Create a test primer result
        primer_result = {
            'SEQUENCE_ID': 'gene1',
            'PRIMER_PAIR_NUM_RETURNED': 2,
            'DATABASE': 'test.fa',
            'PRIMER_LEFT': [{'SEQUENCE': 'CGTTCGTCGTTGGTTACACG'},
                            {'SEQUENCE': 'GCGTTCGTCGTTGGTTACAC'}],
            'PRIMER_RIGHT': [{'SEQUENCE': 'TCTCCGCTTCACTTGGTTCC'},
                             {'SEQUENCE': 'TCCGCTTCACTTGGTTCCAA'}]
        }

        # Run the check_specificity function
        result = check.check_specificity(primer_result)
        # Check the result
        self.assertIsInstance(result, str)

    def test_run(self):
        # Create a test primer results
        primer_results = [
            {
                'SEQUENCE_ID': 'gene1',
                'PRIMER_PAIR_NUM_RETURNED': 2,
                'PRIMER_LEFT': [{'SEQUENCE': 'CGTTCGTCGTTGGTTACACG'},
                                {'SEQUENCE': 'GCGTTCGTCGTTGGTTACAC'}],
                'PRIMER_RIGHT': [{'SEQUENCE': 'TCTCCGCTTCACTTGGTTCC'},
                                 {'SEQUENCE': 'TCCGCTTCACTTGGTTCCAA'}]
            },
            {
                'SEQUENCE_ID': 'gene2',
                'PRIMER_PAIR_NUM_RETURNED': 2,
                'PRIMER_LEFT': [{'SEQUENCE': 'AATGGCGGTGGCAATAAGGA'},
                                {'SEQUENCE': 'TTCCGGTTGGGGGATCTACT'}],
                'PRIMER_RIGHT': [{'SEQUENCE': 'AGTAGATCCCCCAACCGGAA'},
                                 {'SEQUENCE': 'CTTGAGAGCAGCGACAAGGA'}]
            }
        ]

        # Run the run function with a test database and 1 process
        result = check.run(primer_results, 'test.fa', 1)
        # Check the result
        self.assertIsInstance(result, list)
        for res in result:
            self.assertIsInstance(res, str)
