# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 10:32
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : test_design.py

import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from qPrimer import design


class TestDesign(unittest.TestCase):
    def test_design_primers(self):
        # Create a test sequence
        sequence = SeqRecord(Seq("".join(["ATGGATTGGCAAGGACAGAAACTAGCGGAGCAGCTGATGCAGATCCTGCTCTTGATCGCC",
                                          "GCCGTTGTGGCGTTCGTCGTTGGTTACACGACGGCGTCGTTTAGGACGATGATGTTGATT",
                                          "TACGCGGGAGGGGTTGGTGTCACGACGTTGATCACGGTGCCGAACTGGCCATTCTTTAAC",
                                          "CGTCATCCACTCAAGTGGTTGGAACCAAGTGAAGCGGAGAAGCATCCTAAACCGGAAGTC",
                                          "GTTGTTAGCTCGAAGAAGAAGTCATCTAAAAAGTAG"])), id='test')

        # Run the design_primers function
        result = design.design_primers(sequence)
        # Check the result
        self.assertIsInstance(result, dict)

    def test_run(self):
        # Run the run function with a test sequence file and 1 process
        design.run('test_cds.fa', 'test.ini', 'test', 1, 1)
        # Check the result
        with open('test.json') as f:
            result = f.read()
        self.assertIsInstance(result, str)


if __name__ == '__main__':
    unittest.main()
