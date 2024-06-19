# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 9:36
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : visualize.py

import json
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from Bio import SeqIO
import os


def run(primers, seq_file, out_name, genes_num):
    print('Visualize module is running.')
    # Load the primer results from a JSON file
    # primer_results = json.load(open('tests/qPrimer.json'))
    primer_results = json.load(open(primers))

    # Create a Pandas DataFrame from the primer results
    # primer pairs (need to check if the columns are exist)
    columns = ['PENALTY', 'PRODUCT_SIZE', 'PRODUCT_TM', 'EXON_SPAN', 'SNP_SPAN', 'SPECIFICITY']
    columns = list(set(columns).intersection(set(primer_results[0]['PRIMER_PAIR'][0].keys())))
    df_p = pd.concat([pd.DataFrame.from_records(res['PRIMER_PAIR'], columns=columns) for res in primer_results])
    all_primers_num = df_p.shape[0]
    best_primers_num = df_p.query('SPECIFICITY == 1').shape[0] if 'SPECIFICITY' in df_p.columns else 0

    # left and right primers (need to check if the columns are exist)
    columns = ['PENALTY', 'TM', 'BOUND', 'GC_PERCENT', 'END_STABILITY', 'EXON_SPAN', 'SNP_SPAN']
    columns = list(set(columns).intersection(set(primer_results[0]['PRIMER_LEFT'][0].keys())))
    df_f = pd.concat([pd.DataFrame.from_records(res['PRIMER_LEFT'], columns=columns) for res in primer_results])
    df_r = pd.concat([pd.DataFrame.from_records(res['PRIMER_RIGHT'], columns=columns) for res in primer_results])
    if 'EXON_SPAN' not in columns:
        df_f['EXON_SPAN'], df_r['EXON_SPAN'] = 0, 0
    if 'SNP_SPAN' not in columns:
        df_f['SNP_SPAN'], df_r['SNP_SPAN'] = 0, 0

    # Parse the sequence file
    # sequences = SeqIO.parse('tests/test_cds.fa', 'fasta')
    sequences = SeqIO.parse(seq_file, 'fasta')
    seq_num = len(list(sequences))

    # Summary of the primer results
    exon_span = pd.DataFrame.from_records([df_p['EXON_SPAN'].value_counts().to_dict(),
                                           df_f['EXON_SPAN'].value_counts().to_dict(),
                                           df_r['EXON_SPAN'].value_counts().to_dict()])
    exon_span.index = ['Primer pairs', 'Left primers', 'Right primers']
    # print(json.dumps(exon_span.T.reset_index().to_dict(orient='records')))

    snp_span = pd.DataFrame.from_records([df_p['SNP_SPAN'].value_counts().to_dict(),
                                          df_f['SNP_SPAN'].value_counts().to_dict(),
                                          df_r['SNP_SPAN'].value_counts().to_dict()])
    snp_span.index = ['Primer pairs', 'Left primers', 'Right primers']
    # print(json.dumps(snp_span.T.reset_index().to_dict(orient='records')))

    tm_stat_data = [
        ['Type', 'TM'],
        ['Primer pairs', df_p['PRODUCT_TM'].mean()],
        ['Left primers', df_f['TM'].mean()],
        ['Right primers', df_r['TM'].mean()]
    ]

    penalty_stat_data = [
        ['Type', 'Penalty'],
        ['Left primers', df_f['PENALTY'].mean()],
        ['Right primers', df_r['PENALTY'].mean()]
    ]

    bound_stat_data = [
        ['Type', 'Bound'],
        ['Left primers', df_f['BOUND'].mean()],
        ['Right primers', df_r['BOUND'].mean()]
    ]

    gc_stat_data = [
        ['Type', 'GC Percent'],
        ['Left primers', df_f['GC_PERCENT'].mean()],
        ['Right primers', df_r['GC_PERCENT'].mean()]
    ]

    end_stability_stat_data = [
        ['Type', 'End Stability'],
        ['Left primers', df_f['END_STABILITY'].mean()],
        ['Right primers', df_r['END_STABILITY'].mean()]
    ]

    design_summary = {
        'gene_num': len(primer_results),
        'all_primers_num': all_primers_num,
        'best_primers_num': best_primers_num,
        'specificity_rate': round(best_primers_num / all_primers_num * 100, 2),
        'primer_pairs_per_gene': round(all_primers_num / seq_num, 2),
        'success_rate': round(len(primer_results) / seq_num * 100, 2),
        'primer_pair': df_p.to_dict(orient='records'),
        'exon_span': json.dumps(exon_span.T.sort_index().reset_index().to_dict(orient='records')),
        'snp_span': json.dumps(snp_span.T.sort_index().reset_index().to_dict(orient='records')),
        'tm': json.dumps(tm_stat_data),
        'penalty': json.dumps(penalty_stat_data),
        'bound': json.dumps(bound_stat_data),
        'gc': json.dumps(gc_stat_data),
        'end_stability': json.dumps(end_stability_stat_data)
    }

    # Create a Jinja2 environment
    index_dir = os.path.join(os.path.dirname(__file__))
    env = Environment(loader=FileSystemLoader(index_dir))
    template = env.get_template('index.html')

    # Render the template with the primer results
    html = template.render(data=primer_results[:genes_num], summary=design_summary)

    # Save the rendered HTML to a file
    with open(out_name + ".html", 'w', encoding="utf-8") as file:
        file.write(html)

    print('Visualize module is finished.')
