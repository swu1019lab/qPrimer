## Description

The qPrimer package is a collection of tools for designing, checking, annotating, and visualizing qPCR primers.

![qPrimer_pipeline](https://cdn.jsdelivr.net/gh/swu1019lab/md_img/graphical_abstract-2.png)

## Contents
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
    - [design module](#design-module)
    - [annotate module](#annotate-module)
    - [check module](#check-module)
    - [visualize module](#visualize-module)
- [Output](#output)
- [Example](#example)
- [Citation](#citation)

## Dependencies
- [In-Silico PCR (isPcr)](): A tool for checking primer specificity
- Python 3.9+ and the following packages:
    - [biopython](https://biopython.org/): A set of tools for biological computation
    - [pandas](https://pandas.pydata.org/): A fast, powerful, flexible, and easy-to-use open-source data analysis and
      data manipulation library
    - [numpy](https://numpy.org/): The fundamental package for scientific computing with Python
    - [primer3-py](https://github.com/libnano/primer3-py): A Python wrapper for the Primer3 library
    - [jinja2](https://jinja.palletsprojects.com/): A fast, expressive, extensible templating engine

## Installation

```shell
# install the package from pip
pip install qPrimer
```

```shell
# install the package from source
git clone https://github.com/swu1019lab/qPrimer.git
cd qPrimer
# install the package using build and pip commands
pip install build --user
python -m build
pip install dist/qprimer-1.0.4.tar.gz --user
```

## Usage

```shell
# Display help information
qPrimer --help
```

```shell
usage: qPrimer [-h] [--log_file LOG_FILE] [--out_name OUT_NAME] {design,annotate,check,visualize} ...

Run qPrimer package.

positional arguments:
  {design,annotate,check,visualize}
                        Specify the module to run
    design              Run design module
    annotate            Run annotate module
    check               Run check module
    visualize           Run visualize module

optional arguments:
  -h, --help            show this help message and exit
  --log_file LOG_FILE   Path to the log file
  --out_name OUT_NAME   Prefix name of the output file
```

### design module

```shell
usage: qPrimer design [-h] --seq_file SEQ_FILE [--ini_file INI_FILE] [--processes PROCESSES] [--csv]

optional arguments:
  -h, --help            show this help message and exit
  --seq_file SEQ_FILE   Path to the sequence file
  --ini_file INI_FILE   Path to the ini file
  --processes PROCESSES
                        Number of CPUs to use
  --csv                 Output the results in csv format
```

default parameters setting (can ignore upper case and lower case) in the ini file:
```shell
[GLOBAL]
primer_task = generic
primer_num_return = 10
primer_pick_left_primer = 1
primer_pick_right_primer = 1
primer_min_gc = 40
primer_max_gc = 60
primer_opt_gc_percent = 50
primer_opt_size = 22
primer_min_size = 18
primer_max_size = 28
primer_opt_tm = 60
primer_min_tm = 58
primer_max_tm = 64
primer_pair_max_diff_tm = 3
primer_product_size_range = 80-300
primer_explain_flag = 1
primer_thermodynamic_oligo_alignment = 1
PRIMER_SECONDARY_STRUCTURE_ALIGNMENT = 1
primer_annealing_temp = 50
PRIMER_MAX_POLY_X = 5
PRIMER_MAX_SELF_ANY = 4
PRIMER_MAX_SELF_END = 3
PRIMER_PAIR_MAX_COMPL_ANY = 4
PRIMER_PAIR_MAX_COMPL_END = 3
PRIMER_MAX_END_STABILITY = 9

[SEQUENCE]

[PROGRAM]
```

parameters specification:

- primer_task: generic
- primer_num_return: 10, number of primer pairs to return
- primer_pick_left_primer: 1, whether to pick the left primer, 0 or 1
- primer_pick_right_primer: 1, whether to pick the right primer, 0 or 1
- primer_min_gc: 40, minimum GC content, 0-100
- primer_max_gc: 60, maximum GC content, 0-100
- primer_opt_gc_percent: 50, optimal GC content, 0-100
- primer_opt_size: 22, optimal primer length
- primer_min_size: 18, minimum primer length
- primer_max_size: 28, maximum primer length
- primer_opt_tm: 60, optimal melting temperature
- primer_min_tm: 58, minimum melting temperature
- primer_max_tm: 64, maximum melting temperature
- primer_pair_max_diff_tm: 3, maximum difference in melting temperature between the left and right primers
- primer_product_size_range: 80-300, range of product size, e.g., 80-300 or 100-200, 300-400
- primer_explain_flag: 1, whether to include the explanation of the primer
- primer_thermodynamic_oligo_alignment: 1, whether to use thermodynamic oligo alignment
- primer_secondary_structure_alignment: 1, whether to use secondary structure alignment
- primer_annealing_temp: 50, annealing temperature
- primer_max_poly_x: 5, maximum number of consecutive G/C/T/A in the primer
- primer_max_self_any: 4, maximum self-complementarity
- primer_max_self_end: 3, maximum self-complementarity at the 3' end
- primer_pair_max_compl_any: 4, maximum complementarity between the left and right primers
- primer_pair_max_compl_end: 3, maximum complementarity between the 3' end of the left primer and the 3' end of the
  right primer
- primer_max_end_stability: 9, maximum stability of the 3' end of the primer

parameters references:
>- Jeon H, Bae J, Hwang SH, et al. MRPrimerW2: an enhanced tool for rapid design of valid high-quality primers with
  multiple search modes for qPCR experiments. *Nucleic Acids Res*. 2019;47(W1):W614-W622. doi:10.1093/nar/gkz323
>- Lu K, Li T, He J, et al. qPrimerDB: a thermodynamics-based gene-specific qPCR primer database for 147 organisms.
  *Nucleic Acids Res*. 2018;46(D1):D1229-D1236. doi:10.1093/nar/gkx725

custom parameters setting for qPCR and other applications:
> All the parameters from Primer3 can be set in the ini file, and the default parameters are set for qPCR. The
> users can modify the parameters in the ini file to meet the specific requirements of the primer design not only for
> qPCR but also for other applications.

### annotate module

```shell
  
usage: qPrimer annotate [-h] --primers PRIMERS --gtf_file GTF_FILE --snp_file SNP_FILE

optional arguments:
  -h, --help            show this help message and exit
  --primers PRIMERS     Path to the primers file (json format)
  --gtf_file GTF_FILE   Path to the gtf file
  --snp_file SNP_FILE   Path to the snp file
  --processes PROCESSES
                        Number of CPUs to use
```

### check module

```shell
usage: qPrimer check [-h] --primers PRIMERS --database DATABASE

optional arguments:
  -h, --help            show this help message and exit
  --primers PRIMERS     Path to the primers file (json format)
  --database DATABASE   Path to the database file
  --processes PROCESSES
                        Number of CPUs to use
```

### visualize module

```shell
usage: qPrimer visualize [-h] --primers PRIMERS --seq_file SEQ_FILE --genes_num GENES_NUM

optional arguments:
  -h, --help            show this help message and exit
  --primers PRIMERS     Path to the primers file (json format)
  --seq_file SEQ_FILE   Path to the sequence file (fasta format)
  --genes_num GENES_NUM
                        Number of genes to display
```

## Output
The default output format is json, and the output file is named `qPrimer.json`. 
The output file contains the following three main sections:

1. Primer pair
   - **PENALTY**: The penalty score of the primer pair
   - **PRODUCT_SIZE**: The product size of the primer pair
   - **PRODUCT_TM**: The product melting temperature of the primer pair
   - **COMPL_ANY_TH**: The calculated value for the tendency of a primer pair to bind to each other
   - **COMPL_END_TH**: The calculated value for the tendency of the 3'-ENDs of a primer pair to bind to each other.

2. Left primer
   - **PENALTY**: The penalty score of the left primer
   - **SEQUENCE**: The sequence of the left primer
   - **COORDS**: The coordinates of the left primer
   - **TM**: The melting temperature of the left primer
   - **BOUND**: The fraction of primers bound at annealing temperature
   - **GC_PERCENT**: The GC content of the left primer
   - **END_STABILITY**: The stability of the 3' end of the left primer
   - **SELF_ANY_STUCT**: A string representation of the calculated secondary structure
   - **SELF_END_STUCT**: A string representation of the calculated secondary structure
   - **HAIRPIN_STUCT**: A string representation of the calculated secondary structure

3. Right primer
   - **PENALTY**: The penalty score of the right primer
   - **SEQUENCE**: The sequence of the right primer
   - **COORDS**: The coordinates of the right primer
   - **TM**: The melting temperature of the right primer
   - **BOUND**: The fraction of primers bound at annealing temperature
   - **GC_PERCENT**: The GC content of the right primer
   - **END_STABILITY**: The stability of the 3' end of the right primer
   - **SELF_ANY_STUCT**: A string representation of the calculated secondary structure
   - **SELF_END_STUCT**: A string representation of the calculated secondary structure
   - **HAIRPIN_STUCT**: A string representation of the calculated secondary structure

In addition to the json output, the results can also be output in csv format with the `--csv` option.

---
The calculation of the penalty score is based on the following formula:
> 
> penalty = 
> The penalty score is used to evaluate the quality of the primer pair, and the primer pair with the lowest penalty
> score is considered the best primer pair.

The calculation of the stability of the 3' end of the primer is based on the following formula:
>
> stability = 
> The stability score is used to evaluate the stability of the 3' end of the primer, and the primer with the lowest
> stability score is considered the best primer.

The calculation of the secondary structure of the primer is based on the following formula:
>

## Example

```shell
# Design primers with default json output
qPrimer design --seq_file test_cds.fa

# Design primers with json and csv output
qPrimer design --seq_file test_cds.fa --csv

# Design primers with multiple CPUs
qPrimer design --seq_file test_cds.fa --processes 4

# Design primers with specific parameters, like primer length, Tm, GC content, etc.
qPrimer design --seq_file test_cds.fa --ini_file test.ini

# Check primers specificity
qPrimer check \
--primers qPrimer.json \
--database test.fa

# Annotate primers with gene information and SNP information
qPrimer annotate \
--primers qPrimer.json \
--gtf_file genes.gtf \
--snp_file test_snps.bed

# Visualize primers results with a html page
qPrimer visualize \
--primers qPrimer.json \
--seq_file test_cds.fa \
--genes_num 10
```

![qPrimer_report-1](https://cdn.jsdelivr.net/gh/swu1019lab/md_img/qPrimer_report-1.jpeg)
![qPrimer_report-2](https://cdn.jsdelivr.net/gh/swu1019lab/md_img/qPrimer_report-2.jpeg)
![qPrimer_report-3](https://cdn.jsdelivr.net/gh/swu1019lab/md_img/qPrimer_report-3.jpeg)

## Citation

If you use qPrimer and qPrimerDB, please cite the following paper:

Kun Lu†, Tian Li†, Jian He†, Wei Chang†, Rui Zhang, Miao Liu, Mengna Yu, Yonghai Fan, Jinqi Ma, Wei Sun, Cunmin Qu,
Liezhao Liu, Nannan Li, Ying Liang, Rui Wang, Wei Qian, Zhanglin Tang, Xinfu Xu, Bo Lei, Kai Zhang*, Jiana Li*.
qPrimerDB: A thermodynamics-based gene-specific qPCR primer database for 147 organisms.
Nucleic Acids Research. 2018, 46: D1229-D1236.