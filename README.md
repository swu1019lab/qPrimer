## Description
The qPrimer package is a collection of tools for designing, checking, annotating, and visualizing qPCR primers.

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
# --visualize requires the output from --design, --check, and --annotate
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