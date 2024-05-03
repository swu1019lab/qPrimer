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
pip install dist/qprimer-1.0.3.tar.gz --user
```

## Example

```shell
# Design primers
qPrimer --design --seq_file test_cds.fa

# Check primers specificity
qPrimer \
--design --check \
--seq_file test_cds.fa \
--database test.fa

# Annotate primers with gene information and SNP information
qPrimer \
--design --check --annotate \
--seq_file test_cds.fa \
--database test.fa \
--gtf_file genes.gtf \
--snp_file test_snps.bed

# Visualize primers results with a html page
# --visualize requires the output from --design, --check, and --annotate
qPrimer \
--design --check --annotate --visualize \
--seq_file test_cds.fa \
--database test.fa \
--gtf_file genes.gtf \
--snp_file test_snps.bed
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