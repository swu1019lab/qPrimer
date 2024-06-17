# -*- coding: utf-8 -*-
# @Time    : 2024/3/12 9:26
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : setup.py

from setuptools import setup, find_packages
import os
import stat

# Get a list of all files in the 'qPrimer/bin' directory
bin_files = [os.path.join('qPrimer', 'bin', file) for file in os.listdir('qPrimer/bin')]

# Change the permission of all files in the 'qPrimer/bin' directory to 755
for file in bin_files:
    os.chmod(file, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)

setup(
    name='qPrimer',
    version='1.0.4',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'qPrimer = qPrimer.qPrimer:main'
        ]
    },
    data_files=[('qPrimer/bin', bin_files), ('qPrimer', ['qPrimer/index.html'])],
)
