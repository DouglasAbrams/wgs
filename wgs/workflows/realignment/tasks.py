'''
Created on Aug, 2019

@author: Diljot Grewal
'''

import os

import pypeliner
import pysam
import time
from wgs.utils import helpers


def split_by_rg(infile, read1_output, read2_output, tempdir):
    helpers.makedirs(tempdir)

    print("***********")
    print(tempdir)
    print(os.listdir(tempdir))
    print("***********")

    cmd = ['wgs_bamtofastq', infile, tempdir]
    pypeliner.commandline.execute(*cmd)

    print("***********")
    print(tempdir)
    print(os.listdir(tempdir))
    print("***********")

    try:
        readgroups = os.listdir(tempdir)
    except OSError:
        time.sleep(60)
        readgroups = os.listdir(tempdir)

    for readgroup in readgroups:
        os.rename(
            os.path.join(tempdir, readgroup, 'R1.fastq.gz'),
            read1_output[readgroup]
        )

        os.rename(
            os.path.join(tempdir, readgroup, 'R2.fastq.gz'),
            read2_output[readgroup]
        )


def get_read_group(infile):
    bam = pysam.AlignmentFile(infile, mode='rb', check_sq=False)
    header = bam.header['RG']
    return {rg['ID']: rg for rg in header}
