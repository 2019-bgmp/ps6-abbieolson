#!/usr/bin/env python

# ACTUAL RESULTS FROM UNIT TEST FILE

import re

FILE = "/Users/abbiefayeolson/bgmp/bi621/ps6-abbieolson/Unit_test.fa"
KMER_LEN = 49

with open(FILE, "r") as fh:

    phys_len = []
    kmer_cov = []
    contig_num = 0

    for line in fh:
        line = line.strip() # python adds a new line and needs to be stripped

        if line[0] == '>': # if the beginning of the line starts with a '>'
            phys_len.append(re.findall("length_[0-9]+", line)[0][7:]) # findall grabs the specified line but needs to have more grabbed from it, hence specifying the first 7 characters, append to empty array

            kmer_cov.append(re.findall("cov_[0-9]+\.[0-9]+", line)[0][4:])
            for cov in range(len(kmer_cov)):
                kmer_cov[cov] = float(kmer_cov[cov])

    for index, value in enumerate(phys_len):
        phys_len[index] = int(value) + KMER_LEN - 1 # the physical length is "strlen" (ie kmer_len) in the formula "kcnt = strlen" - k + 1"

contig_num += 1 # increment for contig length

max_contig = max(phys_len)
mean_contig = sum(phys_len) / len(phys_len)
mean_cov_depth = sum(kmer_cov) / len(kmer_cov)
total_phys = sum(phys_len)
num_contigs = len(phys_len)

# print(max_contig)
# output: 50

# print(mean_contig)
# output: 30.0

# print(mean_cov_depth)
# output: 385.04105366666664

# print(total_phys)
# output: 240

# print(num_contigs)
# output: 3
