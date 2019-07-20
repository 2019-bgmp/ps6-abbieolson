#!/usr/bin/env python

# don't forget to "conda activate bgmp_py3" on the command line

import re
import matplotlib.pyplot as plt

FILE = "/home/afo/bgmp/Bi621/PS6/KMER49/contigs_49_cov_cutoff_500.fa"
KMER_LEN = 49

# mean depth coverage for contigs, must solve for C
# Ck = C * (L - K + 1) / L (original equation)
# Ck = k-mer coverage
# C = coverage
# L = length of reads
# K = k-mer length
# C = Ck * L / (L - K + 1) (modified equation)

# C = kmer_cov * phys_len / (kmer_len - KMER_LEN + 1)

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

# 5. Calculate the N50 value of your assembly.

# If the position in the entire contig length (item), as it's incremented (total)
# is greater than the sum of the entire physical length divided by 2 (L50),
# then the N50 requirement is perfected.

phys_len_sort = sorted(phys_len)
L50 = sum(phys_len_sort) / 2
total = 0

for item in phys_len_sort:
    total += item
    if total >= L50:
        N50 = item
        break

# 6. Calculate the distribution of contig lengths and bucket the contig lengths into groups of
# 100bp. So, all contigs with lengths between 0 and 99 would be in the 0 bucket, those
# with lengths between 100 and 199 would be in the 100 bucket, etc.

# 7. Print out the distribution.
# # Contig length Number of contigs in this category
# 0 0
# 100 5324
# 200 3345
# 300 1130
# ...

# Must do floor division and multiply by 100 to make the contigs 100, 200, 300, etc.
# I initialized a dictionary (bucket) and incremented with an if, else. Ben also
# showed me how to use the .update function to overwrite a dictionary.

bucket = {}

keys = []
values = []

for item in phys_len_sort:
    x = (item // 100) * 100
    if x in bucket:
        bucket[x] += 1
    else:
        bucket[x] = 1

for key in bucket:
    keys.append(key)
    values.append(bucket[key])

print("# Contig length", "\t", "Number of contigs in this category")
for key in bucket:
    print(str(key), "\t", str(bucket[key]))

plt.bar(keys, values, width=100)
plt.xlabel("# Contig length")
plt.ylabel("Number of contigs in this category")
plt.yscale("log")
plt.title("Contig Distribution, K-mer size 49, cc auto, min contig 500 bp")
plt.savefig("contigs_49_cov_cutoff_500.png")
