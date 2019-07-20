#!/usr/bin/env python

FQ_1 = open("/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1", "r")
FQ_2 = open("/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2", "r")
FQ_UNMATCHED = open("/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched", "r")

seq_cnt = 0
LN = 0
total_nt = 0
gene_len = (40 * 50) # 50 fosmids, each (approximately) 40 Kb long


# When the sequence count line is reached, increment into LN counter. Increment
# a counter for the sequence count (seq_cnt) and make sure to strip the line (line.strip()), as the
# new line default is white space. Increment again into the total nucleotide counter by
# taking the count of the sequence lines.

for line in FQ_1:
    LN += 1
    if LN%4==2:
        seq_cnt += 1
        line = line.strip()
        total_nt += len(line)

for line in FQ_2:
    LN += 1
    if LN%4==2:
        seq_cnt += 1
        line = line.strip()
        total_nt += len(line)

for line in FQ_UNMATCHED:
    LN += 1
    if LN%4==2:
        seq_cnt += 1
        line = line.strip()
        total_nt += len(line)

mean_read_len = total_nt / seq_cnt

expect_cov = total_nt / gene_len


# K-mer coverage:

# Ck = C * (L - K + 1) / L (original equation)
# Ck = k-mer coverage
# C = coverage
# L = length of reads
# K = k-mer length
