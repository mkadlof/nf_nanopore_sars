#!/usr/bin/env python3

"""
To prosty skrypt do generowania pokrycia w oparciu o quality threshold
"""

import sys

import pysam

plik_bam = sys.argv[1]
qual = int(sys.argv[2])
threshold = int(sys.argv[3])
sam_file = pysam.AlignmentFile(plik_bam, "rb")
# quality do ustalenia może na początku delikatne 20 ?
pokrycia = sam_file.count_coverage(contig='MN908947.3', start=0, stop=29903, quality_threshold=qual)
# pokrycia to 4 elementowa krotka z wartościami na danej pozycji A C T i G
with open('quality_mask.bed', 'w') as f:
    for i in range(len(pokrycia[0])):
        suma = pokrycia[0][i] + pokrycia[1][i] + pokrycia[2][i] + pokrycia[3][i]
        if suma < threshold:
            f.write(f'MN908947.3\t{i}\t{i + 1}\n')
