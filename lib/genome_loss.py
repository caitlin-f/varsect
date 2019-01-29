""" Find stretches of low coverage. Can be used on the pipeline low coverage
output, however, minimum MAPQ is set at 30. Therefore, mpileup rerun with no
minimum MAPQ.
samtools mpileup -aa SAMPLE.bam | awk '($4 < 4)' > SAMPLE.mpileup"""

import sys

file = sys.argv[1]

regions = {}
prev_pos = None
start = None
current_chrom = ""
with open(file) as input:
    for line in input:
        chrom, position, char, cov, seq, det = line.strip().split("\t")
        pos = int(position)
        if regions.get(chrom) is None:
            regions[chrom] = []
        if start is None:
            start = pos
        elif current_chrom == chrom and pos > (prev_pos + 5): # coverage over 5bp length will 'break-up' deletion
            if prev_pos - start > 1000:
                regions[chrom].append((start, prev_pos))
            start = pos
            current_chrom = chrom
        elif current_chrom != chrom:
            if prev_pos - start > 1000:
                regions[current_chrom].append((start, prev_pos))
            start = pos
            current_chrom = chrom
        prev_pos = pos
    if prev_pos - start > 1000: # when read end of line
        regions[chrom].append((start, prev_pos))

print(regions)
