import sys
import numpy as np

chrNames = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
]

# how long we consider each chromosome (a rough upper bound)
chrLengths = [
    250000000,244000000,199000000,192000000,181000000,172000000,160000000,147000000,
    142000000,136000000,136000000,134000000,116000000,108000000,103000000,91000000,
    82000000,79000000,60000000,64000000,49000000,52000000,156000000,60000000
]

# compute offsets to embed al the chromsomes in a linear sequence
chrOffsets = dict()
for i in xrange(len(chrLengths)):
    chrOffsets[chrNames[i]] = sum(chrLengths[0:i])

def get_bin_pos(chr, pos):
    global chrOffsets
    return int((chrOffsets[chr]+pos)/1000)

# mark all binned that are touched with 1
binValues = np.zeros(int(sum(chrLengths)/1000), dtype=np.int32)
for line in sys.stdin:
    parts = line.split()
    if parts[0] in chrOffsets:
        startPos = get_bin_pos(parts[0], int(parts[1]))
        endPos = get_bin_pos(parts[0], int(parts[2]))
        for i in xrange(startPos, endPos+1):
            binValues[i] = 1

# print the output of all binary bins
np.savetxt(sys.stdout, binValues, fmt="%d")