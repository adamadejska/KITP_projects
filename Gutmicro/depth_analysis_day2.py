# Gutmicro KITP project Summer 2021 day 2 tasks
# We're working with a (smaller) text file with the matrix of coverage depth after metagenomic alignment
# This file answers questions about the file such as:
# What fraction of the ref genome is not present in a typical sample. How does this vary across species?
# Are some genes present in almost all samples?
# Can you identify the “typical”coverage for each sample? By eye? By some automatic method?
# What would you do to filter out sites not “present”in a given sample? 
# What tradeoffs are you thinking about re: the criteria you chose?

import bz2
import gzip
import matplotlib.pyplot as plt
import numpy as np

# Read in the matrix
filename = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/snps/Bacteroides_vulgatus_57955/snps_depth_small.txt.bz2'
snps_file = bz2.BZ2File(filename,"r")

"""
# What fraction of the ref genome is not present in a typical sample. How does this vary across species?
species_avg_depth = {}

# Calculate average depth across the samples (columns) for each position in the ref genome (rows)

line = str(snps_file.readline()) # header
items = line.strip().split("\\t")
samples = np.array([item.strip() for item in items[1:]])

for line in snps_file:
    line = str(line)[:-3].split('\\t')    # Get rid of trailing \\n and separators \\t
    pos = int(line[0].split('|')[1])
    species = line[0].split('|')[0][2:]

    # Calculate average depth
    depths_per_position = np.array([int(x) for x in line[1:]])
    avg = np.mean(depths_per_position)

    # Add all the info to dict
    if species not in species_avg_depth.keys():
        species_avg_depth[species] = [avg]
    else:
        species_avg_depth[species].append(avg)

# Create a graph of avg depth across the genome for each species.
for s, avg in species_avg_depth.items():
    title = 'Species ' + s + ' avg depth across the genome'
    plt.plot(range(0,len(avg)), avg)
    plt.title(title)
    plt.xlabel('genome position')
    plt.ylabel('average depth')
    plt.tick_params(labelbottom=False)
    plt.show()
"""

# Are some genes present in almost all samples?
# To answer this question, we need to use a MIDAS database file rep_genomes/Bacteroides_vulgatus_57955/genome.features.gz
# Need some kind of 'typical' coverage that would serve as a threshold to indicate if the site is really there or not.

# Set a threshold (need a formal calculations)
typical_threshold = 20

# Read in the depth matrix and look for sites that are above the threshold for 'almost all' samples
almost_all = 0.8

# Read the matrix file.
line = str(snps_file.readline()) # header
items = line.strip().split("\\t")
samples = np.array([item.strip() for item in items[1:]])

pos_in_all = []
for line in snps_file:
    line = str(line)[:-3].split('\\t')    # Get rid of trailing \\n and separators \\t
    pos = int(line[0].split('|')[1])
    species = line[0].split('|')[0][2:]

    # Calculate average depth
    depths_per_position = np.array([int(x) for x in line[1:]])

    # Check if position has a coverage greater than threshold for almost all samples.
    cov_more_than_threshold = depths_per_position >= typical_threshold
    if sum(cov_more_than_threshold) / float(len(samples)) >= almost_all:
        pos_in_all.append(pos)

print(len(pos_in_all))

# Now that we have positions that are present in almost all samples, look up what they code in the genome features file.
# gene_id scaffold_id     start   end     strand  gene_type       functions
features_file = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/rep_genomes/Bacteroides_vulgatus_57955/genome.features.gz'
features = gzip.open(features_file, 'rb')
header = features.readline()

genes_in_all = set()
# Read the file, get the start and end positions of the gene and check if the positions we identified previously
#   fall inside any of the genes. 
for line in features:
    line = str(line).split('\\t')
    gene_id = line[0][2:]
    start, end = int(line[2]), int(line[3]) 
    for i in pos_in_all:
        if i >= start and i <= end:
            genes_in_all.add(gene_id)
        elif i > end or i < start:           # to help with run time
            continue


print(len(genes_in_all))