# Gutmicro KITP project Summer 2021 day 2 tasks
# We're working with a (smaller) text file with the matrix of coverage depth after metagenomic alignment
# This file answers questions about the file such as:
# What fraction of the ref genome is not present in a typical sample. How does this vary across species/ species?
# Are some genes present in almost all samples?
# Can you identify the “typical”coverage for each sample? By eye? By some automatic method?
# What would you do to filter out sites not “present”in a given sample? 
# What tradeoffs are you thinking about re: the criteria you chose?

import bz2
import matplotlib.pyplot as plt
import numpy as np

# Read in the matrix
filename = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/snps/Bacteroides_vulgatus_57955/snps_depth_small.txt.bz2'
snps_file = bz2.BZ2File(filename,"r")
# What fraction of the ref genome is not present in a typical sample. How does this vary across species?
species_avg_depth = {}

# Calculate average depth across the samples (columns) for each position in the ref genome (rows)

line = str(snps_file.readline()) # header
items = line.strip().split("\\t")
samples = np.array([item.strip() for item in items[1:]])

for line in snps_file:
    #print(line)
    line = str(line)[:-3].split('\\t')    # Get rid of trailing \\n and separators \\t
    pos = int(line[0].split('|')[1])
    species = line[0].split('|')[0][2:]

    # Calculate average depth
    #print(line)
    depths_per_position = np.array([int(x) for x in line[1:]])
    avg = np.mean(depths_per_position)
    #print(depths_per_position)

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