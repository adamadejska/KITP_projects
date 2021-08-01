# Gutmicro KITP project Summer 2021 day 4 tasks
# We are working with a subset of SNP data: annotated_snps_small.txt.bz2 
# We are trying to describe the distribution of SNPs in each sample.
# First, we plot within-host SFSs for the Bacteroides and Alistipes samples and answer
#   questions about Can you pick out some examples of single vs multi-colonization by eye?
#   Do you notice any differences between the two species?
# After that, we will discuss what makes the sample simple vs complex and what is the 
# ratio pN/pS and what does it mean?

import bz2
import matplotlib.pyplot as plt
import numpy as np


def make_sfs(filename, samples_ids, name):
    #This function makes Site-Frequency-Spectrum graphs given a SNP file and id numbers of samples 

    # Read the file
    for i in samples_ids:

        # The SNPs file is ~120 Mb big, might be better to read it without saving the whole matrix into the memory. 
        snps_file = bz2.BZ2File(filename,"r")

        # First task: plot within-host Site-Frequency-Spectrum for each type of bacteria samples
        line = str(snps_file.readline()) # header
        items = line.strip().split("\\t")
        samples = np.array([item.strip() for item in items[1:]])

        mut_freq = []
        for line in snps_file:
            # Each value in the matrix is [Ai, Di] -> Ai = number of alt alleles at site i sample j, Di = depth for site i sample j
            line = str(line)[:-3].split('\\t')    # Get rid of trailing \\n and separators \\t
            info = line[0][2:]                    # ID | site i | gene ID | type of mutation (1D:nonsynonymous, 4D: synonymous) | ? | ?
            
            alt_num = float(line[i].split(',')[0])
            depth = int(line[i].split(',')[1])
            if depth == 0 or alt_num == 0:
                continue

            # Save frequencies
            mut_freq.append(alt_num/depth)

        # Plot SFSs for each indicated sample (ignore values above 0.5 since those aren't minority alleles)
        bins = np.arange(0, 0.5, 0.01)
        plt.hist(mut_freq, bins=bins, edgecolor='black')
        plt.title('SFS for sample ' + samples[i-1] + ' in ' + name + ' file')
        plt.xlabel('Alt frequency (alt num / depth) (bin size=0.01)')
        plt.ylabel('Frequency (log scale)')
        plt.yscale('log')
        plt.show()

    return 

"""
# First, we make SFS graphs for Alistipes putredinis
# Read in the matrix
filename = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/snps/Alistipes_putredinis_61533/annotated_snps_small.txt.bz2'
samples_ids = range(1, 39, 1)      # Samples we will make the SFSs for.
name = 'Alistipes_putredinis_61533'
make_sfs(filename, samples_ids)
"""

# Next, make SFS graphs for the second bacteria  Bacteroides vulgatus
filename = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/snps/Bacteroides_vulgatus_57955/annotated_snps_small.txt.bz2'
samples_ids = range(1, 39, 1)      # Samples we will make the SFSs for.
name = 'Bacteroides_vulgatus_57955'
make_sfs(filename, samples_ids, name)