# Gutmicro KITP project Summer 2021 day 7 tasks
# Now we are going to detect genetic differences between samples (Quasi-phasable samples)
# We calculate alt freq for each site for each population (fi1,fi2).
# Because of sampling noise that has a Binomial distribution around the true frequency
#   and because sampling noise is independent for each sample, we will get differences in the consensus 
#   that are false positives.
# We want to find sites where truely consensus 1 != consensus 2, but we can't expect full sweeps so
#   a threshold is chosen (80%)
# To distinguish between genotype 1 and 2, we check the freq for each site i. If fi > 80%, we assign it genotype 2
#   and when fi < 20 we assign it genotype 1 -> results for each bacteria in within_host_changes.txt 


# Use the within_host_changes.txt summary file to check the total numbers of changes
# what kind of numbers are right here, what numbers do you get out of looking at the number of genetic differences within people over time versus in twins 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Read in the within_host_changes table
# Cohort  Sample_t0       Sample_t1       Species NumTestedSNPs   NumSNPChanges   NumTestedGenes  NumGeneChanges
file_path = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/within_host_changes.txt'

df = pd.read_csv(file_path, sep='\t')

# Divide the data into differences between time points (same individual, different samples)
# Subject/TwinID	SampleID	Study	Country	Continent	Timepoint/Twin
metadata_path = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/sample_metadata.txt'
meta_df = pd.read_csv(metadata_path, sep='\t', index_col='SampleID')

# Check which rows in within_host changes are from the same person.
same_person_changes = []
for index, row in df.iterrows():
    s1, s2 = row.Sample_t0, row.Sample_t1
    if meta_df.loc[s1]['Subject/TwinID'] == meta_df.loc[s2]['Subject/TwinID']:
        same_person_changes.append(int(row.NumSNPChanges))


# What are the numbers of SNP changes  in the dataset for the same person at different time points?
# Make a logspace bins since we are interested in small changes in the genome (while a lot of changes are very big)
bins = np.logspace(np.log10(0.1),np.log10(6000.0), 100)
plt.hist(df.NumSNPChanges, bins=bins, edgecolor='black')
#plt.yscale('log')
plt.xscale('log')
plt.ylabel('Frequency')
plt.xlabel('Number of SNP changes (logspace bins)')
plt.title('Person-specific frequency of genetic changes within different time points compared in the file')
plt.show()
