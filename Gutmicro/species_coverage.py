# Gutmicro KITP project Summer 2021
# Use metadata file and species module output to play around w/ species abundance distributions w/in and across hosts
# How many species have coverage >= 10x in average person?(repeat for 1x, 20x, ...)
# What genera do these species tend to come from?
# Input file: /home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/species/count_reads.txt.bz2

from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import seaborn as sns
from taxonomy_ranks import TaxonomyRanks


input_f = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/species/count_reads.txt.bz2'

# Create a dataframe from the reads counts, col = participants, row = bacterial species
coverage_df = pd.read_csv(input_f, sep='\t', index_col='species_id')
print(coverage_df.head())
# Get a general idea of what the value distribution looks like.
#fig, ax = plt.subplots(figsize=(20,20))  
#sns.heatmap(coverage_df, norm=LogNorm(), ax=ax)
#plt.show()

#How many species have coverage >= 10x in average person?(repeat for 1x, 20x, ...)
print('How many species have coverage >= 10x in average person? A: %d' %(sum(coverage_df.mean(axis=1)>=10)))   # Answer: 177

# 1x
print('How many species have coverage >= 1x in average person? A: %d' %(sum(coverage_df.mean(axis=1)>=1)))     # Answer: 315

# 20x
print('How many species have coverage >= 20x in average person? A: %d' %(sum(coverage_df.mean(axis=1)>=20)))   # Answer: 143


# What genera do these species tend to come from? Try on the mean coverage >= 10.
bacteria_cov_10 = coverage_df.index[coverage_df.mean(axis=1)>=10]

cleaned_names = []
# Clean the names
for i in bacteria_cov_10:
    name = i.split('_')[:2]
    name = ' '.join(name)
    cleaned_names.append(name)

errors = []
genera = []
for i in cleaned_names:
    try: 
        rank_taxon = TaxonomyRanks(i)
        rank_taxon.get_lineage_taxids_and_taxanames()
        taxid = list(rank_taxon.lineages.keys())[0]
        genus = rank_taxon.lineages[taxid]['genus'][0]
        genera.append(genus)
    except:
        # Some names in the database have misspellings (catch those errors)
        errors.append(i)

# Show the pie chart of the genera
f = Counter(genera)
labels = f.keys()
vals = f.values()
plt.pie(vals, labels=labels, textprops={'fontsize': 8})
plt.title('Genera names the species tend to come from when avg coverage >= 10')
#plt.show()

# How many of these species does an average pair of people share? What about samples from the same person over time?