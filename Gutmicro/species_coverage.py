# Gutmicro KITP project Summer 2021
# Use metadata file and species module output to play around w/ species abundance distributions w/in and across hosts
# How many species have coverage >= 10x in average person?(repeat for 1x, 20x, ...)
# What genera do these species tend to come from?
# Input file: /home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/species/count_reads.txt.bz2

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

input_f = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/species/count_reads.txt.bz2'

# Create a dataframe from the reads counts, col = participants, row = bacterial species
coverage_df = pd.read_csv(input_f, sep='\t', index_col='species_id')

#How many species have coverage >= 10x in average person?(repeat for 1x, 20x, ...)
print('How many species have coverage >= 10x in average person? A: %d' %(sum(coverage_df.mean(axis=1)>=10)))   # Answer: 177

# 1x
print('How many species have coverage >= 1x in average person? A: %d' %(sum(coverage_df.mean(axis=1)>=1)))     # Answer: 315

# 20x
print('How many species have coverage >= 10x in average person? A: %d' %(sum(coverage_df.mean(axis=1)>=20)))   # Answer: 143