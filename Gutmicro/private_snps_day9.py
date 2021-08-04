# Gutmicro KITP project Summer 2021 day 9 tasks
# For the last analysis, we are focusing on the private SNPs each bacteria has in each host.
# Questions we will answer:
# Use the data in the private SNVs folder to get a sense of the typical number of private marker 
#   SNVs in a given person. Does it vary a lot across species and hosts?
# Main task: make a 2D scatter plot showing the “level of private marker SNV sharing” 
#   vs the number of genetic differences between time points for each of the HMP samples. 
#   We previously saw a big separation between “evolution”and “strain replacement” examples on 
#   the x-axis (# of genetic diffs). Is there also a big separation on the y-axis 
#   (level of private marker SNV sharing)?


from collections import Counter
import gzip
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.core.reshape.melt import wide_to_long


"""
# First, read in the private SNPs file to get a sense of the typical number of markers / person.
# Let's use Bacteroides vulgatus as an example.
# header: contig, location, gene_name, var_type, host_id, list_of_samples_where_snp_was_present
# Each row has a SNP that was present at high freq in ONLY a single host
file_path = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/private_snvs/Bacteroides_vulgatus_57955.txt.gz'

private_snps_hosts = []
with gzip.open(file_path, 'r') as f:
    f.readline()    # header
    for line in f:
        line = str(line)[:-3].split(', ')    # Get rid of trailing \\n and separators \\t
        host = line[4]
        private_snps_hosts.append(host)

# Plot a bar graph for the number of private SNPs for each host.
count_hosts = Counter(private_snps_hosts)

# Sort by number of private SNPs
l = sorted(count_hosts.items(), key = lambda kv:(kv[1], kv[0]),reverse=True)
names = [i[0] for i in l]
vals = [i[1] for i in l]

# Plot
plt.bar(names, vals)
plt.xlabel('Name of a host')
plt.xticks(rotation=90, size=3)
plt.ylabel('Number of private SNPs')
plt.title('How many private SNPs each host has for Bacteroides vulgatus?')
plt.show()
"""

# Next task, Distinguishing local evolution & strain replacement  
#   make a 2D scatter plot showing the “level of private marker SNV sharing” 
#   vs the number of genetic differences between time points for each of the HMP samples.

# header: contig, location, gene_name, var_type, host_id, list_of_samples_where_snp_was_present
file_path = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/private_snvs/Bacteroides_vulgatus_57955.txt.gz'

# Check which hosts are the HMP samples
# Subject/TwinID	SampleID	Study	Country	Continent	Timepoint/Twin
metadata_path = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/sample_metadata.txt'
meta_df = pd.read_csv(metadata_path, sep='\t', index_col='SampleID')

# Check how many private SNPs are shared between two time points for each host
file_path = '/home/ada/Desktop/KITP_tutorials/kitp_2021_microbiome_data/private_snvs/Bacteroides_vulgatus_57955.txt.gz'

host_timepoints = {}
with gzip.open(file_path, 'r') as f:
    f.readline()    # header
    for line in f:
        line = str(line)[3:-3].split(', ')    # Get rid of the trailing \\n 
        time_points = line[-1].split(';')
        host = line[-2]
        pos = int(line[1])

        # Check if host is from the HMP cohort
        row = meta_df.loc[meta_df['Subject/TwinID'] == host]
        if row['Study'][0] == 'hmp':

            # Make a map for each host. {host : {t1: [pos1, pos2], t2: [pos2, pos3]}}
            if host not in host_timepoints.keys():
                host_timepoints[host] = {}
                for i in time_points:
                    host_timepoints[host][i] = [pos]
            else:
                for i in time_points:
                    try:                # Some time points might not have been added yet.
                        host_timepoints[host][i].append(pos)
                    except KeyError:
                        host_timepoints[host][i] = [pos]

# Check the similarity between time points (for simplicity, we will consider just 2 points for each 
#   host)
within_host_similarity = {}
for k, v in host_timepoints.items():
    time_list = list(v.keys())
    set1 = set(v[time_list[0]])
    set2 = set(v[time_list[-1]])

    # Check how many items are in common
    common = len(set1.intersection(set2))
    all_pos = len(set1.union(set2))

    similarity = float(common)/all_pos
    within_host_similarity[k] = similarity

print(within_host_similarity)