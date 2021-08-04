# Gutmicro KITP project Summer 2021 day 9 tasks
# For the last analysis, we are focusing on the private SNPs each bacteria has in each host.
# Questions we will answer:
# Use the data in the private SNVs folder to get a sense of the typical number of private marker 
#   SNVs in a given person. Does it vary a lot across species and hosts?
# Main task: makea 2D scatter plot showing the “level of private marker SNV sharing” 
#   vs the number of genetic differences between timepointsfor each of the HMP samples. 
#   We previously saw a big separation between “evolution”and “strain replacement”examples on 
#   the x-axis (# of genetic diffs). Is there also a big separation on the y-axis 
#   (level of private marker SNV sharing)?


from collections import Counter
import gzip
import matplotlib.pyplot as plt
import numpy as np


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