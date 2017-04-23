##
# This script will extract all of the genes located in the completed clusters
# from clustering_pipeline.py. This is usefull for determinging what genes
# are not included in any of your clusters.
##

import sys, os

#print '''Assumes that all lines with genes start with AT. Only input is the 
#directory where the completed clusters are located.
#'''

completed = sys.argv[1] # directory with all of the completed clusters

try:
	gene_prefix = sys.argv[2]
except IndexError:
	gene_prefix = 'AT'

files = os.listdir(completed) # List of all files in the completed directory.

genes = set([]) # List of all genes in the clusters.

for cluster in files:
    
    for line in open(completed+'/'+cluster, 'r'):
        
        if line.startswith(gene_prefix):
            
            tab = line.strip().split('\t')
            
            name = tab[0] # Name of the gene in the 
            
            genes.add(name)

for i in genes:
    print i

