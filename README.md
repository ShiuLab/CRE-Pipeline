# CRE-Pipeline

#Scripts by Alex Seddon and Cheng Zou
# 1st step Clustering

install.packages('e1071')

python clustering_pipeline.py -m [expression matrix] --max [Maximum cluster size] --min [Minimum cluster size] -f ["k" for k-means, "c" for fuzzy c-means]

python get_genes_all_clusters.py [Directory of completed clusters] > [Desired output name]

python cluster_enrichment.py -m [Directory of raw expresion files] -f [File containing a list of files to use in the expression directories] -c [location of directory with completed clusters]

python Test_Fisher.py [5col file] 1 

python FastaManager.py -f getseq2 -fasta [FASTA contiaing all of the promoter regions] -name [List of gene names in your cluster]

# 2nd step Motif Finders

Running on single cluster
Running on multiple clusters

# Convert to TAMO format

module load TAMO

MDscan, Your files need to end with -n-md.out if NOT
for i in {4..15}; do mv epi_lat.MDscan.w$i epi_lat.MDscan-$i-md.out; done

