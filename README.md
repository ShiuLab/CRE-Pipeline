# CRE-Pipeline

#Scripts by Alex Seddon and Cheng Zou
# Coexpression clustering

install.packages('e1071')

python clustering_pipeline.py -m [expression matrix] --max [Maximum cluster size] --min [Minimum cluster size] -f ["k" for k-means, "c" for fuzzy c-means]

python get_genes_all_clusters.py [Directory of completed clusters] > [Desired output name]

python cluster_enrichment.py -m [Directory of raw expresion files] -f [File containing a list of files to use in the expression directories] -c [location of directory with completed clusters]

python Test_Fisher.py [5col file] 1 

python FastaManager.py -f getseq2 -fasta [FASTA contiaing all of the promoter regions] -name [List of gene names in your cluster]

# Motif finding

Running on single cluster

Running on multiple clusters

# Convert to TAMO format and concatenating files

module load TAMO

MDscan

If your files do not end with -n-md.out:

for i in {4..15}; do mv epi_lat.MDscan.w$i epi_lat.MDscan-$i-md.out; done

else:

python MDscan_tamo.py [Directory where MDscan output is located]

python MDscan_cat.py [Directory of the TAMO formated files] [Common TAMO folder]

MEME

python MEME_tamo.py [Full path to MEME files directory]

python MEME_cat.py [Full path to TAMO converted MEME] [Common TAMO folder]

MotifSampler

python MS_tamo.py [Directory of Motifsampler files]

python MS_cat.py [Directory of TAMO converted Motifsampler files] [Common TAMO folder]

AlignAce

python AA_tamo.py [Directory of where AlignACE files are located]

python AA_cat.py [Directory of TAMO converted AlignACE files] [Common TAMO folder] [Prefix of the AlignACE files]

YMF

python YMF_sumup_and_cat.py [Directory of YMF files] [Common TAMO folder] [Prefix for all of the YMF files in the YMF file directory]

Weeder

python weeder_sumup_and_cat.py [Directory of weeder files] [Common TAMO folder]

# KL clustering of motifs and UPGMA

# Mapping and enrichment

python 1_TAMO_split.py [TAMO from KL distance merging] 190

python 2_create_mapping_cc.py [Path to directory with split TAMO file] [location of the FASTA file containing sequence to be mapped]

nohup python qsub_hpc.py -f queue -c run_mapping_cc -w 90 -n 245 -u [User name] -m 2 &

python 4_motif_enrichment_5col.py -d [Directory with mapping files] -g [List of genes of interest] -s [FASTA file used for mapping] -o [Name of output] 

python Test_Fisher.py [5col file] 1

python 5_get_significant_motifs.py -in [Results of the Fisher test] -a [P-value threshold, 0.05 by default] -map_dir [Full path to the directory of the mapping files] -tamo [TAMO file used for mapping]

# PCC distance merging

# Location merging
