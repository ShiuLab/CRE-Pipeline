# pCRE-Pipeline

#Scripts by Alex Seddon and Cheng Zou

This protocol will show you how to run the motif finders step by step. Motif finding is a two step process. The first step is to use six motif finding packages to identify novel motifs in the promoter regions of your clusters. The second step is to combine the motifs that have been identified using TAMO to eliminate redundant motifs identified by the motif finders. Alex Seddon

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

# Convert to TAMO format and concatanate files

module load TAMO

MDscan

If your files do not end with -n-md.out: for i in {4..15}; do mv epi_lat.MDscan.w$i epi_lat.MDscan-$i-md.out; done

otherwise use the following

python MDscan_tamo.py [Directory where MDscan output is located]

(using 1.3_MD_scan.py)

python MDscan_cat.py [Directory of the TAMO formated files] [Common TAMO folder]

MEME

python MEME_tamo.py [Full path to MEME files directory]

(using 1.2_MEME2tamo.py and 2_split_length.py)

python MEME_cat.py [Full path to TAMO converted MEME] [Common TAMO folder]

MotifSampler

python MS_tamo.py [Directory of Motifsampler files]

(using 1.4_MotifSampler.py)

python MS_cat.py [Directory of TAMO converted Motifsampler files] [Common TAMO folder]

AlignAce

python AA_tamo.py [Directory of where AlignACE files are located]

(using 1.1_AA2tamo.py)

python AA_cat.py [Directory of TAMO converted AlignACE files] [Common TAMO folder] [Prefix of the AlignACE files]

YMF

python YMF_sumup_and_cat.py [Directory of YMF files] [Common TAMO folder] [Prefix for all of the YMF files in the YMF file directory]

(using 2.transform_tm.py)

Weeder

python weeder_sumup_and_cat.py [Directory of weeder files] [Common TAMO folder]

(using 2.transform_tm.py)

# KL clustering of motifs and UPGMA

-Create one TAMO file with all of the motifs (i.e. a TAMO that is not divided up by size)

-python 5.1_motif_clustering_step1.py -m all.tm -d 0.03964 --dfunc entropyrange to cluster very similar motifs (i.e. Any motif with the same one-letter code is clustered and merged. This reduced the motifs from ~80000 to ~6000)

-Combine all the merged TAMOs back together using cat *.0.03964-cl > After_cl_all.tm_d0.03964. 

-Run the UPGMA script to merge all the motifs in the larger TAMO: python UPGMA.py -m /mnt/research/ShiuLab/6_coronatine_CRE/5_motif_merging/After_cl_all.tm_d0.03964 -d 0.03964 --dfunc entropyrange >log 
This will take a long time, so I would recomend submitting it as a job.

# Mapping and enrichment

python 1_TAMO_split.py [TAMO from KL distance merging] 190

python 2_create_mapping_cc.py [Path to directory with split TAMO file] [location of the FASTA file containing sequence to be mapped]

nohup python qsub_hpc.py -f queue -c run_mapping_cc -w 90 -n 245 -u [User name] -m 2 &

python 4_motif_enrichment_5col.py -d [Directory with mapping files] -g [List of genes of interest] -s [FASTA file used for mapping] -o [Name of output] 

python Test_Fisher.py [5col file] 1

python 5_get_significant_motifs.py -in [Results of the Fisher test] -a [P-value threshold, 0.05 by default] -map_dir [Full path to the directory of the mapping files] -tamo [TAMO file used for mapping]

# PCC distance merging

# Location merging
