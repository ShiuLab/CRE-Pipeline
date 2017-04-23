# CRE-Pipeline

#Scripts by Alex Seddon and Cheng Zou
# 1st step Clustering

install.packages('e1071')

python clustering_pipeline.py -m [expression matrix] --max [Maximum cluster size] --min [Minimum cluster size] -f ["k" for k-means, "c" for fuzzy c-means]
