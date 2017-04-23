# This script is used to create clusters that are less than a certain size. For
# example, you may want to cluster an expression matrix, and you want clusters
# that are larger than 9 genes, but smaller than 60 genes. This script will 
# Create initial clusters. The clusters that are bellow the minimum size of 9
# genes will be thrown out, and then the clusters larger than 60 genes will be
# clustered recursively to smaller clusters.

import sys, os

function = 'kmeans_clustering.R'

print '''~ Parameters
-m    := The initial matrix for clustering. Assumes that it has row and column
         labels.
--max := The maximum cluster size.
--min := The minimum cluster size.
-f    := Specify if you want k-means (k) or fuzzy c-means (c) clustering.
         K-means is the default.
'''

# Get the directory where the script is located.
script_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])

# Get the parameters
for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-m':
        matrix      = sys.argv[i+1]
    if sys.argv[i] == '--max':
        max_size    = int(sys.argv[i+1])
    if sys.argv[i] == '--min':
        min_size    = int(sys.argv[i+1]) - 1 
    if sys.argv[i] == '-f':
        function    = sys.argv[i+1]
        if   function == 'k':
            function = 'kmeans_clustering.R'
        elif function == 'c':
            function = 'fuzzy_cmeans_clustering.R'
        else:
            pass

print function
function = script_dir + '/' + function

# Do the initial clustering. The last parameter is y, we need to cange the 
# rownames of this matrix in R.
os.system("R --vanilla --slave --args %s 100 y < %s" % (matrix, function))

# This function figures out the size of the clusters using wc.
def cluster_size():
    # Get the cluster sizes and loaded them in to the fDict
    os.system('wc * > tmp')
    wc = open('tmp', 'r')
    
    lines = []
    for l in wc:
        lines.append(l.strip().split(' '))
    
    fDict = {} # The dictionary contianing the files and the motifs in the files.
    for i in lines:
        
        x = []
        
        for n in i:
            if not n == '':
                x.append(n)
        
        fDict[x[3].strip()] = int(x[0])-1 #file name as a key, linked to cluster size
    
    return(fDict) 
        
# Create directories for the completed clusters and for the thrown out clusters.
os.system('mkdir completed thrown_out reclustered')

reclustered = 10 # This is an arbitrary initial value.

# This while loop will run untill no more reclustering is performed.
while reclustered > 0:

    # Get the cluster size dictionary.
    fDict = cluster_size()
    
    # Reset the reclustering counter
    reclustered = 0
    
    # For each of the cluster files, check there size, and move them accordingly
    for cluster in fDict:
        # If the cluster is the same as the orginal matrix, skip it.
        if cluster == matrix.split('/')[-1]:
            print "skipped %s" % cluster
            
        # Pass if the cluster is actually a directory or some other type of file,
        # other than a cluster.
        elif not cluster.startswith(matrix.split('/')[-1]):  
            print "skipped %s" % cluster
            
        # Through out the cluster if it is too small
        elif fDict[cluster] < min_size:
            os.system("mv %s thrown_out" % cluster)
            print "%s thrown out because it has %s genes" % (cluster, fDict[cluster])
            
        # Keep the cluster if it is within the minimum and maximum size.
        elif fDict[cluster] < max_size:
            os.system("mv %s completed" % cluster)
            print "%s is complete because it has %s genes" % (cluster, fDict[cluster])
            
        # If the cluster is too large, recluster it.
        elif fDict[cluster] >= max_size:
            os.system("mv %s reclustered" % cluster)
            # 
            os.system("R --vanilla --slave --args reclustered/%s 60 y < %s" % (cluster, function))
            os.system("mv reclustered/%s.* ./" % cluster)
            print "%s is reclustered because it has %s genes" % (cluster, fDict[cluster])
            reclustered += 1

os.system("rm tmp")          
print "Done!"
