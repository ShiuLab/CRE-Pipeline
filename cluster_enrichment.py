##
# This script is used to create a 5col file used for checking the enrichment of 
# a specific stress condition within a cluster.
# Assumptions:
# The column names in the expression matricies are the same as the raw expression
# file name, except that the file should start with contrast_.
##

################################################################################
# Define all functions

# This function takes the file of focus conditions, and converts it into a list.
def load_conditions(cond):

    conditions = [] # List of focus conditions 
    
    # This loop takes each condition in the file and puts it in the list.
    for line in open(cond, 'r'):
        conditions.append(line.strip())
    
    return(conditions)

# This function will create two sets. The first is the set of all genes, and the
# second is the set of all genes considered differentially expressed based on 
# its fold change and p-value.
def get_reg_sets(exprs, conditons):
    
    # remove trailing / from the expression directory.
    exprs = exprs.rstrip('/')
    
    reg_genes = [] # List of regulated genes
    genes     = [] # List of all genes
    
    # This loop checks the expression of all of the genes in each of the 
    # enrichment conditions. It will place regulated genes in reg_genes. All
    # genes are placed in genes, regardless of expression.
    for i in conditions:
    
        # open the raw expression file.
        file = open(exprs+'/'+i, 'r')
        
        # Check the each gene in the file for differential expression.
        for line in file:
            if line.startswith('AT'):
                tab = line.strip().split('\t')
                
                name = tab[0] # gene name
                
                # Get the log2(FC)
                if v == False:
                    M    = float(tab[1]) 
                elif v == True:
                    M    = abs(float(tab[1]))
                
                p = float(tab[2]) # adjusted p-value
                
                # Add the gene to the genes list, regardless of its expression.
                genes.append(name)
                
                # If the fc and p value are above and bellow the threshold, 
                # respectively, call differential expression.
                if M >= fc and p <= pt:
                    reg_genes.append(name)
    
    # Convert the lists to sets
    reg_genes = set(reg_genes)
    genes     = set(genes)
    
    return(genes, reg_genes)

# This function loads the genes in a cluster into a set.
def load_cluster(cluster, complete):
    
    clus_List = [] # List of genes in the cluster
    
    file = open(complete + '/' + cluster, 'r')
    
    # This loop will get the name of all genes in the cluster file.
    for line in file:
        if line.startswith('AT'):
            clus_List.append(line.strip().split('\t')[0])
    
    
    # Convert the list into a set.
    clus_genes = set(clus_List)
    return(clus_genes)
#
################################################################################

################################################################################
# Get the parameters for the function.
import sys, os

print '''~ Parameters
-m   := Directory of raw expression values
-c   := Directory of completed clusters
-f   := File with names of focus expression files
--fc := The log2 fold change for differential expression (default = 1)
-p   := p-value for differential expression (default = 0.05)
-v   := Take the absolute value of log2(fc) as differentially expressed. 
-o   := Output name (default = fisher_table)
'''

fc  = 1.0  # log2(fc) for calling expression
pt  = 0.05 # adjusted p-value threshold for calling expression
v   = False
out = "fisher_table"

for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-m':
        exprs       = sys.argv[i+1] # Expression matrix
    if sys.argv[i] == '-c':
        complete    = sys.argv[i+1] # Directory of completed clusters.
    if sys.argv[i] == '-f':
        cond        = sys.argv[i+1] # File with conditions that you are looking
                                    # for enrichment in.
    if sys.argv[i] == '--fc':
        fc          = float(sys.argv[i+1])
    if sys.argv[i] == '-p':
        pt          = float(sys.argv[i+1])
    if sys.argv[i] == '-v':
        v           = True
    if sys.argv[i] == '-o':
        out         = sys.argv[i+1]
# 
################################################################################

################################################################################
# Run the script (meat and potatoes of the script).

# Load the conditions
conditions = load_conditions(cond)

# Get the sets of all genes and regulated genes
genes, reg_genes = get_reg_sets(exprs, conditions)

# Get the list of all files in the completed clusters directory.
complete_clusters = os.listdir(complete)

output = open(out, 'w')

for cluster in complete_clusters:
    
    # Load the genes in the cluster
    clus_genes = load_cluster(cluster, complete)
    
    # Calculate a, b, c, and d for the Fisher exact test table
    #                | Responsive | Not Responsive |cd 
    #                |------------|----------------|
    #     in cluster |      a     |        b       |
    # ---------------|------------|----------------|
    # out of cluster |      c     |        d       |
    #                -------------------------------
    
    a = str(len(reg_genes  & clus_genes)) 
    b = str(len(clus_genes - reg_genes)) 
    c = str(len(reg_genes  - clus_genes))
    d = str(len(genes - (reg_genes | clus_genes)))
    
    outline = '\t'.join([cluster, a, b, c, d])+'\n'
    
    output.write(outline)
