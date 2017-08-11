# This script will create command lines for the submission of the mapping 
# jobs.

import sys, os

def get_args():
    
    global TAMO_dir, FASTA_file, TAMO_list, map_p, map_thresh, ATfreq, GCfreq
    
    ##
    # Get the required arguments
    TAMO_dir = sys.argv[1].rstrip('/') # Directory with tamo files
    FASTA_file = sys.argv[2] # Directory with FASTA file for mapping.
    TAMO_list = os.listdir(TAMO_dir) # List of all the TAMO files.
    #
    ##
    
    ##
    # Set defaults
    map_p      = '5'
    map_thresh = '0.9'
    ATfreq     = '0.33'
    GCfreq     = '0.17'
    #
    ##
    
    ##
    # Get any changes to the defaults.
    for i in range(3, len(sys.argv)):
        
        if sys.argv[i] == '-p':
            map_p       = sys.argv[i+1]
        if sys.argv[i] == '-t':
            map_thresh  = sys.argv[i+1]
        if sys.argv[i] == '-at':
            ATfreq      = sys.argv[i+1]
        if sys.argv[i] == '-gc':
            GCfreq      = sys.argv[i+1]
    #
    ##
    
def main():

    get_args()
    output = open('run_mapping_cc', 'w')
    for TAMO_file in TAMO_list:

        output.write("module load TAMO; module load motility;"+\
                     " python /mnt/home/seddonal/scripts/6_motif_mapping/"+\
                     "_tamo_mapping_scorePvalue.py %s/%s %s %s %s %s %s\n" % (TAMO_dir, TAMO_file, map_p, map_thresh, ATfreq, GCfreq, FASTA_file))
                     
    output.close()

if __name__ == '__main__':
    main()
