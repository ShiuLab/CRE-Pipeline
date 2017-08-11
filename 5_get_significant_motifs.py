import sys
from TAMO.MotifTools import txt2motifs
from TAMO.MotifTools import Motif, save_motifs

##
# This function should parse out features that are significant. It requires a
# table with features and associated p-values. The user specifies an alpha,
# and the script gives an output of the fignificant features.
##

print '''
-in      := input file
-a       := alpha level
-index   := The location of the adjusted p-vaqlue or q-value
-output  := The name of the output file
'''

##
#
def parse_args():

	global FET_file_name, alpha, index, output_name
	
	for i in range(0, len(sys.argv)):
		if sys.argv[i] == '-in':
			FET_file_name = sys.argv[i+1]
		if sys.argv[i] == '-a':
			alpha       = float(sys.argv[i+1])
		if sys.argv[i] == '-index':
			index       = int(sys.argv[i+1])
		if sys.argv[i] == '-output':
			output_name = sys.argv[i+1]

##
# This function extracts the name of a motif from a mapping file.
def get_motif_name(map_file_path):
	
	map_file = open(map_file_path, 'r') # open the map file.
	
	# The first line in the mapping file contains the motif name.
	line = next(map_file)
	tab = line.strip().split('\t')
	
	motif = tab[1]
	
	return(motif)
#
##

##
# This function looks at the oupout of a fisher exact test table, identifies which
# mapping file is for a motif that is significantly enriched. The function than
# pulls the motif name from the mapping file and adds it to the significant
# motif dictionary.
def get_significant_motifs(FET_file_name, alpha, map_dir):
	FET_file = open(FET_file_name, 'r')
	
	sig_motif_dict = {} # Dictionary containing the names of the significant 
						# motifs.
	
	for line in FET_file:
		
		# Get tab delimited line.
		tab = line.strip().split('\t')
		map_file_name = tab[0]
		
		pvalue = float(tab[7]) # adj p-value for motif
		
		direction = tab[5] # Is the enrichment greater (+) or less than (-)
		
		# If the motif is significant, add it to the dictionary.
		if pvalue <= alpha and direction == '+':
			
			map_file_path = map_dir + '/' + map_file_name
			motif_name    = get_motif_name(map_file_path)
			
			# Add the motif to the dictionary
			sig_motif_dict[motif_name] = None
	
	return(sig_motif_dict)
#
##

##
# create a motif index dictionary for a motif list
def create_motif_index(motif_list):

	motif_index_dict = {}
	
	for i in range(0, len(motif_list)):
	
		motif_name = str(motif_list[i]).split(' ')[0]
		motif_index_dict[motif_name] = i
		
	return motif_index_dict
#
## 

##
# 
def create_significant_motif_TAMO(input_TAMO_file_name, sig_motif_dict):
	
	# Load the moitfs
	motif_list = txt2motifs(input_TAMO_file_name)
	
	# Create motif index dictionary
	motif_index_dict = create_motif_index(motif_list)
	
	
	# Get the list of significant motifs from the motif_list
	sig_motif_list = []
	for motif in sig_motif_dict:
		motif_index = int(motif_index_dict[motif])
		sig_motif_list.append(motif_list[motif_index])
	
	save_motifs(sig_motif_list, input_TAMO_file_name+'.sig_enriched')
#
##

##
# Parse arguments.
def parse_args():

	global FET_file_name, alpha, map_dir, input_TAMO_file_name
	
	##
	# Set defaults
	alpha = 0.05
	#
	##

	for i in range(0, len(sys.argv)):
		if sys.argv[i] == '-in':
			FET_file_name = sys.argv[i+1]
		elif sys.argv[i] == '-a':
			alpha       = float(sys.argv[i+1])
		elif sys.argv[i] == '-index':
			index       = int(sys.argv[i+1])
		elif sys.argv[i] == '-tamo':
			input_TAMO_file_name = sys.argv[i+1]
		elif sys.argv[i] == '-map_dir':
			map_dir  = sys.argv[i+1]	

def main():
	
	parse_args()
	
	sig_motif_dict = get_significant_motifs(FET_file_name, alpha, map_dir)
	print "There are", len(sig_motif_dict), "significantly enriched motifs."
	create_significant_motif_TAMO(input_TAMO_file_name, sig_motif_dict)

if __name__ == "__main__":
	main()

