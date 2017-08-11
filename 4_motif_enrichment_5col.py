##
# This script creates the 5col for the Fisher exact test script.
##

import sys, os

##
# checks if the gene in a line is a gene of interest or not.
def check_gene(line, goi_dict, a_dict, b_dict):
	
	gene = line.strip().split('\t')[0] # Get the name of the gene from the line.
	
	
	if gene in goi_dict:
		a_dict[gene] = None
	else:
		b_dict[gene] = None
		
	return a_dict, b_dict
#
##

##
# This script will get the number of genes in the FASTA file used for mapping.
# It can be set to include or exclude plastid and mitochondrial genes
def get_gene_number_from_fasta(FASTA_file_name, only_nuclear = 1):
	
	FASTA_file = open(FASTA_file_name, 'r')
	genes_dict = {}
	
	for line in FASTA_file:
		
		# If the line starts with >AT, then it is a gene name.
		if line.startswith('>AT'):
			
			if only_nuclear == 1:
				
				if line.startswith('>ATM') or line.startswith('>ATC'):
					pass
				else:
					genes_dict[line.strip()] = None
			
			elif only_nuclear == 0:
				genes_dict[line.strip()] = None
			
			else:	
				raise TypeError('Expected 1 or 0 for only_nuclear in '+\
								'get_gene_number_from_fasta.')			
	
	return len(genes_dict)
#
##


##
# Load the genes of interest into a dictionary
def load_goi(goi_file_name):
	goi_dict = {}
	goi_file = open(goi_file_name, 'r')
	for line in goi_file:
		goi_dict[line.strip()] = 1
	return(goi_dict)
#
##

##
# Create and write the 5col file for the Fisher exact tables
# The table is defined as follows
#                     | Regulated gene (goi) | Non regulated gene (non goi)
# Has motif           |          a           |              b
# --------------------|----------------------|------------------------------
# Does not have motif |          c           |              d                   
                       
def create_table(map_dir, goi_dict, goi_total, non_goi_total, out_name):
	
	map_file_list = os.listdir(map_dir)
	
	output = open(out_name, 'w')
	
	# This loop counts the genes of interest in a mapping file, generates
	# the numbers for the contingency table, and then writes them the output.
	
	n = 0 # keeps track of the number of files looked at so far.
	for map_file_name in map_file_list:
		
		a_dict = {} # dict of goi in mapping file.
		b_dict = {} # dict of non-goi in mapping file.
		
		map_file = open(map_dir+'/'+map_file_name, 'r') # open the mapping file

		try:
			next(map_file) # Skip the first line (header).
		# Skip a file if it is empty 
		except StopIteration:
			pass 
		else:
			# This loop counts the number any type of genes with a motif
			for line in map_file:
				
				a_dict, b_dict = check_gene(line, goi_dict, a_dict, b_dict)
					
			a = len(a_dict)
			b = len(b_dict)
		
			# Calculate c and d based on the values of 
			c = str(goi_total - a)
			d = str(non_goi_total - b)
		
			# Write the line in to the 5col file.
			out_line = '\t'.join([map_file_name, str(a), str(b), c, d])+'\n'
			#print out_line
			output.write(out_line)
		
		n += 1
		if n % 1000 == 0:
			print n, "mapping files so far."
		
	output.close()
#
##

##
#
def help():
	print """Parameters:
-d := Mapping directory
-g := File of regulated genes
-s := FASTA file used for mapping
--nuclear := (1) only use nuclear genes, (0) use all genes.
-o := Name of the output file
	"""
#
##

##
# Parse argument
def parse_args():
	
	global map_dir, goi_file_name, FASTA_file_name, out_name, only_nuclear
	
	##
	# Set Defaults
	only_nuclear = 1
	#
	##
	
	for i in range(1, len(sys.argv)):
		
		if sys.argv[i] == '-d':
			map_dir     = sys.argv[i+1]
		elif sys.argv[i] == '-g':
			goi_file_name = sys.argv[i+1]
		elif sys.argv[i] == '-s':
			FASTA_file_name = sys.argv[i+1]
		elif sys.argv[i] == '--nuclear':
			only_nuclear  = int(sys.argv[i+1])
		elif sys.argv[i] == '-o':
			out_name      = sys.argv[i+1]
		elif sys.argv[i] == '-h' or sys.argv[i] == '--help':
			help()
			sys.exit()
#
##

##
# Main function
def main():
	
	parse_args()
	
	goi_dict = load_goi(goi_file_name) # load all of the regulated genes
	
	
	total_gene_num = get_gene_number_from_fasta(FASTA_file_name, only_nuclear)
	print "There are", total_gene_num, "genes loaded from the FASTA."
	
	goi_total     = len(goi_dict) # Total number of regulated genes
	print "There are", goi_total, "regulated genes loaded."
	
	non_goi_total = total_gene_num - goi_total # Total of non regulated
	print "There are", non_goi_total, "non-regulated genes."
	
	create_table(map_dir, goi_dict, goi_total, non_goi_total, out_name)
#
##

if __name__ == "__main__":
	main()
