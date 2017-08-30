##
# This script creates the 5col for the Fisher exact test script.
##

import sys, os
from time import time

##
# checks if the gene in a line is a gene of interest or not.
def check_gene(line, goi_dict, a_dict, b_dict):
	gene = line.strip().split('\t')[0] # Get the name of the gene from the line.
	if gene in goi_dict:
		a_dict[gene] = None
	else:
		b_dict[gene] = None
		
	return a_dict, b_dict

##
# This script will get the number of genes in the FASTA file used for mapping.
# It can be set to include or exclude plastid and mitochondrial genes
def get_fasta_gene_num(FASTA_file_name, only_nuclear, oid):
	
	FASTA_file = open(FASTA_file_name, 'r')
	genes_dict = {}
	
	for line in FASTA_file:
		
		if line.startswith('>'):
			if only_nuclear == 1:
				if oid == "":
					print "You need to provide organellar gene identifier prefixs."
					print "Eg. for Arabidopsis, ATM,ATC"
					sys.exit(0)
				else:
					oid = oid.split(',')
					isorganellar = 0
					for j in oid:
						if line.startswith(j):
							isorganellar = 1
					if isorganellar:
						continue
					else:
						genes_dict[line.strip()] = None
			elif only_nuclear == 0:
				genes_dict[line.strip()] = None
			else:	
				raise TypeError('Expected 1 or 0 for only_nuclear in '+\
								'get_fasta_gene_num.')			
	
	return len(genes_dict)

##
# Load the genes of interest into a dictionary
def load_goi(goi_file_name):
	goi_dict = {}
	goi_file = open(goi_file_name, 'r')
	for line in goi_file:
		goi_dict[line.strip()] = 1
	return(goi_dict)

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
	
	map_file_list.sort()
	
	# Get a list of out.pvalue files
	MFL2 = []
	for i in map_file_list:
		if i.endswith('.out.pvalue'):
			MFL2.append(i)
			
	print "Total %i files" % len(MFL2)	
	fileNum = 0
	for map_file_name in MFL2:
		if fileNum % 100 == 0:
			print ' %i x100' % (fileNum/100)
		fileNum += 1
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
				# Make it faster, build the dict directly.
				#a_dict, b_dict = check_gene(line, goi_dict, a_dict, b_dict)
				
				# Make it faster, don't split.
				#gene = line.split('\t')[0] # Get gene name from the lines.
				gene = line[:line.find('\t')]
				if gene in goi_dict:
					a_dict[gene] = 1
				else:
					b_dict[gene] = 1
				
			a = len(a_dict)
			b = len(b_dict)
		
			# Calculate c and d based on the values of 
			c = str(goi_total - a)
			d = str(non_goi_total - b)
		
			# Write the line in to the 5col file.
			out_line = '\t'.join([map_file_name, str(a), str(b), c, d])+'\n'
			#print out_line
			output.write(out_line)
		
	output.close()

##
#
def help():
	print """Parameters:
		-d  Mapping directory.
		-g  File with a list of DEGs.
		-s  FASTA file used for mapping.
		-n  (1) only use nuclear genes, (0) use all genes.
		-o  Name of the output file.
		-c  If want to parse DEG name based on chr, provide delimiter here.
		-oid organellar gene ID prefixs, separated by ','.
	"""

##
# Parse argument
def parse_args():
	
	global map_dir, goi_file_name, FASTA_file_name, out_name, only_nuclear, oid
	
	##
	# Set Defaults
	only_nuclear = 0
	oid = c_delim = ""
	
	for i in range(1, len(sys.argv)):
		
		if sys.argv[i] == '-d':
			map_dir     = sys.argv[i+1]
		elif sys.argv[i] == '-g':
			goi_file_name = sys.argv[i+1]
		elif sys.argv[i] == '-s':
			FASTA_file_name = sys.argv[i+1]
		elif sys.argv[i] == '-n':
			only_nuclear  = int(sys.argv[i+1])
		elif sys.argv[i] == '-o':
			out_name      = sys.argv[i+1]
		elif sys.argv[i] == '-oid':
			oid      = sys.argv[i+1]
		elif sys.argv[i] == '-h' or sys.argv[i] == '--help':
			help()
			sys.exit()

##
# Main function
def main():
	
	parse_args()
	
	goi_dict = load_goi(goi_file_name) # load all of the regulated genes
	
	print goi_file_name
	gene_num  = get_fasta_gene_num(FASTA_file_name,only_nuclear,oid)
	goi_total = len(goi_dict) # Total number of regulated genes
	ngoi_total = gene_num - goi_total # Total of non regulated
	print " %i genes total" % gene_num 
	print " %i DEG"         % goi_total
	print " %i non-DEG"     % ngoi_total
	
	create_table(map_dir, goi_dict, goi_total, ngoi_total, out_name)

if __name__ == "__main__":
	main()
