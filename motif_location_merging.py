'''Functions used to identify the percentage of times that motifs overlap.

The script contains a set of functions that allow you to calculate the overlap
of all motifs from a set of mapping files. The functions are based on scripts 
written by Cheng Zou. They have been inserted into this script with some 
significant changes. Namely, the overlap scripts now accepts motifs mapped by
motility, and not MAST.

Usage:
python motif_location_merging.py [Function]

Functions:
	sim_matrix: Use this function to identify the number of overlaps between 
	each motif. Requires -i, -m, -g, -t
	
	run_as_job: Runs sim_matrix as a job. The walltime for the job can be 
	estimated by running sim_matrix to determine the average time necessary
	to run 5000 combinations. Your total combinations will equal n choose 2, 
	where n is the number of motifs. The walltime could thus me estimated as
	 ~~(n choose 2) * (average time to run 5000 combinations) / 5000 = walltime.
	Requires: -i, -m, -g, -t, -wall
	
	create_matrix: Takes the SIM file and turns it into a matrix. Requires: -s
	
	remove_duplicates: Removes duplicate motifs from a TAMO file. Should be 
	performed on any TAMO file before proceeding with the motif location 
	merging. Requires: -t
	
	mapped_TAMO: Creates a TAMO containing only the motifs that have been mapped
	to the genome. The output will end with '.mapped.tm'. Should be performed
	on all TAMO files before proceeding with motif location merging. This 
	function will also do the same processes as remove_duplicates, so you only
	need to use this function. Requires: -t, -i
	
	run_UPGMA: Run clustering on the matrix created with create_matrix. The 
	function will create a job that will be submitted with the R script used for
	clustering motifs. Requires: -mat, -h. Options: -w (working directory is
	current directory by default), -wall (wall time default is 120), -mem 
	(memory default is 124 gb).
	
Parameters:
	-g := Gene list. The list of genes that you are intersted in using to 
		identify overlap.
	-i := File that gives the mapping file for each motif.
	-m := Directory where the mapping files are located.
	-mat := Matrix created from the SIM file.
	-mem := memory requirement for running jobs.
	-s := SIM file.
	-t := TAMO file containing all the motifs.
	-w := Working directory
	-wall := walltime
	
Example Protocol:
This is a typical set of commands for merging motifs based on location.

##
# Get a TAMO of only the mapped motifs. If this step is not performed, then
# your matrix will not match the indicies in your TAMO file.

python motif_location_merging.py mapped_TAMO -t TAMO_file.tm -i motif_mapping_index


##
# Run sim_matrix to find the estimate of the wall time for running sim_matrix
# as a job. Use crtl Z to kill the process after it has printed out an estimated 
# walltime. If you let the process run to long, HPCC will kill it before it is
# done.

python sim_matrix -i motif_mapping_index -m motif_mapping_directory -g gene_list
-t TAMO_file.mapped.tm
ctrl Z

##
# With your walltime estimate, you can now run the above script as a job using 
# the following command.

python run_as_job -i motif_mapping_index -m motif_mapping_directory -g gene_list
-t TAMO_file.mapped.tm -wall 300
'''

import sys, time, os
from TAMO import MotifTools

class motif_mapping_methods(object):

	def __init__(self, mapping_index_file):
		self.mapping_index_file = mapping_index_file
		
	def get_dict(self):
		mdic  = {} # Dic of motifs = {motif name: mapping file}

		for motif_line in open(self.mapping_index_file, 'r'):
	
			# For each line that is not empty
			if motif_line.strip('\n') != '':
	
				# split the line. Index 0 is the motif name, and 1 is the map
				# file name.
				motif_line_list = motif_line.strip().split('\t')

				mdic[motif_line_list[1]] = motif_line_list[0]
		
		return (mdic)

class motif_mapping(motif_mapping_methods):

	'''Creates a dictionary that links motifs to their mapping file.
	
	Also creates a list of all motifs that are found in a mapping file.
	'''
	def __init__(self, mapping_index_file):
		self.mapping_index_file = mapping_index_file
		self.dict = self.get_dict()
		self.list = self.dict.keys()

def xuniqueCombinations(items, n):

	'''Creates list of tuples of all unique combinations of n elements in items.
	
	This function is used to pair up all of the motifs into all possible 
	combinations.
	
	[(start, stop), ...]
	'''

	if n==0:
		yield []
	else:
		for i in xrange(len(items)):
			for cc in xuniqueCombinations(items[i+1:],n-1):
				yield [items[i]]+cc

def tuple_compare(x, y):

	'''Compares the coordinates of two tuples, x and y.
	
	If x starts after y, 1 is returned
	If x and y start at the same spot, 1 is returned if x ends after y, and
	0 is returned if x and y are completely the same
	If x starts before y, -1 is returned.
	
	The script allows a list of tuples to be sorted so that the earlier tuples
	are placed toward the front of the list.
	
	example: (200,209) (200,208)-> (200,208) (200,209)
	'''
	
	# If x start is bigger than y start
	if x[0]>y[0]:
		return 1
		
	# If the starts are in the same spot
	elif x[0]==y[0]:
		
		# If x end is bigger then y end
		if x[1]>y[1]:
			return 1
		
		# If the ends are in the same spot
		elif  x[1]==y[1]:
			return 0
	
		else:
			return -1
	else: # x<y
		return -1

def remove_duplicate_motifs(TAMO_file):
	'''Removes motifs that are duplicates in a TAMO.
	
	Returns a TAMO file that contains only one of each motif in another TAMO
	file.
	'''
	
	# Acquire motifs from TAMO_file
	ml = MotifTools.txt2motifs(TAMO_file)
	
	# Look at each motifs. It has been seen before, it is excluded from the new
	# TAMO file.
	tmp_ml = []
	motif_list = []
	for m in range(len(ml)):
		motif = ml[m].oneletter
		if not motif in motif_list:
			motif_list.append(motif)
			tmp_ml.append(ml[m])
	ml = tmp_ml
	MotifTools.save_motifs(ml,TAMO_file+'.no_dups.tm')
	
def TAMO_of_mapped_motifs(TAMO_file, motif_index_file):
	
	'''Creates a TAMO file of motifs that are in a mapping index
	
	Creates a TAMO containing only the motifs that have been mapped
	to the genome. The output will end with '.mapped.tm'. Should be performed
	on all TAMO files before proceeding with motif location merging. This 
	function will also do the same processes as remove_duplicates, so you only
	need to use this function.
	'''
	
	# Acquire motifs from TAMO_file
	ml = MotifTools.txt2motifs(TAMO_file)
	
	# Get a list of all motifs in the mapping directory
	mapped_motif_list = motif_mapping(motif_index_file).list
	
	
	# Check remove motifs in ml that are not in the mapping directory.
	tmp_ml = []
	for m in range(len(ml)):
		
		# Get IUPAC representation of motif
		motif = ''
		for letter in ml[m].oneletter:
			if letter == '.': motif = motif + 'N'
			else: motif = motif + letter
		
		if motif in mapped_motif_list:
			tmp_ml.append(ml[m])
	
	MotifTools.save_motifs(tmp_ml, TAMO_file + '.mapped_motif.tm')
	
	# create a new mapping file
	

def create_motif_order_dict(TAMO_file):

	'''Creates a dictionary that links motifs names to their number in a TAMO.
	
	it returns the following dictionary: mdic = {motif name: motif number}
	'''
	
	ml = MotifTools.txt2motifs(TAMO_file)
	mdic = {}
	for m in range(len(ml)):
		motif = ml[m].oneletter
		tmp_motif = ''
		for letter in motif:
			if letter == '.':
				tmp_motif = tmp_motif + 'N'
			else:
				tmp_motif = tmp_motif + letter
		motif = tmp_motif
		mdic[motif] = str(m)
	return mdic
	

def extract_motif_location_information(map_file_path, gl_dic):
		
	'''Extracts information about location of each motif.
	
	For each gene that the motif maps to, the script creates a list of tuples
	of the coordinates on that genes promoter. This function is used as part
	of motif_gene_mast().
	'''
	# Look at each line in the mapping filem
	
	# Open the file and skip the first line
	map_file = open(map_file_path, 'r')
	next(map_file)
	
	# Get the map information from each subsequent line.
	for line in map_file:
	
		# Get the required information from the line
		line_list = line.strip('\n').split('\t')
		gene_name = line_list[0]
		start	  = line_list[1]
		stop	  = line_list[2]
		p_value   = line_list[6]
		
		if abs(int(start)-int(stop))+1 == 6:
			threshold = 0.005
		else:
			threshold = 1e-5
		
		if float(p_value) < threshold and gl_dic.has_key(gene_name):
			
			# Put the gene in the temperary dictionary if it is not their
			# already.
			if not temp_dic.has_key(gene_name):
				temp_dic[gene_name]={}
		
			# Add the left border and the p_value to the gene in temp_dic.
			temp_dic[gene_name][int(start)]=float(p_value)

			# If the motif_map_dic dies not have the gene, create a key for it.
			if not motif_map_dic.has_key(gene_name):
				motif_map_dic[gene_name]= []
				
			# Add the location information to the gene
			motif_map_dic[gene_name].append((int(start),int(stop)))


def Motif_gene_mast(map_file_name, mapping_dir_root, gl_dic): 
	
	'''This function links up motifs that are on the same gene.
	
	It takes the motif dictionary as input, identifies all ot the genes that
	each motif maps to, along with its coordinates, and outputs this information
	to a dictionary.
	
	The product of this script is a dictionary with each gene the motif is 
	mapped, and a list of the coordinates on the gene.
	
	motif_map_dic = {gene:[(start, end), ...], ...}
	'''
	
	global motif_map_dic, temp_dic

	motif_map_dic = {} 
	temp_dic  = {} # gene[200][pvaue]
	mapping_dir_root = '/'+mapping_dir_root.strip('/')
	
	#print map_file_name
	
	# Look at each motif in motif_map_dict
	map_file_path = mapping_dir_root + '/' + map_file_name
	extract_motif_location_information(map_file_path, gl_dic)
		
	# For each motif in the motif_dic
	# for motif in motif_dic.keys():
	
	# For each gene that the motif maps to.
	for gene in motif_map_dic.keys():
	
		list= motif_map_dic[gene] # List of all the coordinate tuples for the gene
		list.sort(tuple_compare)
		
		# If there is more than one coordinate tuple in the list.
		if len(list)>1:
		
			index_list = range(len(list))
			end_of_list = False
			while end_of_list == False:
				k = index_list[0]
				try: list[k+1]
				except IndexError:
					end_of_list == True
					break
				else:
					if abs(list[k][0]-list[k][1])>16:
						list.pop(k+1)
						index_list = range(len(list))
					else:
						index_list.pop(0)

			
			##
			# Remove duplicate motifs
			end_of_list = False
			index_list = range(len(list))
			while end_of_list == False:
	
				k = index_list[0]
	
				try: list[k+1]
				except IndexError:
					end_of_list == True
					break
				else:
					if list[k] == list[k+1]:
						list.pop(k+1)
						index_list = range(len(list))
					else:
						index_list.pop(0)
			
			
			##
			# Choose the most significant motif of motifs that overlap.
			end_of_list = False
			index_list = range(len(list))
			while end_of_list == False:
				k = index_list[0]
	
				try:
					list[k+1]
				except IndexError:
					end_of_list == True
					break
				else:
					if list[k][1] >= list[k+1][0]:
			
						if temp_dic[gene][list[k][0]] <= temp_dic[gene][list[k+1][0]]:
							list.pop(k+1)
						else:
							list.pop(k)
			
						index_list = range(len(list))
					else:
						index_list.pop(0)
			
		motif_map_dic[gene]=list
	return motif_map_dic

def overlapping(list1,list2):
	
	'''Identifies coordinate tuples that overlap inbetween two lists of tuples.
	'''
	
	# list1 and list2 are the lists of all the cordinates for motif1 and motif2.
	# Each coordinate is a tuple of the start and stop position.
	tlist = []
	n=0
	
	for i in list1:
		while list2[n][0]<i[0] and list2[n][1]<i[0]:
			n+=1
			if n==len(list2):
				break
		if n==len(list2):
				break
		#print i,list2[n]
		if list2[n][0]>=i[0] and list2[n][0]<= i[1]:
			#print "y"
			tlist.append(1)
		elif list2[n][0]<=i[0] and list2[n][1]>= i[0]:
			#print "y"
			tlist.append(1)
		else:
			continue

	return tlist.count(1)

def len_value(dic):
	
	'''Returns the sum of lenght of all value for all keys in a dictionary.
	
	In the case of the motif_dict, this would be the total number of instances
	of a motif in the genome.
	'''
	
	n=0
	for i in dic.keys():
		n+=len(dic[i])
	return n
	
def sim_matrix(motif_index_file, mapping_dir_root, gene_list_file, TAMO_file):	
	'''Determines the number of overlaps between each unique pair of motifs.
	
	Creates every possible combination of motifs in a list of motifs. It then
	goes through each combination and identifies everytime the motifs overlap.
	
	The output is a SIM file with three columns. The first colum is a comma
	seperated combination of the motifs, the next is the total number of
	overlaps, and the last is the number of times that each motifs appears in 
	the genome.
	'''
	# Load the TAMO file
	TAMO_index = create_motif_order_dict(TAMO_file)
	
	##
	# Convert the gene_list_file to a dictionary
	gl_dic = {}
	for line in open(gene_list_file, 'r'): 
		gl_dic[line.strip()] = "1"

	print "output: %s_%s.SIM" % (motif_index_file,gene_list_file.split("/")[-1])
	
	##
	# build mlist and mdic
	motif_map_index = motif_mapping(motif_index_file)
	mlist = motif_map_index.list
	mdic  = motif_map_index.dict
	
	##
	# Check that the TAMO file contains only the motifs that are in the mapping
	# files.
	create_new_tamo = False
	for motif in TAMO_index.keys():
		if not motif in mlist:
			print 'The TAMO file you submited contains contains motifs that \
are not present in the mapping files.'
			create_new_tamo = True

	if create_new_tamo == True:
		TAMO_of_mapped_motifs(TAMO_file, motif_index_file)
		TAMO_index = create_motif_order_dict(TAMO_file + '.mapped_motif.tm')
		print 'A new TAMO file named', TAMO_file + '.mapped_motif.tm.'
		print 'the new TAMO will now be used. It should also be used in \
subsequent analysis.'

	for motif_line in open(motif_index_file, 'r'):
	
		# For each line that is not empty
		if motif_line.strip('\n') != '':
	
			# split the line. Index 0 is the motif name, and 1 is the map
			# file name.
			motif_line_list = motif_line.strip().split('\t')
		
			mlist.append(motif_line_list[1])
			mdic[motif_line_list[1]] = motif_line_list[0]
	#
	##
	
	##
	# Create motif_dic, which contains the genes and locations of all of the
	# motifs in mlist.
	# motif_dic = {motif name:{gene:[(start, end), ...], ...}, ...}
	motif_dic = {}
	n = 0
	for motif in mlist:
		n += 1
		map_file_name = mdic[motif]
		motif_dic[motif] = Motif_gene_mast(map_file_name, mapping_dir_root, gl_dic)
	print 'Motif location information extracted'
	#
	##
	
	##
	# Create dic1 and dic2. They contain all the combinations of motifs 
	# as a set of tuples. dic1 will contain a count of how many times each
	# combination of motifs overlaps. dic2 contains a count of all the
	# instances of each motif in a pair.
	dic1 = {} # {(motfi1, motif2): # of overlaps}

	
	for i in list(xuniqueCombinations(mlist, 2)):
		dic1[tuple(i)] = 0
	print 'All motif combinations are in tuples.'
	#
	##
	
	##
	# Calculate the values for dic1. Look at each tuple of motif combinations
	
	combinations = len(dic1.keys())
	print combinations
	n = 0
	total_time = time.time()
	oup=open("%s_%s.SIM" % (motif_index_file,gene_list_file.split("/")[-1]),"w")
	for ituple in dic1.keys():
		
		##
		# Count the number of tuples seen so far. Print a current count and
		# time for every 5000 tuples.
		n += 1
		if n%5000 == 0:
			print n, 'combinations checked so far.'
			
			estimated_wall_time = combinations*(time.time()-total_time)/n/60
			
			print 'Estimated wall time:', estimated_wall_time, 'min'
		
		##
		
		# Pulls out motifs in the tuple
		motif_1 = ituple[0]
		motif_2 = ituple[1]
		overlaps = 0
	

		# If motifs are not in motif_dic, the tuple has 0 overlaps.
		if not motif_dic.has_key(motif_1) or not motif_dic.has_key(motif_2):
			overlaps=0
		
		
		# Otherwise, count the overlaps.
		else:
		
			# Filter to find genes that overlap in the motifs.
			for og in filter(motif_dic[motif_1].has_key,motif_dic[motif_2].keys()): # og: overlapped genes
		
				# Overlapping returns the number of times that the motifs overlap
				# on the same coordinates. This is placed in the index of the
				# tuple in the dictionary.
				overlaps += overlapping(motif_dic[motif_1][og], motif_dic[motif_2][og])
				
		
		##
		# Count the total number of times that each motif appears in the genome
		# total_ins = apperances of motif1 + app of motif 2 - overlaps
		if not motif_dic.has_key(motif_1) and motif_dic.has_key(motif_2):
			total_ins = len_value(motif_dic[motif_2])
		elif not motif_dic.has_key(motif_2) and motif_dic.has_key(motif_1):
			total_ins = len_value(motif_dic[motif_1])
		elif motif_dic.has_key(motif_2) and motif_dic.has_key(motif_1):
			total_ins = len_value(motif_dic[motif_1]) + len_value(motif_dic[motif_2])- overlaps
		else:
			total_ins = 0

		motif_1_index = TAMO_index[motif_1]
		motif_2_index = TAMO_index[motif_2]		
		oup.write("%s,%s\t%s\t%s\n" % (motif_1_index, motif_2_index, overlaps, total_ins))
	
	#oup=open("%s_%s.SIM" % (motif_index_file,gene_list_file.split("/")[-1]),"w")
	#for i in dic1.keys():
	#	oup.write("%s,%s\t%s\t%s\n" % (i[0],i[1],dic1[i],dic2[i]))
	oup.close()
	print time.time() - total_time

def run_as_job(motif_index_file, mapping_dir_root, gene_list_file, TAMO_file, walltime, memory):
	'''This script allows you to run the SIM matrix script as a job.
	
	motif_location_merging.py run_as_job [motif map index] [mapping root 
	directory] [gene list] [wall time]
	
	You should estimate the walltime first by running the script normally and
	finding the average time it takes to run 5000 combinations, and then 
	multiplying this by the total number of combinations.
	
	The total number of combinations will be at most n choose 2, where n is the
	number of motifs that you are checking.
	'''
	motif_index_file = os.path.abspath(motif_index_file)
	mapping_dir_root = os.path.abspath(mapping_dir_root)
	gene_list_file   = os.path.abspath(gene_list_file)
	TAMO_file        = os.path.abspath(TAMO_file)
	
	##	
	# Load the motif index file
	#motif_index_file = sys.argv[2]
	print 'motif index:', motif_index_file

	script_dir = os.path.abspath(__file__) # path to pcc_merge_CC.py script
	script_dir = '/'.join(script_dir.split('/')[:-1]) # path to script directory
	
	##
	# Load mapping directory root, which is the root to the path for all of
	# the mapping files.
	#mapping_dir_root = sys.argv[3]
	print 'mapping directory root:', mapping_dir_root
	
	
	# Load the gene list file
	#gene_list_file = sys.argv[4]
	print 'gene list:', gene_list_file
	
	#TAMO_file = sys.argv[5]
	
	#walltime = sys.argv[6]
	
	# Load Walltime
	h = str(int(walltime) / 60)
	m = str(int(walltime) % 60)
	
	oup = open('merge_location_%s.sh' % motif_index_file.split('/')[-1], 'w')
	
	oup.write('#!/bin/sh --login\n\n')
	oup.write('#PBS -q main\n')
	oup.write('#PBS -l nodes=1:ppn=1,walltime=%s:%s:00,mem=%sgb\n' % (h, m, mem))
	#job.write("#PBS -l nodes=1:ppn=1,walltime=%s:%s:00,mem=%sgb\n" % (w_hour, w_min, mem))

	oup.write("module load TAMO; python %s/motif_location_merging.py sim_matrix -i %s -m %s -g %s -t %s" % \
	         (script_dir, motif_index_file, mapping_dir_root, gene_list_file,\
	          TAMO_file))
	oup.close()
	
	os.system("qsub merge_location_%s.sh" % (motif_index_file.split('/')[-1]))

def create_matrix(sim_file_path):
	'''Makes a matrix out of the SIM file.
	
	SIM files are four column files, with each line representing one cell in a
	matrix. This function goes through the SIM file, and creates a matrix
	that can then be run through R to do clustering.
	
	The function takes the information in the SIM file to calculate the 
	percentage of non-overlap. For example, if two motifs overlap 65% of the
	time, the percentage of non-overlap would be 35%. This can be used as a
	distance for purposes of clustering.
	'''
	
	mDic = {}
	
	sim_file = open(sim_file_path, 'r')
	
	for line in sim_file:
		
		# split the line
		tab = line.strip().split('\t')
		
		# extract the information from the line
		motif1, motif2 = tab[0].split(',')
		
		motif1, motif2 = (int(motif1), int(motif2))
		# Calculate the percentage of overlap
		try:
			overlap_percentage = str(100-float(tab[1])/float(tab[2])*100)
		except ZeroDivisionError:
			overlap_percentage = '100.0'
		
		# Create keys for motifs if non exist
		if not motif1 in mDic:
			mDic[motif1] = {}
		if not motif2 in mDic:
			mDic[motif2] = {}
		
		# Put the overlap percentages in to the dictionary
		mDic[motif1][motif2] = overlap_percentage
		mDic[motif1][motif2] = overlap_percentage
	
	
	sim_file.close()
	
	# Create the output file
	motif_overlap_matrix_file = open(sim_file_path + '.matrix', 'w')
	
	motif_list = mDic.keys()
	
	motif_list.sort()
	#print motif_list
	print len(motif_list)
	
	for n in range(0, motif_list[-1]+1):
		if not n in motif_list: print 'motif', n, 'is not in motif_list.'
	
	# check that all motfs are included
	print 'It is', len(motif_list) == motif_list[-1]+1, 'that all motifs are \
included in the matrix.'
	
	for motif1 in motif_list:
	
		# Create the output line that will be writen to the output file.
		output_line = ''
		
		# Get the overlap dict for the motif
		motif_overlap_dict = mDic[motif1]
		
		for motif2 in motif_list[:motif_list.index(motif1)]:
			output_line = output_line + '\t-'
		
		output_line = output_line + '\t0.0'
		
		
		# Go through the motif_list again to find overlap scores and add them to
		# the output_line.
		for motif2 in motif_list[motif_list.index(motif1)+1:]:
		
			if motif2 in motif_overlap_dict:
				output_line = output_line + '\t' + motif_overlap_dict[motif2]
			elif not motif2 in motif_overlap_dict:
				output_line = output_line + '\t100.0'
		
		motif_overlap_matrix_file.write(output_line.strip('\t') + '\n')
	
	motif_overlap_matrix_file.close()

def run_UPGMA(matrix, height, wdir, walltime = 120, mem = '124'):
	
	'''This script creates a tree for motifs based on their PCC distance.
	
	The function works on the output of the combine_distance_matrix and 
	check_square_matrix function. The tree that is output is used by the 
	merge_runs function.
	
	The construction of the tree is performed by an R script, and it is submited
	as a job.
	
	This function is based on a script writen by Cheng Zou, and then converted
	to a function by Alex Seddon.
	'''
	
	# Get the path to pcc_merge_CC.py directory
	script_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])
	
	w_hour = str(walltime / 60)
	w_min  = str(walltime % 60)
	
	job = open("%s_%s.UPGMA.sh" % (matrix, height), 'w')
	job.write("#!/bin/sh --login\n\n")
	job.write("#PBS -q main\n")
	job.write("#PBS -l nodes=1:ppn=1,walltime=%s:%s:00,mem=%sgb\n" % (w_hour, w_min, mem))
	
	job.write("R --vanilla --slave --args %s/%s  %s < %s/UPGMA_final.R> %s.Rout" % (wdir,matrix,height,script_dir,matrix))
	
	job.close()
	
	os.system("qsub %s_%s.UPGMA.sh" % (matrix, height))

def trim_motif(TAMO_file, cut = 0.4):

	'''Trims the motifs in TAMO_file, eliminating low-information flanks.'''
	
	testmotifs = MotifTools.load(TAMO_file)
	file=TAMO_file+"_"+str(cut)+".trim"

	new_mlist=[]
	for motif in testmotifs:
		m = motif.trimmed(cut)
		new_mlist.append(m)
	save_motifs(new_mlist,file)


def pick_chunk_score(wdir, TAMO_file, target, genome):

	'''Trims and returns the top motif in a cluster.
	
	This script takes in the TAMO file from the motifs in a single cluster. It
	trims the low-information ends from each motifs. It then indentifies the
	motif that is most significantly represented in the target genes in your
	genome. If no motif is significantly represented, then a blank top motif
	file is created.
	'''
	os.system("cd %s" % wdir)
	os.chdir(wdir)

	script_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1]) # path to pcc_merge_CC.py script
	
	##
	# step 1 trim tamo to eliminate low information flanking sequence
	trim_motif(TAMO_file, 0.1)
	
	##
	# step 2 Group Specificity Score" from the Church lab
	# python MotifMetrics.py [Genes of interest] -genome [FASTA of promoter sequence] -t [Trimmed TAMO of cluster motifs] 
	# MotifMetrics.py checks if the motifs appear disproportionatly to the 
	# targets compared to the rest of the genes.
	os.system("python %s/MotifMetrics.py %s -genome %s -t %s_0.1.trim -spec > %s_0.1.trim_Cout" % (script_dir,target,genome,TAMO_file,TAMO_file))

	##
	# Gets the motif that is most significantly represented in your target genes
	# Returns "None" if none of the motifs has a p-value above 0.001.
	topm=parse_out_pcs("%s_0.1.trim_Cout" % TAMO_file)
	print "topm",topm
	
	##
	# Writes the top motif to its own directory.
	if topm!="None":

		newdic={}
		ml=MotifTools.txt2motifs("%s_0.1.trim" % TAMO_file)

		for m in ml:
			
			if m.oneletter == topm:
				newdic[m.oneletter] = m

		save_motifs(newdic.values(),"%s.TOP" % TAMO_file)
		os.system("rm %s_0.1.trim" % TAMO_file)
		os.system("rm %s_0.1.trim_Cout" % TAMO_file)
	
	##
	# Writes a blank document if there was no top motif.
	else:
		oup=open("%s.TOP" % TAMO_file,"w")
		oup.close()

def merge_runs_cc(TAMO_file, wdir, ancestor, target, genome, clustering_file):
	
	'''This script is used to merge motifs with the PCC matrix of all motifs.
	
	The script was originally written by Cheng Zou, and then converted to a 
	function by Alex Seddon.
	'''
	
	print "Here are the parameters you specified in this run "
	print "-tamo		%s" % TAMO_file
	print "-wdir		%s" % wdir
	print "-h		height to cut the tree, %s" % height
	print "-ancestor	%s" % ancestor
	print "-target	%s" % target
	print "-genome	%s" % genome
		
	if TAMO_file == '' or wdir == '':
		print __doc__
		sys.exit()

	os.system("cd %s" % wdir)

	os.chdir(wdir)
	
	# Get the directory where the script is located.
	script_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])

	cl_dic={}
	n=0
	
	# The file, clustering file, is inorder of the motifs that appear
	# in the TAMO_file. If two motifs have the same number, they are considered 
	# a part of the same cluster.
	# This loop pulls the clustering information out of this file and creats
	# the dictionary cl_dic = {cluster_index:{motif_index:'1'}}
	for line in open("%s" % (clustering_file),"r"):
		
		# Gets the clusterindex of this motif
		cl=line.strip()
		
		# Adds the cluster index if it has not been 
		if not cl_dic.has_key(cl):
			cl_dic[cl]={}
			
		cl_dic[cl][n]="1" # Adds the motif to that cluster
		n+=1 # Increases the motif index for the next motif

	#print cl_dic
	
	ml=MotifTools.txt2motifs(TAMO_file)
	old=[] # List of motifs that are the sole members of a cluster.
	
	# I think I can divide up this portion of the code to create a series 
	print ancestor,ancestor==0
	
	cc_output = open('merge_runs_cc', 'w')
	
	if ancestor==0:
	
		# This loop Looks at each cluster and attempts to merge the motifs
		# in the cluster if there are multiple motifs.
		for i in cl_dic.keys():
			
			print i,cl_dic[i]
			
			# If there are multiple motifs in the cluster, it merges the motifs
			if len(cl_dic[i])>1:
			
				# Adds all of the motifs in the cluster to an object called 
				# mlist.
				mlist=[] 
				for j in cl_dic[i]:
					mlist.append(ml[j])
					
				# Saves these motifs to there own TAMO file.
				save_motifs(mlist,"%s_sub_%s.tm" % (TAMO_file,i))
				
				cc_output.write('module load TAMO; python %s/motif_location_merging.py merge_runs_no_ancestor -t %s/%s -i %s -target %s -genome %s\n' % (script_dir, wdir, TAMO_file, i, target, genome))
				
			# If there is only one motif in the cluster, it leaves it alone, 
			# And adds it to old
			else:
				key=cl_dic[i].keys()[0]
				old.append(ml[key])

	if ancestor==1:
	
		# This loop Looks at each cluster and attempts to merge the motifs
		# in the cluster if there are multiple motifs.	
		for i in cl_dic.keys():
		
			print i,cl_dic[i]

			# If there are multiple motifs in the cluster, it merges the motifs			
			if len(cl_dic[i])>1:

				# Adds all of the motifs in the cluster to an object called 
				# mlist.
				mlist=[]
				for j in cl_dic[i]:
					mlist.append(ml[j])
				
				# Saves these motifs to there own TAMO file.				
				save_motifs(mlist,"%s_sub_%s.tm" % (TAMO_file,i))
				
				cc_output.write('module load TAMO; module load STAMPmotif; python %s/motif_location_merging.py merge_runs_ancestor -t %s/%s -i %s -target %s -genome %s\n' % (script_dir, wdir, TAMO_file, i, target, genome))
				
			else:
				key=cl_dic[i].keys()[0]
				old.append(ml[key])

	# Combine together the motifs that are in there own cluster.
	#os.system("cat %s_sub_*_sum.tm.tf.tm.TOP > %s_sub_new.tm" % (TAMO_file,TAMO_file))
	save_motifs(old,"%s_sub_old.tm" % (TAMO_file))
	#os.system("cat %s_sub_old.tm %s_sub_new.tm > %s_P1.tm" % (TAMO_file,TAMO_file,TAMO_file))

def merge_runs_no_ancestor(wdir, TAMO_file, i, target, genome):
	
	# I am fairly certain that this process of converting to TF and
	# then returning it to TAMO format is only for keeping the names 
	# consistent. I need to verify this suspicion, AS.
	tamo2tf(wdir, "%s/%s_sub_%s.tm" % (wdir, TAMO_file, i))
	os.system("cat  %s/%s_sub_%s.tm.tf > %s/%s_sub_%s_sum.tm.tf" % (wdir,TAMO_file,i,wdir,TAMO_file,i))
	tf2tamo(wdir, "%s_sub_%s_sum.tm.tf" % (TAMO_file,i))

	# Gets the top motif in the cluster.
	pick_chunk_score(wdir, '%s/%s_sub_%s_sum.tm.tf.tm' % (wdir,TAMO_file,i), target, genome)
	
	# Removes the files that were created.
	os.system("rm %s_sub_%s_sum.tm.tf.tm" % (wdir,TAMO_file,i))
	os.system("rm %s/%s_sub_%s_sum.tm.tf" % (wdir,TAMO_file,i))
	os.system("rm -R %s/%s_sub_%s.tm.tf_ST*" % (wdir,TAMO_file,i))


def merge_runs_ancestor(wdir, TAMO_file, i, target, genome):
	
	'''Merges the motifs within a single cluster.
	
	This function will identify motifs that are within the JASPER data set that
	are similar to the motifs within the cluster.
	'''

	# Merges the motifs in the same cluster using STAMP
	print "%s_sub_%s.tm" % (TAMO_file, i)
	tamo2tf("%s_sub_%s.tm" % (TAMO_file, i))
	
	# Gets the JASPER motifs that best match the motifs from within
	# the cluster.
	print "STAMP -tf  %s_sub_%s.tm.tf  -sd /home/chengzou/bin/STAMP/ScoreDists/JaspRand_PCC_SWU.scores -go  1000 -ge 1000 -cc PCC -align SWU -out %s_sub_%s.tm.tf_STout -chp > %s_sub_%s.tm.tf_STout.log" % (TAMO_file,i,TAMO_file,i,TAMO_file,i)
	os.system("STAMP -tf  %s_sub_%s.tm.tf  -sd /home/chengzou/bin/STAMP/ScoreDists/JaspRand_PCC_SWU.scores \
	 -go  1000 -ge 1000 -cc PCC -align SWU -out %s_sub_%s.tm.tf_STout -chp > %s_sub_%s.tm.tf_STout.log" % (TAMO_file,i,TAMO_file,i,TAMO_file,i))
	parse_out_STAMP_job(TAMO_file, i)
	
	# combines the JASPER motifs with the cluster motif and then
	# converts them all to one TAMO file
	os.system("cat  %s_sub_%s.tm.tf %s_sub_%s.tm.tf_SToutFBP.txt.mod %s_sub_%s.tm.tf_STout_tree_clusters.txt > %s_sub_%s_sum.tm.tf" % (TAMO_file,i,TAMO_file,i,TAMO_file,i,TAMO_file,i))
	tf2tamo("%s_sub_%s_sum.tm.tf" % (TAMO_file,i))

	# Gets the top motif within the TAMO file.
	pick_chunk_score(wdir, '%s_sub_%s_sum.tm.tf.tm' % (TAMO_file,i), target, genome)
	
	# Removes any files created in the processing.
	os.system("rm %s_sub_%s_sum.tm.tf.tm" % (TAMO_file,i))
	os.system("rm %s_sub_%s_sum.tm.tf"    % (TAMO_file,i))
	os.system("rm -R %s_sub_%s.tm.tf_ST*" % (TAMO_file,i))

def check_top_files(wdir):

	'''Check that each cluster has a TOP file'''

	# Gets all of the files in the working directory.
	wdir_file_list = os.listdir(wdir)
	
	# Find which clusters that have a command line in the merge_runs_cc and 
	# determines which job file this cluster coresponds to.
	merge_runs_cc_file = open('%s/merge_runs_cc' % wdir, 'r')
	
	# Dictionary contianing the cluster information
	# cluster_cc_dict = {cluster number: (job number, command line)}
	cluster_cc_dict = {}
	job_number = 1
	for line in merge_runs_cc_file:
		cluster = line.strip().split('-i')[1].split(' ')[1]
		cluster_cc_dict[cluster] = (str(job_number), line)
		job_number += 1
	
	
	# Identify which clusters have a TOP file.
	
	cluster_TOP_dict = {} # Dictionary containing clusters with TOP files
	
	for file_name in wdir_file_list:
		if file_name.endswith('.TOP'):
			
			cluster = file_name.strip().split('_sub_')[1].split('_')[0]
			
			cluster_TOP_dict[cluster] = 1
			
	
	print 'There are', len(cluster_cc_dict), 'clusters in merge_runs_cc,'
	print 'and', len(cluster_TOP_dict), 'of the clusters have a TOP file.'
		
	output = open('failed_merge_runs_cc', 'w')
	for cluster in cluster_cc_dict:
		if not cluster in cluster_TOP_dict:
			print 'job'+cluster_cc_dict[cluster][0]+'.sh'
			output.write(cluster_cc_dict[cluster][1])			
	output.close()

def main(): 
		
	# Get the function
	try:
		function = sys.argv[1]
	except IndexError:
		print __doc__
		sys.exit()
	
	# Set defaults
	TAMO_file = ''
	wdir      = os.path.abspath(os.curdir)
	distance  = 1
	ancestor  = 1
	height    = 20.0
	genome    = "/home/chengzou/project/Motif_run/_tamo_rawdata/all_on_array_withpromoter_seq.fa"
	mem       = '124'
	walltime  = 120
	write_failed_cc = False
	
	# Get arguments
	for n in range(2, len(sys.argv), 2):
		if sys.argv[n]   == '-w':
			wdir          = sys.argv[n+1].rstrip('/') 
		elif sys.argv[n] == '-t':
			TAMO_file     = sys.argv[n+1]
		elif sys.argv[n] == '-i':
			motif_index_file = sys.argv[n+1]
		elif sys.argv[n] == '-m':
			mapping_dir_root = sys.argv[n+1]
		elif sys.argv[n] == '-g':
			gene_list_file   = sys.argv[n+1]
		elif sys.argv[n] == '-wall':
			walltime         = sys.argv[n+1]
		elif sys.argv[n] == '-s':
			sim_file         = sys.argv[n+1]
		elif sys.argv[n] == '-mat':
			matrix           = sys.argv[n+1]
		elif sys.argv[n] == '-h':
			height           = sys.argv[n+1]
		elif sys.argv[n] == '-mem':
			mem              = sys.argv[n+1]
		elif sys.argv[n] == '-c':
			clustering_file  = sys.argv[n+1]
		elif sys.argv[n] == '-a':
			ancestor         = sys.argv[n+1]
		elif sys.argv[n] == '-target':
			target           = sys.argv[n+1]
		elif sys.argv[n] == '-genome':
			genome           = sys.argv[n+1]
		
		
	

		else:
			print sys.argv[n], 'is an unknown flag.' 
			
			cont = raw_input('Do you wish to procede (y/n/help)')
			if cont == 'y':
				pass
			elif cont == 'n':
				sys.exit()
			elif cont == 'help':
				print __doc__
	
	if function == 'sim_matrix':
		sim_matrix(motif_index_file, mapping_dir_root, gene_list_file, TAMO_file)
	elif function == 'run_as_job':
		run_as_job(motif_index_file, mapping_dir_root, gene_list_file, TAMO_file, walltime, mem)
	elif function == 'create_matrix':
		create_matrix(sim_file)
	elif function == 'remove_duplicate':
		remove_duplicate_motifs(TAMO_file)
	elif function == 'mapped_TAMO':
		TAMO_of_mapped_motifs(TAMO_file, motif_index_file)
	elif function == 'run_UPGMA':
		run_UPGMA(matrix, height, wdir, walltime = 120, mem = '124')
	elif function == 'merge_runs_cc':
		merge_runs_cc(TAMO_file, wdir, ancestor, target, genome, clustering_file)
		

	elif function == 'help':
		print '/'.join(os.path.abspath(__file__).split('/')[:-1])
		print __doc__
		sys.exit()
		
	else:

		cont = raw_input('This is an unknown function, would you like some\
 help? (y/n)')
		if cont == 'y':
			print __doc__
			sys.exit()
		elif cont == 'n':
			print "Goodbye"
			sys.exit()
		else:
			print __doc__
			print '\nThought you would like to see the help, anyways.'
			sys.exit()

if __name__ == "__main__":
	main()

