# This script will extract the motif name from each file and then output it a
# 2col file with the file name in the first file and the motif name in the
# second line. The motif name will be modified so that . are replaced with N.

import sys, os

map_dir = sys.argv[1] # the name of the directory with the mapping files
output_motif_list  = open(sys.argv[2], 'w') # Open an output file.
output_mapping_index = open(sys.argv[2]+'.index', 'w')

map_files = os.listdir(map_dir) # The list of mapping files.


# This motif 
for map_file in map_files:
	if map_file.endswith('.pvalue'):
		map = open(map_dir+map_file, 'r') # Opens the map file
		try:
			motif = map.readline().strip().split('\t')[1] 
		except IndexError:
			print map_file
			
			
		else:
			# The following loop 
			tmp_motif = ''
			for base in motif:
				if base == '.':
					base = 'N'
				tmp_motif = tmp_motif+base
			motif = tmp_motif
			output_motif_list.write(motif+'\n')
			output_mapping_index.write(map_file+'\t'+motif+'\n')

		map.close()

output_motif_list.close()
output_mapping_index.close()

