import os
import sys

directory = sys.argv[1]


def isoutfile(filename):
    if "." in filename:
        if filename[filename.rindex('.'):] == '.out':
            return True
        else:
            return False
    else:
	return False


for fil in filter(isoutfile, os.listdir(directory)):
    os.system("python /mnt/home/seddonal/scripts/5_motif_merging/1.2_MEME2tamo.py "+directory+fil)
    os.system("python /mnt/home/seddonal/scripts/5_motif_merging/2_split_length.py "+directory+fil+".tm")
