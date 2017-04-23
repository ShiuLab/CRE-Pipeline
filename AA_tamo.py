import os
import sys

directory = sys.argv[1]
if not directory.endswith("/"):
	directory= directory+"/"
def isoutfile(filename):
    if filename.endswith('.out'):
        return True
    else:
        return False

for fil in filter(isoutfile, os.listdir(directory)):
    os.system("python /mnt/home/seddonal/scripts/5_motif_merging/1.1_AA2tamo.py "+directory+fil)
    
