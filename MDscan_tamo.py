import os
import sys

directory = sys.argv[1]


if not directory.endswith("/"):
	directory=directory+"/"


def isoutfile(filename):
    if "." in filename:
        if filename[filename.rindex('.'):] == '.out':
            return True
        else:
            return False
    else:
            return False


for fil in filter(isoutfile, os.listdir(directory)):
    print fil
    #os.system("/home/zou/bin/Python-2.5.2/python /home/zou/project/cis_cold/_script/1.3_MD_scan.py "+directory+fil)
    os.system("python /mnt/home/seddonal/scripts/5_motif_merging/1.3_MD_scan.py "+directory+fil)
