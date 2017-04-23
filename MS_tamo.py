import os
import sys

directory = sys.argv[1]

if not  directory.endswith("/"):
	directory= directory + "/"

def isoutfile(filename):
    if "." in filename:
        if filename[filename.rindex('.'):] == '.out':
            return True
        else:
            return False
    else:
	return False

for fil in filter(isoutfile, os.listdir(directory)):
    os.system("python 1.4_MotifSampler.py "+directory+fil)
