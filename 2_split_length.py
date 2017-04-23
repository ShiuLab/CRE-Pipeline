import sys, re, os, math, getopt


from TAMO.Clustering.MotifCompare import *
from TAMO.Clustering              import MotifCompare
from TAMO                         import MotifTools
from   TAMO.MotifTools import Motif, print_motifs,save_motifs



tm_file=sys.argv[1]
motifs=  MotifTools.txt2motifs(tm_file)
dic={}
for m in motifs:
	if not dic.has_key(len(m.oneletter)):
		dic[len(m.oneletter)]=[]
	dic[len(m.oneletter)].append(m)
print len(motifs),len(dic.keys())
for i in dic.keys():
	file = tm_file[:-6] + str(i) + ".out.tm"
	print file
	save_motifs(dic[i],file)
