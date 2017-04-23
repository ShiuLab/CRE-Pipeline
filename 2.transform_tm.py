from   TAMO.MotifTools import Motif, print_motifs,save_motifs
from   TAMO.MotifTools import Motif_from_text
import sys,os,re

inp=sys.argv[1]
list=[]
for line in open(inp,"r"):
	m = Motif_from_text(line.strip())
	list.append(m)
save_motifs(list,"%s.tm" % inp)
