from TAMO import MotifTools
from   TAMO.MotifTools import Motif, save_motifs
import sys,os,re,math


if __name__ == "__main__":
	file=sys.argv[1]
	by=int(sys.argv[2]) # how many motifs per file
	ml=MotifTools.txt2motifs(file)
	
	total=len(ml)/by
	for i in range(total):
		print i
		print i*by+by,file+'_n%s' % i
		save_motifs(ml[i*by:i*by+by],file+'_n%s' % i)
	print 	total*by,len(ml),file+'_n%s' % (total)
	save_motifs(ml[total*by:len(ml)],file+'_n%s' % (total))



