
from   TAMO.MotifTools import Motif, save_motifs
import sys,os,re

inp=sys.argv[1]
os.system("rm %s.tm" % inp)
class MyWriter:

    def __init__(self, stdout, filename):
        self.stdout = stdout
        self.logfile = file(filename, 'a')

    def write(self, text):
        self.stdout.write(text)
        self.logfile.write(text)

    def close(self):
        self.stdout.close()
        self.logfile.close()
    def flush(self):
        self.stdout.flush()
        

def AA_out(file):
	motifs=[]
	inp=open(file,"r")
	lines=inp.readlines()
	i = 0
	while i < len(lines):
		line = lines[i]
		toks = line.split()
		#print toks
		#if '-i' in toks and not self.fastafile:
		#	idx = toks.index('-i')
		#	self.fastafile = toks[idx+1]
		if (len(toks) > 0 and toks[0] == 'Motif'):
			seqs = []
			while 1:
				i = i + 1
				line = lines[i]
				toks = line.split()
				if toks[0][0] == '*' or toks[0].startswith("-"): break
				seqs.append(toks[0])
			#print "input_seq",seqs
			M = Motif(seqs)
			i = i + 1
			line = lines[i]
			toks = line.split()
			MAP = float(toks[2])
			if MAP > 1000: MAP = 0  #likely to be an AlignACE error
			M.MAP = MAP
			#print "M",MAP
			#self.motifs.append(M)
			motifs.append(M)
			#print motifs
		i = i + 1
	nmotifs = len(motifs)
	dic={}
	for i in motifs:
		lenn=len(i.oneletter)
		if not dic.has_key(lenn):
			dic[lenn]=[]
		dic[lenn].append(i)
		#print dic
	for i in dic.keys():
		#print dic[i]
		save_motifs(dic[i],'%s-%s.tm' % (file,i))
		
AA_out(inp)
