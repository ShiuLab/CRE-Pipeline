"""
likelihood score use as MAP, since the bigger the likeliscore,
the Motif is better.

"""

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
        

def MS_parse(file):
        'Parse MDscan file'
        motifs=[]
        inp=open(file,"r")
        lines=inp.readlines()
        alloutput = '\n'.join(lines)
        premotifs = alloutput.split('\n#id:')
        #print 'lens',len(premotifs)
        for pm in premotifs:
            #print "pm",pm
            sublines = pm.split('\n')
            score, seednum = 0,0
            seqs = []
            for line in sublines:
                if line.strip().startswith("box"):
                    toks    = line.split()
                    #print 'toks',toks
                    score   = float(toks[-1])
                    seednum = int(toks[-7])
                elif not line.startswith("#") and line:
                    #print line.split()
                    seqs.append(line.split()[-1][1:-2])
            #print "SEQS: ",seqs
            if seqs:
                m = Motif(seqs)
                m.MAP = score
                m.seednum = seednum
                motifs.append(m)
        nmotifs = len(motifs)
        save_motifs(motifs,'%s.tm' % file)

MS_parse(inp)

