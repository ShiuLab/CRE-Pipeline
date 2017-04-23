
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
        

def MD_parse(file):
        'Parse MDscan file'
        
        motifs    = []
        inp       = open(file,"r")
        lines     = inp.readlines()
        alloutput = '\n'.join(lines)
        premotifs = alloutput.split('\nMotif ')  # Creates a list of all of the
        										 # info for each motif.
        
        # This loop parses the infor for each individual motif in premotifs
        for pm in premotifs:
        	
        	# Skips the header, or adds 
            if not pm.startswith("*"):
                pm="Motif "+pm
                #print "pm",pm
                
            sublines = pm.split('\n') # Splits the pm into its component lines
            score, seednum = 0,0
            seqs = []

			
			# This loop looks at each individual line in sublines.
            j = 0
            while j < len(sublines):
            
                line = sublines[j]
                if line.startswith("Motif"):
                    toks    = line.split()
                    #print "toks",toks
                    score   = float(toks[5][:-1])
                    seednum = int(toks[7][:-1])
                
                if line.find('>') == 0:
                    j+=2
                    line=sublines[j]
                    #print "add_line",score,line
                    seqs.append(line)
                j+=1
            #print "SEQS: ",seqs
            if seqs:
                m = Motif(seqs)
                m.MAP = score
                m.seednum = seednum
                motifs.append(m)
        nmotifs = len(motifs)
        save_motifs(motifs,'%s.tm' % file)
MD_parse(inp)

