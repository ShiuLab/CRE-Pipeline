from   TAMO.MotifTools import Motif, print_motifs,save_motifs
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
        
def Meme_out(file):
        'Parse MEME file'
        motifs=[]
        inp=open(file,"r")
        lines=inp.readlines()
        fastafile=''
        probes=[]
        i = 0
        num = 0
        background = {}
        for i in range(len(lines)):
            toks = lines[i].split()
            #print toks
            for tok in toks: tok.strip()
            if len(toks) == 2 and not fastafile and toks[0] == 'DATAFILE=':
                fastafile = toks[1]
            if len(toks) < 4: continue
            if toks[0] == 'Background' and toks[2] == 'frequencies':
                toks = lines[i+1].split()
                for j in [0,2,4,6]:
                    background[toks[j]] = float(toks[j+1])
            #print lines[i],
            if toks[0] == "MOTIF":
                m = Motif([],background)
                num    = int(toks[1])
                m.evalue = float(toks[13])
            elif (toks[0] == 'Motif' and int(toks[1]) == num and
                  toks[3] == 'BLOCKS'):
                offset = 3
                j      = 0
                while lines[i+offset+j][0:2] != '//':
                    #print lines[i+offset+j],
                    _toks = lines[i+offset+j].split()
                    seq = _toks[-2]
                    seq = re.sub('X','N',seq)
                    m.seqs.append(seq)
                    j = j + 1                
            elif (toks[0] == 'Motif' and int(toks[1]) == num and
                  toks[3] == 'probability'):
                offset = 3
                j      = 0
                counts = []
                while lines[i+offset+j][0:10] != '------------'[0:10]:
                    probs = lines[i+offset+j].split()
                    _d = {}
                    for key,val in zip(['A', 'C', 'G', 'T'], probs):
                        _d[key] = float(val)
                    counts.append(_d)
                    j = j + 1
                #for c in counts:
                    #print c
                m.compute_from_counts(counts,0.00001)
                motifs.append(m)
            elif (toks[0] == 'Sequence' and toks[1] == 'name' and
                  toks[2] == 'Weight'):
                offset = 2
                j      = 0
                while lines[i+offset+j][0:5] != '*****':
                    _toks = lines[i+offset+j].split()
                    probes.append(_toks[0])
                    if len(_toks) > 5:
                        probes.append(_toks[3])
                    j = j + 1
        #writer = MyWriter(sys.stdout, '%s.tm' % file)
        #sys.stdout = writer
        save_motifs(motifs,'%s.tm' % file)
if __name__=="__main__":
	file=sys.argv[1]
	Meme_out(file)
