import sys, re, os, math, getopt
GLOBALS = {}

from TAMO.Clustering.MotifCompare import *
from TAMO.Clustering              import MotifCompare
from TAMO                         import MotifTools
from   TAMO.MotifTools import Motif, print_motifs


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
        
    
def main():
    os.system("rm %s.dm" % getarg('file'))
    motifs = getarg('motifs')
    oup=open("%s.dm" % getarg('file'),"w")
    motifs = getarg('motifs')
    n=0
    #Dmat = computeDmat(motifs,VERBOSE=1,DFUNC=DFUNC)
    #oup=open("%s.dm" % getarg('file'),"w")

    for i in range(len(motifs)):
        print i      
        for j in range(0,i):
        	oup.write("-\t")
        
        Dmat = computeDpairs(motifs[i+1:len(motifs)],[motifs[i]],VERBOSE=0,DFUNC=DFUNC)
        oup.write('0\t%s\n' % ("\t".join(map(lambda x: "%0.2f" % x,Dmat[0].values()))))
        n+=1
    oup.close()

# This function is not in the normal TAMO distribution on HPCC. I found it on
# calculon in /home/zou/lib/python2.5/site-packages/TAMO/Clustering/MotifCompare.py.
# My guess is that it is a custom written script by Cheng, because I cannot 
# find any evidence of the function in any other copy of TAMO or its documentation.
# Alex Seddon 27 Jun 2012

def pccrange(self,other,Srange,Orange):
    from numpy import array
    from scipy.stats import pearsonr
    #import statistics
    '''Utility function: compute diff of '''
    '''self and other from self+Sstart, width of other'''
    POW     = math.pow
    Dtot    = 0
    divisor = float(len(Srange))
    
    #print Srange,Orange
    #print len(Srange),len(Orange)
    '''Computes distance'''
    pcc=0
    
    for si, oi in zip(Srange,Orange):
        Sl=[]
        Ol=[]
        for L in ACGT:
            Sl.append(POW(2,self.logP[si][L]))
            Ol.append(POW(2,other.logP[oi][L]))
        #print Sl,Ol,statistics.correlation(Sl, Ol, method = "Pearson")
        pcc+= pearsonr(Sl, Ol)[0]
    _dist=1-pcc/divisor
    """
    for L in ACGT:
        Sl=[]
        Ol=[]
        for si, oi in zip(Srange,Orange):
            Sl.append(POW(2,self.logP[si][L]))
            Ol.append(POW(2,other.logP[oi][L]))
        print Sl,Ol,statistics.correlation(Sl, Ol, method = "Pearson")
        pcc+= statistics.correlation(Sl, Ol, method = "Pearson")
    _dist=pcc/divisor

        D = 0
        for L in ACGT:
            D = D + POW( POW(2,self.logP[si][L]) - POW(2,other.logP[oi][L]), 2 )
        col_dist  = math.sqrt(D)/math.sqrt(2.0)
        Dtot      = Dtot + col_dist
    _dist= Dtot/divisor
    """
    #DBG
    #print self[min(Srange),max(Srange)+1], other[min(Orange),max(Orange)+1], _dist
    #print self.oneletter[min(Srange):max(Srange)+1]
    #print other.oneletter[min(Orange):max(Orange)+1], '%6.4f'%_dist
    #print 
    return _dist

# This is another custom function written by Cheng and found in the same place
# as pccrange. Alex Seddon, 

def computeDpairs(motifs1,motifs2,VERBOSE=0,DFUNC=DFUNC):
    #print motifs1,motifs2
    N = len(motifs1)
    N2 =len(motifs2)
    dmat = {}
   
    #if VERBOSE: print "              |%s|"%('-'*(N-1))  #Pretty status bar
    #if VERBOSE: print "Computing ...  ",
    Nhalf = int(.293*N)
    for i in range(N2):
        dmat[i]={}
        #print i
        if VERBOSE:
            if i == Nhalf: sys.stdout.write('|'); sys.stdout.flush()
            else:          sys.stdout.write('.'); sys.stdout.flush()
        A = motifs2[i]
        #print A,range(N2)
        for j in range(N):
            B = motifs1[j]
            #print "B", B
            #D = minshortestoverhangdiff(A,B,OVLP(A,B),DFUNC=DFUNC)
            D=minshortestoverhangdiff(A,B,5,DFUNC=DFUNC) # use the minimal overlap=6,which is default
            dmat[i][j] = D              #Symmetric
            #print D
    if VERBOSE: print                               #End of Pretty status bar
    return dmat
    
def parse_opts():
    global GLOBALS
    global DFUNC
    short_opts = 'i:'
    long_opts = ['dfunc=']
    try:   opts, args = getopt.getopt(sys.argv[1:], short_opts,long_opts)
    except getopt.GetoptError:
        print getopt.GetoptError.__dict__
        usage()
    if not opts: usage()
    print opts
    GLOBALS['args'] = args
    GLOBALS['motifs'] = []
    DFUNCtxt = None
    for opt,value in opts:
        if opt == '-i':                   GLOBALS['motifs'] = MotifTools.txt2motifs(value)
        if opt == '-i':                   GLOBALS['file'] = value
        if opt == '--dfunc':              DFUNCtxt = value


    # Deal with DFUNC and DMAX
    if DFUNCtxt == 'pccrange':
        DFUNC = pccrange
    else:
        if DFUNCtxt == 'NCB':
            #_DFUNC = MotifCompare.negcommonbits
            _DFUNC = MotifCompare.negcommonbitsrange
        elif DFUNCtxt:
            try:
                exec ("_DFUNC = MotifCompare.%s"%DFUNCtxt)
            except:
                usage("No such distance metric: %s"%DFUNCtxt)
        if _DFUNC:  set_dfunc(_DFUNC)

def usage(txt=''):
    '''
    Place information about command-line arguments here
    '''
    if txt: print "Error: %s"%txt
    print 'Usage: %s -m motifs'%(
        re.sub('^.*/','',sys.argv[0]))
    print ""

    print '         [-d threshold]        Distance threshold.  Default is 0.2'
    print '         [--dfunc <function>]  Distance metric.  Examples are NCB, and diffrange.  Any'
    print '                               "range" function in MotifCompare.py is acceptible'
    print ''
    print 'Some useful examples:'
    print ''
    print '   UPGMA.py -m <motiffile.tamo> -d  0.2  --dfunc diffrange (default)'
    print '   UPGMA.py -m <motiffile.tamo> -d -8.0 --dfunc NCB       (Count common bits, negated for minimization)'

    
    sys.exit(1)

       
def getarg(varname):
    global GLOBALS
    if GLOBALS.has_key(varname):   return GLOBALS[varname]
    else:                          return None


def set_dfunc(_dfunc,):
    '''
    set_dfunc(dfunc, dmax)

    Set the distance/divergence/difference metric and threshold, and synchronize
    between this module and MotifCompare

    Examples:
    set_dfunc(MotifCompare.negcommonbitsrange, -8.0)
    set_dfunc(MotifCompare.diffrange, 0.23)
    '''
    global DFUNC
    DFUNC = _dfunc
    MotifCompare.DFUNC = _dfunc



   
if __name__ == '__main__': 
    parse_opts()
    main()  
