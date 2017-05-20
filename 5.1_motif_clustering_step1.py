########################################################
# The original Motif clustering requirs lots of memery.
# Here I seperated them by the one letter representation
# 080801:too long list to cat,write one by one
########################################################
"""
~/bin/Python-2.5.2/python /home/zou/project/cis_cold/4_TAMO_test/_script/5.1_motif_clustering_step1.py
 -m all_15.out.tm -d 0.03964 --dfunc  entropyrange
"""



import sys, re, os, math, getopt,time
GLOBALS = {}

from TAMO.Clustering.MotifCompare import *
from TAMO.Clustering              import MotifCompare
from TAMO                         import MotifTools
from   TAMO.MotifTools import Motif, save_motifs


def check_running(user):
	while 1:
		print "1"
		done=1
		os.system("qstat > TMP_check_running_log")
		inp = open("TMP_check_running_log")
		for i in inp.readlines():
			if user in i and "tmp" in i:
				done=0
		print done
		os.system("rm TMP_check_running_log")
		if done==1:
			break
		time.sleep(30)
	return "done"

def main():
    motifs = getarg('motifs')
    print 'There are', len(motifs), "motifs in the TAMO file."
    #print type(motifs) #motifs is a list
    #build a dic by the one letter
    dic_one={}
    
    # Place all motifs with the same one letter code into the same dictionary.
    for i in range(0,len(motifs)):
        oneletter=motifs[i].oneletter.upper()
        if oneletter.startswith("."):
            oneletter="N" + oneletter[1:]
        if not dic_one.has_key(oneletter):
            dic_one[oneletter]=[]
        dic_one[oneletter].append(i)
    print len(motifs)
    
    # Cluster all the motifs that have the same one-letter code.
    command_txt=open("clustering_command","w")
    ori_motifs=[]
    for key in dic_one.keys():
        print dic_one[key]
        if len(dic_one[key])>1:
            mlist=[]
            for j in dic_one[key]:
                mlist.append(motifs[j])
            #print mlist
            if len(dic_one[key])>1000:
                save_motifs(mlist,"%s.tmp.big" % key)
                save_motifs([mlist[0]],"%s.tmp.big-cl" % key)
            else:
                save_motifs(mlist,"%s.tmp" % key)
                #print "/home/zou/bin/Python-2.5.2/python /home/zou/lib/python2.5/site-packages/TAMO/Clustering/UPGMA.py  -m  %s.tmp -d  %s  --dfunc %s >log\n" % (key,DMAX,DFUNCtxt)
                command_txt.write("module load TAMO; module load motility; python UPGMA.py  -m  %s/%s.tmp -d  %s  --dfunc %s >log\n" % (directory, key,DMAX,DFUNCtxt))
        else:
            for j in dic_one[key]:
               ori_motifs.append(motifs[j])
    save_motifs(ori_motifs,"ori_motifs")
    command_txt.close()
    os.system("python qsub_hpc.py -f quene -c clustering_command -u %s -n 230" % user)
    ###some times there are too many "-cl" 080801
    oup=open("after_cluster.tm","w")
    if check_running(user)=="done":
        for file in os.listdir(os.getcwd()):
            if file.endswith("-cl"):
                for line in open(file,"r"):
                    oup.write("%s" % line)
                #os.system("rm %s" % file)
        oup.close()

    os.system("cat ori_motifs after_cluster.tm > After_cl_%s_d%s" % (getarg('file'),DMAX))
    os.system("rm ori_motifs")
    os.system("rm after_cluster.tm")
    #os.system("rm clustering_command")
    #os.system("rm -R *.tmp")
    os.system("rm *.big")
    for i in range(0,10):
        os.system("rm -R tmp%s*" % i)

def set_dfunc(_dfunc, _dmax):
    '''
    set_dfunc(dfunc, dmax)

    Set the distance/divergence/difference metric and threshold, and synchronize
    between this module and MotifCompare

    Examples:
    set_dfunc(MotifCompare.negcommonbitsrange, -8.0)
    set_dfunc(MotifCompare.diffrange, 0.23)
    '''

    global DFUNC, DMAX
    DFUNC = _dfunc
    MotifCompare.DFUNC = _dfunc
    DMAX  = _dmax


def parse_opts():
    global GLOBALS
    global DFUNCtxt, DMAX, user, directory
    short_opts = 'm:d:'
    long_opts  = ['dfunc=', 'user=', 'dir=']
    try:   opts, args = getopt.getopt(sys.argv[1:], short_opts, long_opts)
    except getopt.GetoptError:
        print getopt.GetoptError.__dict__
        usage()
    if not opts: usage()
    print opts
    GLOBALS['args'] = args
    GLOBALS['motifs'] = []
    DFUNCtxt = None
    for opt,value in opts:
        if opt == '-m':                   GLOBALS['motifs'] = MotifTools.txt2motifs(value)
        if opt == '--dfunc':              DFUNCtxt          = value
        if opt == '-d':                   DMAX              = float(value)
        if opt == '-m':                   GLOBALS['file']   = value
        if opt == '-d':                   GLOBALS['dis']    = value
        if opt == '--user':               user              = value
        if opt == '--dir':                directory         = value.rstrip()

    # Deal with DFUNC and DMAX
    if DFUNCtxt == 'NCB':
        #_DFUNC = MotifCompare.negcommonbits
        _DFUNC = MotifCompare.negcommonbitsrange
    elif DFUNCtxt:
        try:
            exec ("_DFUNC = MotifCompare.%s"%DFUNCtxt)
        except:
            usage("No such distance metric: %s"%DFUNCtxt)
    #if _DFUNC:  set_dfunc(_DFUNC,DMAX)


def getarg(varname):
    global GLOBALS
    if GLOBALS.has_key(varname):   return GLOBALS[varname]
    else:                          return None


if __name__ == '__main__': 
    parse_opts()
    #file=getarg('file')
    #os.system("cd %s" % file[:file.rfind("/")])
    #os.chidir("%s" % file[:file.rfind("/")])
    print "#" + ' '.join([x.replace(' ','\ ') for x in sys.argv])
    main()                                  



