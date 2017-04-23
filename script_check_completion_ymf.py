
#
# Check completion of YMF runs, requires:
# fa_dir as sys.argv[1] : Fasta directory specified in batch YMF run
# ou_dir as sys.argv[2] : YMF output directory specified in batch YMF run
# 

import sys, os

def find_not_done(fa_dir,ou_dir):
	wlist = ["6","7","8","9","10","10-spacer"]
	flist = os.listdir(fa_dir)
	not_done = {} # {fasta:[w]}
	for i in flist:
		OI = "%s/%s" % (ou_dir,i) # ith output dir
		print i
		try:
			fOI = os.listdir(OI)    # files in ith output dir
			for j in wlist:
				fOJ = os.listdir("%s/%s" % (OI,j))
				if "results" not in fOJ:
					print " not done:",j
					if i not in not_done:
						not_done[i] = [j]
					else:
						not_done[i].append(j)
		except OSError:
			print " no_output"
			not_done[i] = wlist
	return not_done
	
def fasta_to_dict(fa):
	# Read the fasta file, join all lines, then split with ">". 1st element empty
	flines = "".join(open(fa).readlines()).split(">")[1:]
	fdict = {} # {name:seq}
	for i in flines:
		nl = "\n"
		if i.find('\r') != -1:
			if i.find('\n') != -1:
				nl = "\r\n"
			else:
				nl = "\r"
		i = i.split(nl) # 1st element is name, following are sequences
		name = i[0]
		seq  = "".join(i[1:])
		fdict[name] = seq
	return fdict
			
#-------------------------------------------------------------------------------

fa_dir   = sys.argv[1] # Fasta directory specified in batch YMF run
ou_dir   = sys.argv[2] # YMF output directory specified in batch YMF run
not_done = find_not_done(fa_dir,ou_dir) # {fasta_name:[window_sizes]}

# Check if the FASTA file contains N, if so, modify the sequences, if not, just
# copy them over.
incomp_fasta = "_seq_ymf_incomplete"

if not os.path.isdir(incomp_fasta):
	os.system("mkdir %s" % incomp_fasta)
	
for i in not_done:
	print i
	fd = fasta_to_dict("%s/%s" % (fa_dir,i))
	oup = open("%s/%s" % (incomp_fasta,i),"w")
	for j in fd:
		if fd[j].find("N") == -1:
			oup.write(">%s\n%s\n" % (j,fd[j]))
		else:
			print " found N:",j
		
oup = open("fasta_config","w")
for i in not_done:
	oup.write("%s\t%s\n" % (i,",".join(not_done[i])))
oup.close()

