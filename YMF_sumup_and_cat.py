import os,sys
import re
"""


"""


dir=sys.argv[1].rstrip("/")
des_dir=sys.argv[2].rstrip("/")
prefix=sys.argv[3]

os.system("cd %s" % dir)
os.chdir(dir)


dic={}
for file in os.listdir(dir):
	if not file.endswith("fa") and file.startswith(prefix):
		for line in open(file,"r"):
			if not line.startswith("The"):
				m=re.search("\d",line)
				end=m.start()
				motif=line[:end].strip()
				#print line,motif
				ML=len(motif)
				if not dic.has_key(ML):
					dic[ML]={}
				dic[ML][motif]="1"

for i in dic.keys():
	oup=open("%s/all_YMF_out_%s" % (des_dir,i),"w")
	oup.write("%s\n" % "\n".join(dic[i].keys()))
	oup.close()

	os.system("python 2.transform_tm.py \
		%s/all_YMF_out_%s" % (des_dir,i))
