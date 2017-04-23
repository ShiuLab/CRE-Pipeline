import os,sys
import re
dir=sys.argv[1].rstrip("/")
des_dir=sys.argv[2].rstrip("/")

os.system("cd %s" % dir)
os.chdir(dir)
dic={}
for file in os.listdir(dir):
	if file.endswith(".mix"):
		for line in open(file,"r"):
			
			if 1:
				m=re.search("\d",line[line.find(")")+1:])
				#print m
				end=m.start()
				motif=line[line.find(")")+1:line.find(")")+1+end].strip()
				#print line,motif
				ML=len(motif)
				if not dic.has_key(ML):
					dic[ML]={}
				dic[ML][motif]="1"
			"""
			L=line.split()
			motif=L[0]
			ML=len(motif)
			if not dic.has_key(ML):
				dic[ML]={}
			dic[ML][motif]="1"
			"""
			

for i in dic.keys():
	oup=open("%s/all_Weeder_out_%s" % (des_dir,i),"w")
	oup.write("%s\n" % "\n".join(dic[i].keys()))
	oup.close()
	os.system("python /mnt/home/seddonal/scripts/5_motif_merging/2.transform_tm.py \
		%s/all_Weeder_out_%s" % (des_dir,i))
