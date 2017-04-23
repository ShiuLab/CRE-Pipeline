import os,sys
dir=sys.argv[1]
des_dir=sys.argv[2]
dir=dir.rstrip("/")
des_dir=des_dir.rstrip("/")
prefix=sys.argv[3]

os.system("cd %s" % dir)
os.chdir(dir)

##############file dic
file_dic={}	
for file in os.listdir(dir):
	if file.endswith(".tm") and file.startswith(prefix):
		#print file
		len=int(file[file.rfind("-")+1:file.rfind(".tm")])
		run=int(file[file.find("-")+1:file.find("-",file.find("-")+1)])
		#if len<=15 and run<=3:
		if len<=15 and run<=1:
			if not file_dic.has_key(len):
				file_dic[len]={}
			file_dic[len][file]="1"
		
	
#for j in range(6,16):
for j in range(6,16):
	print j
	if file_dic.has_key(j):
		oup=open("AA_%s.out.tm" % j,"w")
		oup.close()
		for file in file_dic[j].keys():
			#print j,file
			os.system("cat %s>> AA_%s.out.tm" % (file,j))
		os.system("mv AA_%s.out.tm %s/AA_%s.out.tm" % (j,des_dir,j))
