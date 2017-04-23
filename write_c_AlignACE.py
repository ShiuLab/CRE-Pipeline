import os,sys



### You have to install  under /mnt/home/chengzou/bin/alignace2004/AlignACE
###

files_path =sys.argv[1].rstrip("/")
#program_dir=sys.argv[2] #  /mnt/home/chengzou/bin/alignace2004/AlignACE
program_dir="/home/chengzou/bin/alignace2004/AlignACE"
oup=open("commands_AA" ,"w" )
for file in os.listdir(files_path):
	if file.endswith(".fa"):
		for W in range(1,11):
			for L in range(6,16):
				oup.write("%s -i %s/%s -seed %s -numcols %s -gcback 0.333 >%s/%s-%s-L%s-AA.out \n" 
% (program_dir,files_path,file,W,L,files_path,file[:-3],W,L))

oup.close()