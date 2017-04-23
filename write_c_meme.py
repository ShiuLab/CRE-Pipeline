import os,sys

### You have to install meme under /mnt/home/%s/bin/meme
###

files_path = sys.argv[1]
files_path = files_path.rstrip("/")
meme_path  = sys.argv[2] # /mnt/home/chengzou/bin/meme/bin/meme
maxsize    = sys.argv[3]

oup=open("commands_meme_%s" % files_path[files_path.rfind("/")+1:] ,"w" )

os.system("cd %s\n" % files_path)
os.chdir(files_path)
for file in os.listdir(files_path):
	if file.endswith(".fa"):
		oup.write("%s %s/%s \
-dna -mod anr -nmotifs 10 -nsites 100 -minw 6 -maxw 15  -maxiter 20 -maxsize %s -revcomp -text \
>%s/%s.out\n" % (meme_path,files_path,file,maxsize,files_path,file[:-3]))

oup.close()
