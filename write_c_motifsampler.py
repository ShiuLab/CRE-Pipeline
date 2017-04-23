import os,sys

### You have to install meme under /mnt/home/%s/bin/meme
###

files_path =sys.argv[1].rstrip("/")
#program=sys.argv[2] #  /mnt/home/chengzou/bin/MDscan/MDscan.linux
program='/home/chengzou/bin/MotifSamplerx86-64'
oup=open("commands_MS" ,"w" )
for file in os.listdir(files_path):
	if file.endswith(".fa"):
		for W in range(6,16):
			oup.write("%s -f %s/%s -b /mnt/home/chengzou/bin/AT_bg_3order -w %s -n 10 > %s/%s.MS-%s.out \n" 
% (program,files_path,file,W,files_path,file[:-3],W))

oup.close()