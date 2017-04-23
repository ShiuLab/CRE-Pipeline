import os,sys

# Written by Alex Seddon
# You have to install MDScan in a dir you have write access to
# 06/13/12 Modified by Shiu. 

files_path = sys.argv[1].rstrip("/")
mdscan     = sys.argv[2] # Need to specify path + name of executable
out_dir    = sys.argv[3] # Output directory
os.system("mkdir %s" % out_dir)

oup=open("commands_MD" ,"w" )
for file in os.listdir(files_path):
	if file.endswith(".fa"):
		for W in range(6,16):
			oup.write("%s -i %s/%s -s 30 -r 5 -w %s -o %s/%s-%s-md.out \n" % \
							(mdscan,files_path,file,W,out_dir,file[:-3],W))

oup.close()
