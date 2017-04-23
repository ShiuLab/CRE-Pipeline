import os,sys

files_path =sys.argv[1].rstrip("/")
program_dir=sys.argv[2] 
#program_dir='/home/chengzou/bin/Weeder1.3/weederlauncher.out'

# See if user has specified the species code
try:
    species = sys.argv[3]
except IndexError:
    species = 'AT' # sys.argv[3] is empty, then set arabidopsis code as default.

oup=open("commands_Wd" ,"w" )

for file in os.listdir(files_path):
	if file.endswith(".fa"):
		oup.write("%s %s/%s %s medium M S T10 \n" \
% (program_dir,files_path,file,species))

oup.close()
