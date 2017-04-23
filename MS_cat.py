import os,sys
dir=sys.argv[1]
des_dir=sys.argv[2]
dir=dir.rstrip("/")
des_dir=des_dir.rstrip("/")
for j in range(6,16):
	print j
	print ("find %s -name \*MS-%s.out.tm -print0 | sort -z | xargs -0 -n2 cat > %s/MEME_%s.out.tm" % (dir,j,des_dir,j))
	os.system("find %s -name \*MS-%s.out.tm -print0 | sort -z | xargs -0 -n2 cat > %s/MS_%s.out.tm" % (dir,j,des_dir,j))
