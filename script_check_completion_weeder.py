#
# Check weeder run completion
#

import sys,os

run_dir  = sys.argv[1]
cmd_file = sys.argv[2]

cmd = open(cmd_file).readlines()
oup = open(cmd_file+".RERUN","w")
for i in cmd:
	L  = i.split(" ")
	fa = ""
	for j in L:
		if j[-3:] == ".fa":
			fa = j
			break
	print fa
	ofiles = [0,0,0]
	if os.path.isfile("%s.html" % (fa)):
		ofiles[0] = 1
	if os.path.isfile("%s.mix" % (fa)):
		ofiles[1] = 1
	if os.path.isfile("%s.wee" % (fa)):
		ofiles[2] = 1
	
	if 0 in ofiles:
		print " incomplete",ofiles
		oup.write(i)
