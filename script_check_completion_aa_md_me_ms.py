
#
# Check completion
#

import sys,os

if len(sys.argv) != 4:
	print "script cmd_file output_dir outfile_index_in_cmd"
	print "Quit!\n"
	sys.exit(0)

cmd = open(sys.argv[1]).readlines() # command lines
odir= sys.argv[2]        # output directory
otok= int(sys.argv[3])   # output token in the cmd line, zero indexing

oup = open("%s.RERUN" % sys.argv[1],"w")
c   = 0
for i in cmd:
  oname = i.strip().split(" ")[otok].split("/")[-1]
  if not os.path.isfile("%s/%s" % (odir,oname)):
    #print " not found:",oname
    oup.write(i)
    c += 1
 
print "Need to rerun %i" % c
