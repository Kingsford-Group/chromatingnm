import os, glob
import re, sys

# salmon-v0.6.0/bin/
# salmon quant -i transcripts_index -r reads.fa -o transcripts_quant
# script-run-salmon-paired defaults to the following:
# salmon quant -i transcripts_index -l IU/ISR -1 reads1.fa -2 reads2.fa -o transcripts_quant

read_dir = sys.argv[1]
outhead =sys.argv[2]

unique_reads =  set([f.split(".")[0].split("_")[0] for f in glob.glob(os.path.dirname(read_dir) +"/*.fasta.gz")])
#print unique_reads

for read in unique_reads:
    outname = outhead +os.path.basename(read)
    #print outname
    files = sorted(glob.glob("{}*".format(read)))
    #print files
    if len(files) == 1:
        #pass
        os.system("./script-run-salmon-single.sh {} {}".format(files[0], outname))
    elif len(files) == 2:
        #pass
        os.system("./script-run-salmon-paired.sh {} {} {}".format(files[0], files[1], outname))
    else:
        print "This should not happen"
        exit(1)
