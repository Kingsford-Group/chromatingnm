import csv
import itertools
import numpy as np
import glob

genelocdict = {}
genetranscriptdict = {}
with open('fulltranscriptref.txt','rb') as transcriptfile:
    reader = csv.reader(transcriptfile,delimiter = '\t')
    next(reader, None)
    ambiguousgenes = []
    for transcript in reader:
        if len(transcript[2]) <=2:
            if transcript[1] not in genelocdict:
                genelocdict[transcript[1]] = [transcript[2], (int(transcript[3]), int(transcript[4]))]
                genetranscriptdict[transcript[1]] = [transcript[0]]
#                print transcript[1], genedict[transcript[1]]
#                pass
            elif transcript[2] != genelocdict[transcript[1]][0]:
                ambiguousgenes.append(transcript[1])
            else:
                genelocdict[transcript[1]].append((int(transcript[3]), int(transcript[4])))
                genetranscriptdict[transcript[1]].append(transcript[0])
#                print genetranscriptdict[transcript[1]]
    for badgene in ambiguousgenes:
        del genelocdict[badgene]
        del genetranscriptdict[badgene]

with open('fullgeneref.txt','wb') as genefile:
    writer = csv.writer(genefile, delimiter='\t')
    for gene in genelocdict:
#        print gene, genedict[gene]
        transcriptlocs = list(itertools.chain(*genelocdict[gene][1:]))
        genestart = min(transcriptlocs)
        geneend = max(transcriptlocs)
#        print genestart, geneend
#        print [gene, genedict[gene][0], genestart, geneend]
        writer.writerow([gene, genelocdict[gene][0], genestart, geneend])
#        break


for gene in genetranscriptdict:
    expvalues = np.zeros(212)
#    print genetranscriptdict[gene]
    for transcript in genetranscriptdict[gene]:
#        print transcript
        transcriptfile = glob.glob('transcriptexpfiles/'+transcript+'*')[0]
        with open(transcriptfile,'rb') as tfile:
            transexp = tfile.read().splitlines()
            transexp = [float(x) for x in transexp]
        expvalues += transexp
#        print expvalues
    with open('geneexpfiles/'+gene+'.expressionvals.txt', 'wb') as gfile:
        writer = csv.writer(gfile)
        for expval in expvalues:
            writer.writerow([expval])
