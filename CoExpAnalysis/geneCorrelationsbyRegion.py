import sys
import csv
import glob
import intervaltree as it
import numpy as np

genedict = {}
with open('fullgeneref.txt','rb') as gfile:
    reader = csv.reader(gfile, delimiter = '\t')
    for gene in reader:
        if gene[1] == 'X': gene[1] = '23'
        if gene[1] not in genedict:
            #print gene[1]
            genedict[gene[1]] = [[gene[0], gene[2], gene[3]]]
        else:
            genedict[gene[1]].append([gene[0], gene[2], gene[3]])

regiondict = {}
for chr in genedict:
    if chr != 'MT' and chr != 'Y':
        print 'parsing chromosome '+chr
        genetree = it.IntervalTree()
        if chr == 'X': chr = '23'
        for gene in genedict[chr]:
            genetree[int(gene[1]):int(gene[2])] = gene[0]
        #print genetree
        with open('genecorrs/genecorrelations_chr'+chr+'_usingmidpt.txt','wb') as newfile:
            writer = csv.writer(newfile,delimiter='\t')
            writer.writerow(['Interval1start', 'Interval1end', 'Interval2start','Interval2end','GNMweight','Gene1','Gene2','GeneDistance','GeneCorr'])
            with open(sys.argv[1]+chr+'.txt','rb') as odfile:
                reader = csv.reader(odfile, delimiter = '\t')
                for region in reader:
                    int1genes = sorted(genetree[int(region[0]):int(region[1])])
                    int2genes = sorted(genetree[int(region[2]):int(region[3])])
                    #print region
                    #print int1genes
                    #print int2genes
                    for gene1 in int1genes:
                        gene1midpt = (gene1.begin + gene1.end)/2
                        if int(region[0]) < gene1midpt < int(region[1]):
                            geneexpfile = 'geneexpfiles/'+gene1.data+'.expressionvals.txt'
                            gene1exp = np.genfromtxt(geneexpfile)
                            for gene2 in int2genes:
                                gene2midpt = (gene2.begin + gene2.end)/2
                                if int(region[2]) < gene2midpt < int(region[3]):
                                    geneexpfile = 'geneexpfiles/'+gene2.data+'.expressionvals.txt'
                                    gene2exp = np.genfromtxt(geneexpfile)
                                    corrval = np.corrcoef(gene1exp,gene2exp)[0,1]
                                    regionkey = (int(chr), int(region[0]), int(region[1]), int(region[2]), int(region[3]))
                                    if regionkey not in regiondict and ~np.isnan(corrval):
                                        regiondict[regionkey] = [1,corrval]
                                    elif ~np.isnan(corrval):
                                        regiondict[regionkey][0] += 1
                                        regiondict[regionkey][1] += corrval
                                    #print corrval
                                    writer.writerow([region[0], region[1], region[2], region[3], region[4], gene1.data, gene2.data, abs(gene2midpt - gene1midpt), corrval])

#ODDlist = []
for key in regiondict:
    regiondict[key][1] = regiondict[key][1]/regiondict[key][0]
    #sanity check
    if regiondict[key][1] >= 1:
        print 'problem: ', key
    #ODDlist.append([key[0], key[1], key[2], key[3], key[4], regiondict[key][0], regiondict[key][1]])

# ODDlist = np.asarray(ODDlist)
# sort by correlation values
# ODDlist = ODDlist[ODDlist[:,6].argsort()[::-1]]
# print ODDlist[:5,:]

# with open('ODDregions_sortedbyCoexp.txt','wb') as newfile:
#    writer = csv.writer(newfile, delimiter = '\t')
#    for row in ODDlist:
#        writer.writerow(row)
