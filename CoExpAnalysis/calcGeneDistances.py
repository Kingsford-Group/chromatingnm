import csv
import numpy as np

genedict = {}
with open('fullgeneref.txt','rb') as gfile:
    reader = csv.reader(gfile, delimiter='\t')
    for gene in reader:
        if len(gene[1]) <= 2 and gene[0] not in genedict:
            genedict[gene[0]] = [gene[1], int(gene[2]), int(gene[3])]
        elif len(gene[1]) <= 2 and gene[0] in genedict:
            print gene[0]+' duplicate entry'

genepairdict = {}
with open('allgenepairdistances.txt','wb') as distfile:
    writer = csv.writer(distfile, delimiter='\t')
    writer.writerow(['Gene1', 'Gene2', 'Chr', 'MidptDistance', 'Coexpression'])
    for gene1 in genedict:
        for gene2 in genedict:
            if (gene1, gene2) not in genepairdict and (gene2, gene1) not in genepairdict and genedict[gene1][0] == genedict[gene2][0]:
                gene1expfile = '/mnt/disk17/RNAtoTADs/geneexpfiles/'+gene1+'.expressionvals.txt'
                gene1exp = np.genfromtxt(gene1expfile)
                gene2expfile = '/mnt/disk17/RNAtoTADs/geneexpfiles/'+gene2+'.expressionvals.txt'
                gene2exp = np.genfromtxt(gene2expfile)
                gene1midpt = (genedict[gene1][1] + genedict[gene1][2])/2
                gene2midpt = (genedict[gene2][1] + genedict[gene2][2])/2
                genepairdict[(gene1,gene2)] = [genedict[gene1][0], abs(gene1midpt - gene2midpt)]
                writer.writerow([gene1, gene2, genedict[gene1][0], abs(gene1midpt - gene2midpt), np.corrcoef(gene1exp,gene2exp)[0,1]])
