import csv
from scipy import stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plotly.plotly as py
import glob
#import matplotlib.patches as mpatches
#import sys

oddgenedict = {}
oddgenevals = []
for oddfilename in glob.glob('/mnt/disk17/RNAtoTADs/genecorrs/*usingmidpt.txt'):
    chrnum = oddfilename.split('_')[1]
    if chrnum == 'chrX': chrnum = 'chr23'
    #print chrnum
    #print oddfilename
    with open(oddfilename,'rb') as oddfile:
        reader = csv.reader(oddfile, delimiter='\t')
        next(reader, None)
        for oddgp in reader:
            #print oddgp
            #break
            oddgenedict[(oddgp[5], oddgp[6])] = [int(oddgp[7]),float(oddgp[8]),int(chrnum[3:])]
            oddgenevals.append([int(oddgp[7]), float(oddgp[8]), int(chrnum[3:])])
oddgenevals = np.asarray(oddgenevals)
#print len(oddgenedict)
#print 'max odd gene dist: ', max(oddgenevals[:,0])
#print 'min odd gene dist: ', min(oddgenevals[:,0])
mingenedist = min(oddgenevals[:,0])

randgenedict = {}
randgenevals = []
with open('/mnt/disk17/RNAtoTADs/allgenepairdistances.txt','rb') as allgpfile:
    reader = csv.reader(allgpfile, delimiter='\t')
    next(reader, None)
    for genepair in reader:
        if (genepair[0], genepair[1]) not in oddgenedict and (genepair[1], genepair[0]) not in oddgenedict and int(genepair[3]) >= mingenedist:
            randgenedict[(genepair[0], genepair[1])] = [int(genepair[3]), float(genepair[4])]
            randgenevals.append([int(genepair[3]), float(genepair[4])])
randgenevals = np.asarray(randgenevals)
#print 'max rgp distance:', max(randgenevals[:,0])

x = 0
delta = 25000000
maxn = 5
#maxn = int(np.ceil(max(oddgenevals[:,0])/delta))
#print maxn

#for n in xrange(1,maxn):
    #print ' '
    #print 'Distribution for n = '+str(n)
    #lowerdistbound = x+(n-1)*delta
    #upperdistbound = x+n*delta
    #plot histogram of chr numbers used in this dist range
    #ax0 = plt.figure()
    #chrnumsinrange = oddgenevals[(lowerdistbound <= oddgenevals[:,0]) & (oddgenevals[:,0] < upperdistbound),2]
    #for chrnum in xrange(1,24):
        #chrnums = chrnumsinrange.tolist()
        #print 'gp from chr '+str(chrnum)+': '+str(chrnums.count(chrnum))
        #bins = np.linspace(9,22,14)
        #ax0.set_title('dist of chrs from 80M-90M range')
        #plt.hist(chrnumsinrange,bins)
        #plt.savefig('chrdistforgenecorrs.pdf', format='pdf', dpi=1000)


plt.figure(figsize=(17.8/2.54, 9.25/2.54))
for n in xrange(1,maxn):
    lowerdistbound = x+(n-1)*delta
    upperdistbound = x+n*delta
    
    ax1 = plt.subplot(2,2,n)
    if lowerdistbound == 0:
        plottitle = '0-'+str(upperdistbound)[:-6]+'M bp'
    else:
        plottitle = str(lowerdistbound)[:-6]+'M-'+str(upperdistbound)[:-6]+'M bp'
    ax1.set_title(plottitle)

    bins = np.linspace(-0.5,1,100)
    rgvalsinrange = randgenevals[(lowerdistbound <= randgenevals[:,0]) & (randgenevals[:,0] < upperdistbound),1] 
    hist1 = ax1.hist(rgvalsinrange, bins, alpha=0.5,label='background gene pairs', normed=True, edgecolor='none')
    plt.xlim((-0.5,1))
    plt.locator_params(axis='y', nbins=5)
    plt.locator_params(axis='x', nbins=6)
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
        

    oddvalsinrange = oddgenevals[(lowerdistbound <= oddgenevals[:,0]) & (oddgenevals[:,0] < upperdistbound),1]
    ax2 = ax1.twinx()
    ax2.hist(oddvalsinrange, bins, color='y', alpha=0.5, label='ODD gene pairs', normed=True,edgecolor='none')
    plt.xlim((-0.5,1))
    plt.locator_params(axis='y', nbins=5)
    plt.locator_params(axis='x', nbins=6)
    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    if n < 3:
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
    
    print ''
    print 'mean of CCDD vals: ', np.nanmean(oddvalsinrange)
    print 'mean of random gp vals: ', np.nanmean(rgvalsinrange)
    print 'number of CCDD vals in range: ', oddvalsinrange.shape[0]
    print 'KS test results:', stats.ks_2samp(oddvalsinrange, rgvalsinrange)

#plt.setp([a.get_xticklabels() for a in ax1[0,:]], visible=False)
plt.savefig('ccdd_coexpression_histogram_figure.pdf', format='pdf', dpi=1000)
    #plt.clf()
