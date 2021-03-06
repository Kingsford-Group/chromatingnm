import glob
import csv

#create dictionary of refFlat data
refdata = {}
with open('refFlat.txt','rb') as reffile:
    reader = csv.reader(reffile, delimiter='\t')
    for row in reader:
        refdata[row[1]] = row
#        print refdata[row[1]]
#        break

with open('generef.txt','wb') as newref:
    refwriter = csv.writer(newref, delimiter='\t')
    with open('leftovergenes.txt','wb') as leftover:
        lowriter = csv.writer(leftover)
        for transfile in glob.glob('geneexpfiles/ENST*'):
            transid = transfile[35:50]
#            print transid
            if transid not in refdata:
                lowriter.writerow([transid])
            else:
                refinfo = [transid, refdata[transid][0], refdata[transid][2], int(refdata[transid][4]), int(refdata[transid][5])]
#                print refinfo
#                break
                refwriter.writerow(refinfo)
