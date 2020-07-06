import os, sys
import pysam
import argparse
import datetime
parser = argparse.ArgumentParser(description='Read fragment information from bam file.')
parser.add_argument('-in','--bamfile', type=str, help='the bam file to read')
parser.add_argument('-out','--outFolder', type=str, help='The folder to output the fragment information')

args = parser.parse_args()

path_in=args.bamfile
os.mkdir(args.outFolder);
samfile = pysam.AlignmentFile(path_in, "rb")
start = datetime.datetime.now()
for read in samfile.fetch():
    if ('H' in read.reference_name) or ('h' in read.reference_name):
        chrm = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
        'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']
    else:
        chrm = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
        '20', '21', '22']
    break

for ch in chrm:
    if ('H' in ch) or ('h' in ch):
        name = ch
    else:
        name = 'chr'+ch
    file = open(args.outFolder+'/'+name+".txt", "w")
    for read in samfile.fetch(ch):
        te=read.cigartuples[0]
        file.write(str(read.qname) + "    " + str(read.pos+1) + "    " + str(read.isize) + "    " + str(read.mpos+1)+"\n")
    file.close()
end = datetime.datetime.now()
print (end-start)
