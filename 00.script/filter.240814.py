import sys
import re
dir=sys.argv[1]
inputTranscript=sys.argv[2]
gene=sys.argv[3]
input=dir+gene+'/'+gene+'.tsv'
output=dir+gene+'/'+gene+'network.tsv'
outRawTemp=dir+'temp/'+gene+'.raw.temp.tsv'
outDelTemp=dir+'temp/'+gene+'.del.temp.tsv'

## match unique transcripts
nowTranscript=''
with open(inputTranscript,'r')as f:
	for line in f:
		tmp = line.strip().split('\t')
		if tmp[0]==gene:
			nowTranscript=tmp[1]
			print('nowTranscript',gene,nowTranscript)
		else:
			continue

## Filter original tsv files based on unique transcripts
with open(input,'r')as f,open(outRawTemp,"w") as f1,open(outDelTemp,"w") as f2:
	maxP = 0.005
	lineNum = 0
	for line in f:
		lineNum += 1
		temp=line.strip().split('\t')
		if lineNum == 1:
			header = [i.split("]")[1].split(":")[0] for i in temp]
			f1.write("\t".join(header)+'\n')
		else:
			VEP_SYMBOL = temp[0].split("|") ##split by |
			CHROM = temp[1]
			POS = temp[2]
			REF = temp[3]
			ALT = temp[4]
			snp = "_".join(temp[1:5])
			VEP_Feature = temp[5].split("|") ##split by |
			VEP_HGVSc = temp[6].split("|") ##split by |
			VEP_HGVSp = temp[7].split("|") ##split by |
			VEP_IMPACT = temp[8].split("|") ##split by |
			VEP_Consequence = temp[9].split("|") ##split by |
			AF = temp[10:17]
			sample = temp[17:]
			for trans in VEP_Feature:
				if trans == nowTranscript:
					index = VEP_Feature.index(trans)
					f1.write("\t".join([VEP_SYMBOL[index],CHROM,POS,REF,ALT,VEP_Feature[index],VEP_HGVSc[index],VEP_HGVSp[index],VEP_IMPACT[index],VEP_Consequence[index],"\t".join(AF),"\t".join(temp[17:])])+'\n')
				elif nowTranscript not in VEP_Feature:
					print("the snp was not annotated by the transcript")
					f2.write(line)
## Label sites
with open(outRawTemp,'r')as f1,open(output, 'w') as f2:
	linenum = 0
	sampleIDs = []
	for line in f1:
		linenum += 1
		temp=line.strip().split('\t')
		genotypes=temp[17:]
		snp = "_".join(temp[1:5])
		VEP_IMPACT = temp[8]
		VEP_Consequence = temp[9]
		if linenum == 1:
			sampleIDs = temp[17:]
			f2.write("\t".join(["gene","CHROM","POS","REF","ALT","VEP_IMPACT","\t".join(sampleIDs)])+"\n")
		elif re.search('0|0',str(genotypes)) or re.search('0|1',str(genotypes)) or re.search('1|0',str(genotypes)) or re.search('1|1',str(genotypes)): ##Remove the location where the genotype is all .|.
			f2.write('\t'.join(["\t".join(temp[0:5]),temp[8],"\t".join(genotypes)])+'\n')
		else:
			continue
