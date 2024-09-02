#coding:utf-8
'''
Convert the filtered nework.tsv and calculate the 0 and 1 in GT of node in each gene.
'''
import sys

if __name__ == "__main__":
	dir=sys.argv[1]#dir 
	gene=sys.argv[2]
	
	input1=dir+gene+'/'+gene+'network.tsv'
	output0=dir+gene+'/'+gene+'trans_case.tsv'#trans_case.tsv
	output1=dir+gene+'/'+gene+'split_nodeList.tsv'#split_nodeList.tsv
	output2=dir+gene+'/'+gene+'gene_snp.tsv'#gene_snp.tsv
	feature=[]
	snp_info={}
	with open(input1,'r')as f:
		num=0
		for line in f:
			num+=1
			if(num==1):
				continue
			else:
				tmp = line.replace('\n', '').split('\t')
				before=(list(tmp[0:6]))
				other=(list(tmp[6:])) #GT
				count_0=0
				count_1=0
				count_2=0
				for i in range(len(other)):
					if(other[i]=='0|1' or other[i]=='1|0'): #
						count_1+=1
					elif(other[i]=='1|1'):
						count_2+=2
					else:
						count_0+=2
				with open(output0,'a')as f:
					tmp1='\t'.join(before)
					f.write(tmp1+'\t'+ str(count_0+count_1)+'\t'+str(count_1+count_2)+'\n')
				variant=list(tmp[1:5])
				if gene not in snp_info.keys():
					snp_info[gene]=[]
					snp_info[gene].append(variant)	
				else:
					snp_info[gene].append(variant)
												
	for key_gene in snp_info.keys():
		all_variant=key_gene
		for i in snp_info[key_gene]:
			variant='_'.join(i)
			with open (output2,'a')as f:
				f.write(str(key_gene)+'\t'+variant+'\n')
			all_variant+='\t'+(variant)
		with open(output1,'a')as f:
			f.write(str(all_variant)+'\n')

					
