# -*- coding:UTF-8 -*-
from __future__ import print_function
import random
import sys
import os

snp_info={}
snp={}
snp_num=-1
dir=sys.argv[1]
gene=sys.argv[2]
sample_size=2*int(sys.argv[3])-1
per_gene=dir+gene+"/"

def homo(list1):
	num=0
	for i in list1:
		if (i%2==0):
			if (i-1) in list1:
				num+=1
			else:
				continue
		if(i%2==1):
			if(i+1) in list1:
				num+=1
			else:
				continue
	return(num)
def com(list1,list2):
	num=0
	for i in list1:
		if (i%2==0):
			if (i-1) in list2 :
				num+=1
			else:
				continue
		if(i%2==1):
			if(i+1) in list2 :
				num+=1
			else:
				continue
	return(num)
if(os.path.isfile(per_gene+gene+'trans_case.tsv')):
	with open(per_gene+gene+'trans_case.tsv', 'r') as f:
		for line in f:
			snp_num=snp_num+1
			snp_random=[]
			each_snp=line.strip().split('\t')
			randomNum=0
			while(randomNum<int(each_snp[7])): 
				random_t=random.randint(0,sample_size)
				random_t+=1
				if random_t not in snp_random :
					snp_random.append(random_t)
					randomNum+=1
				else:
					randomNum-=1
			variant=each_snp[1]+'_'+each_snp[2]+'_'+each_snp[3]+'_'+each_snp[4]
			
			if each_snp[0] not in snp_info.keys():
				snp={}
				snp[variant]=snp_random
				snp_info[each_snp[0]]=snp
				
			else:
				snp=snp_info[each_snp[0]]
				if variant not in snp.keys():
					snp[variant]=snp_random
					snp_info[each_snp[0]]=snp

	for key_gene, snp in snp_info.items():
		if(key_gene==gene):
			print(key_gene,end='')
			nowNodes = list(snp.keys())
			for i in range(len(nowNodes)-1):
				node1 = nowNodes[i]
				list1 = snp[node1]#
				num1=0
				num1=homo(list1)
				if (num1>0):
					print('',end='\t')
					print(str(nowNodes[i])+':'+str(nowNodes[i])+':'+str(num1),end='')
				for j in range(i + 1, len(nowNodes)):
					num=0
					node2 = nowNodes[j]
					list2 = snp[node2]
					num=com(list1,list2)
					if num>0:
						print('',end='\t')
						print(str(nowNodes[i])+':'+str(nowNodes[j])+':'+str(num),end='')

