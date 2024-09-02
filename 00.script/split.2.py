# coding:utf-8

import sys
import os

if __name__ == "__main__":
	input1=sys.argv[1] #dir pempirical
	out=sys.argv[2]  #dir homo_split_gene
	gene=sys.argv[3]
	if(os.path.isfile(input1+gene )):
		with open(input1+gene,'r')as f:
			for line in f:
				tmp=line.strip().split('\t')
				nowlabel=tmp[0]+'__'+tmp[2] #gene__chr_pos_ref_alt
				with open(out+nowlabel,'a')as f:
					f.write(line)
				
