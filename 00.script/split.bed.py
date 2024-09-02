#coding:utf-8

import sys
import os
inputfile=sys.argv[1]
dir=sys.argv[2] #$outputDir/split_bed/
gene=sys.argv[3]

with open (inputfile,'r')as f:
	for line in f:
		gene2=line.strip().split('\t')[-1]
		if gene==gene2:
			if os.path.exists(dir+gene+'/'):
				with open(dir+gene+'/'+gene+'.bed','w')as f:
					f.write(line)
			else:
				os.mkdir(dir+gene+'/')
				with open(dir+gene+'/'+gene+'.bed','w')as f:
					f.write(line)
		
