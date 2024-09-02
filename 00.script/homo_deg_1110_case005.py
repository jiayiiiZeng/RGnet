# coding:utf-8

###
#The path to main_random_1110_case.py in line 34 must be updated according to the actual situation.
###

import numpy as np
import random
import json
import networkx as nx
import os
import sys
import subprocess

def deg(nodeList,nodeNum,gene,gene_tmp,random_num):
	count=nodeNum
	file_name=str(random_num)
	G = nx.Graph()
	for pair in gene_tmp:
		if(pair!=''):
			each_pair=pair.split(':')
			G.add_edge(str(each_pair[0]),str(each_pair[1]), weight=float(each_pair[2]))

	with open(outputdir+gene,'a') as f :
		for m in range(count):
			if nodeList[m] not in G.nodes():
				f.write(gene+'\t'+file_name+'\t'+nodeList[m]+'\t'+str('0')+'\t'+str('0')+'\n')
			else:
				f.write(gene+'\t'+file_name+'\t'+nodeList[m]+'\t'+str(G.degree(nodeList[m]))+'\t'+str(G.degree(nodeList[m],'weight'))+'\n')
	
	return G

def mult(random_num):
	out=subprocess.check_output(["python",'/public/home/zengjiayi/work/RGnet/00.script/main_random_1110_case.py',input_per_gene_dir,nowgene,samplesize])
	out = out.decode()
	tmp=out.split('\t') 
	gene_tmp=tmp[0]
	nodeList=node_tmp[gene_tmp]
	nodeNum=len(nodeList)
	G_random=deg(nodeList,nodeNum,gene_tmp,tmp[1:],random_num)
	edgesNum=G_random.number_of_edges()
	random_edges.append(edgesNum)
	random_source_deg_num=0
	for node in G_random.nodes():
		if(G_random.degree(node)>1):
			random_source_deg_num+=1
	random_node_degree.append(random_source_deg_num)
	
	return G_random,gene_tmp[0]

if __name__ == "__main__":
	outputdir=sys.argv[1]
	input_per_gene_dir=sys.argv[2]
	nowgene=sys.argv[5]
	samplesize=sys.argv[4]
	split_nodeList=input_per_gene_dir+nowgene+'/'+nowgene+'split_nodeList.tsv'
	print(nowgene)
	node_tmp={}
	with open(split_nodeList,'r')as f:
		for line in f:
			tmp=line.strip().split('\t')
			node_tmp[tmp[0]]=tmp[1:]
	random_edges=[]
	random_node_degree=[]
	for random_num in range(int(sys.argv[3])):
		G_random,gene=mult(random_num)
		edgesNum=G_random.number_of_edges()
		random_source_deg_num=0
		for node in G_random.nodes():
			if(G_random.degree(node)>1):
				random_source_deg_num+=1
		with open(outputdir+nowgene+'random.tsv','a')as f:
			f.write(str(random_source_deg_num)+'\t'+str(edgesNum)+'\n')
		
	

