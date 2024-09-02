# coding:utf-8
import numpy as np
import random
import json
import networkx as nx
import matplotlib
matplotlib.use('cairo')  # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from matplotlib.artist import Artist
import seaborn as sns
import community
from igraph import *
import os
import re
import sys
sns.set_style('whitegrid')
import pandas as pd

class GraphArtist(Artist):
	def __init__(self, graph,bbox, palette=None, *args, **kwds):
		Artist.__init__(self)
		if not isinstance(graph, Graph):
			raise TypeError("expected igraph.Graph, got %r" % type(graph))
		self.graph = graph
		self.palette = palette or palettes["gray"]
		self.bbox = BoundingBox(bbox)
		self.args = args
		self.kwds = kwds
	def draw(self, renderer):
		from matplotlib.backends.backend_cairo import RendererCairo
		if not isinstance(renderer, RendererCairo):
			raise TypeError("graph plotting is supported only on Cairo backends")
		self.graph.__plot__(renderer.gc.ctx, self.bbox, self.palette, *self.args, **self.kwds)
	
def readTable(gene):
	nodeList = []
	totalPos = []
	with open(inputfile, 'r') as f:
		lineNum = 0
		for line in f:
			lineNum += 1
			tmp = line.replace('\n', '').split('\t')
			if lineNum == 1:
				sampleIDs = list(tmp[6:]) #GT
				totalList = []
				for i in range(len(sampleIDs)):
					totalList.append([])
			if tmp[0] != gene:
				continue
				
			sampleTmp = list(tmp[6:])
			node=tmp[1] + '_' + tmp[2] + '_' + tmp[3] + '_' + tmp[4]
			label = tmp[5].lower()
			if  node  in snp_p:
				nowLabel = 'P'
			elif  node in snp_benign:
				nowLabel = 'B'
			else:
				nowLabel = 'V'
			nodeList.append(node+'_'+nowLabel)

			for i in range(len(sampleTmp)):
				totalList[i].append(sampleTmp[i])

	return nodeList, sampleIDs, totalList

def buildG0(nodeList, sampleIDs, totalList,gene):
	nodeNum = len(nodeList)
	matrix = np.zeros((nodeNum, nodeNum))
	sampleDict = {}

	for i in range(len(totalList)):
		sampleID = sampleIDs[i]
		colList = totalList[i]	
		for j in range(nodeNum):
			nodeStr = nodeList[j]
			if colList[j] == '0|1' or colList[j] == '1|0' : 
				if sampleID not in sampleDict.keys():
					nodeNumDict = {}
					nodeNumDict[nodeStr] = 1
					sampleDict[sampleID] = nodeNumDict
				else:
					nodeNumDict = sampleDict[sampleID]
					if nodeStr not in nodeNumDict.keys():
						nodeNumDict[nodeStr] = 1
					else:
						nodeNumDict[nodeStr] += 1
					sampleDict[sampleID] = nodeNumDict
			elif colList[j]=='1|1':
				if sampleID not in sampleDict.keys():
					nodeNumDict = {}
					nodeNumDict[nodeStr] = 2
					sampleDict[sampleID] = nodeNumDict
				else:
					nodeNumDict = sampleDict[sampleID]
					if nodeStr not in nodeNumDict.keys():
						nodeNumDict[nodeStr] = 2
					else:
						nodeNumDict[nodeStr] += 2
					sampleDict[sampleID] = nodeNumDict
			else:
				continue
	node_in_graph = []
	for key, nodeNumDict in sampleDict.items():
		nowNodes = list(nodeNumDict.keys())
		for i in range(len(nowNodes)):
			node1 = nowNodes[i]
			if not node1 in node_in_graph:
				node_in_graph.append(node1)
			index1 = nodeList.index(node1)
			n1 = nodeNumDict[node1]
			if(n1==2):
				matrix[index1, index1]+=0.5
			for j in range(i+1, len(nowNodes)):
				node2 = nowNodes[j]
				if not node2 in node_in_graph:
					node_in_graph.append(node2)
				index2 = nodeList.index(node2)
				n2 = nodeNumDict[node2]
				matrix[index1, index2] += 1
				matrix[index2, index1] =matrix[index1, index2]
	G = nx.Graph()
	for i in range(nodeNum):
		for j in range(i,nodeNum):
			if matrix[i][j] != 0:
				G.add_edge(nodeList[i], nodeList[j],weight=matrix[i][j])
	source_deg_num=0
	for node in G.nodes():
		if(G.degree(node)>1):
			source_deg_num+=1
		else:
			continue

	source_RG_types=G.number_of_edges()
	edges = []
	#####################################################
	##igraph
	#####################################################
	g0 = Graph()
	vertices_num = len(node_in_graph)
	node_matrix = np.zeros((vertices_num, vertices_num))
	for i in range(vertices_num):
		for j in range(i,vertices_num):
			node1 = node_in_graph[i]
			node2 = node_in_graph[j]
			node1_index = nodeList.index(node1)
			node2_index = nodeList.index(node2)
			if matrix[node1_index][node2_index] > 0:
				edges.append((i,j))
				
	g0.add_vertices(vertices_num)
	g0.add_edges(edges)
	colorDict = {'P': '#d92722', 'B': '#2e8c7a', 'V': '#c0c0c0'}##dcdcdc,#c0c0c0
	colors = []
	for i in range(vertices_num):
		colors.append(colorDict[node_in_graph[i].split('_')[-1]])
	
	g0.vs["color"] = colors
	layout = g0.layout("fr")
	visual_style = {}
	visual_style["edge_width"]=0.2
	visual_style["edge_color"]="black"
	visual_style["layout"] = layout
	visual_style["edge_curved"] = 0.3
	visual_style["vertex_size"] = 2 
	visual_style["vertex_frame_width"] = 0
	fig=plt.figure(figsize=(1.5748,1.5748))
	ax1 =fig.add_subplot(111)
	graph_artist=GraphArtist(g0,bbox=(1,1,110,110),**visual_style)
	graph_artist.set_zorder(float('inf'))
	ax1.artists.append(graph_artist)
	plt.axis("off")
	fig.savefig(graphdir+gene+'source'+'.pdf')

	if (len(edges)>0):
		largest_cc_nodes = max(nx.connected_components(G), key=len)
		Lcc = G.subgraph(largest_cc_nodes).copy()
		node2NumDict = {}
		lcc_edges = []
		for edge in Lcc.edges():
			n0 = edge[0]
			n1 = edge[1]
			if n0 not in node2NumDict.keys():
				node2NumDict[n0] = len(node2NumDict)
			if n1 not in node2NumDict.keys():
				node2NumDict[n1] = len(node2NumDict)
			lcc_edges.append((node2NumDict[n0],node2NumDict[n1]))
		g1 = Graph()
		g1.add_vertices(len(node2NumDict))
		g1.add_edges(lcc_edges)
		
		colorDict = {'P': '#d92722', 'B': '#2e8c7a', 'V': '#c0c0c0'}
		color = []
		for node in node2NumDict.keys():
			color.append(colorDict[node.split('_')[-1]])

		g1.vs["color"] = color
		layout = g1.layout("fr")
		visual_style = {}
		visual_style["edge_width"]=0.2
		visual_style["edge_color"]="black"
		visual_style["layout"] = layout
		visual_style["edge_curved"] = 0.3
		visual_style["vertex_size"]= 2
		visual_style["vertex_frame_width"] = 0
		fig=plt.figure(figsize=(1.5748,1.5748))
		ax1 =fig.add_subplot(111)
		graph_artist=GraphArtist(g1,bbox=(1,1,110,110),**visual_style)
		graph_artist.set_zorder(float('inf'))
		ax1.artists.append(graph_artist)
		plt.axis("off")
		fig.savefig(graphdir+gene+'lcc'+'.pdf')
	
	return G,source_deg_num,source_RG_types,matrix

def run_proc1(nodeList, sampleIDs,totalList,gene):
	coffList=[]
	G,source_deg_num,source_edges,matrix= buildG0(nodeList, sampleIDs,totalList,gene)
	coffList=[]
	coffList.append('snp\tdegs\twdegs\n')	
	with open(outputdir + gene+'_cas.tsv', 'w') as f:
		f.write('\n'.join(coffList))
	edges = list(G.edges())
	for i in range(len(nodeList)):
		n=0
		for edge in edges:
			if nodeList[i] in edge:
				n+=1
		with open(outputdir + gene+'_cas.tsv', 'a') as f:
			if nodeList[i] not in G.nodes():
				f.write(nodeList[i]+'\t'+str('0')+'\t'+str('0')+'\n')
			else:
				f.write(nodeList[i]+'\t'+str(n)+'\t'+str(G.degree(nodeList[i],'weight'))+'\n')
	return source_deg_num,source_edges

def pempirical(randomlist,source):
	count=0
	for i in range(len(randomlist)):
		if(int(randomlist[i])>source or int(randomlist[i])==source):
			count+=1
	p=count/randomNum
	return(p)			

def p_sites(input_p):
	with open(input_p,'r')as f:
		linenum=0
		snp_p=[]
		for line in f:
			linenum+=1
			if(linenum==1):
				continue
			else:
				head=line.strip().split('\t')[1:5]
				snp_p.append('_'.join(head))
	return snp_p

def benign_sites(benign):
	with open(benign,'r')as f:
		linenum=0
		snp_benign=[]
		for line in f:
			linenum+=1
			if(linenum==1):
				continue
			else:
				head=line.strip().split('\t')[1:5]
				snp_benign.append('_'.join(head))
	return snp_benign

if __name__ == "__main__":
	randomNum=int(sys.argv[1])
	input1=sys.argv[2] #outputDir1/pempirical/
	graphdir=sys.argv[3] #graph
	outputdir=sys.argv[4] #source deg wDegs
	inputdir=sys.argv[5] #source file 
	input_p=sys.argv[6]
	benign=sys.argv[7]
	gene=sys.argv[8]
	inputfile=inputdir+gene+'/'+gene+'network.tsv'
	snp_p=p_sites(input_p)
	snp_benign=benign_sites(benign)
	nodeList, sampleIDs, totalList = readTable(gene)
	source_deg_num,source_edges=run_proc1(nodeList, sampleIDs,totalList,gene)
	random_edges=[]
	random_node_degree=[]
	with open(input1+gene+'random.tsv','r')as f:
		for line in f:
			random=line.strip().split('\t')
			random_node_degree.append(int(random[0])) #count of random node degree>1
			random_edges.append(int(random[1])) #count of random edges
	p_edges=pempirical(random_edges,source_edges)
	p_deg=pempirical(random_node_degree,source_deg_num)

	## draw graph
	plt.figure(1,figsize=(6,2),dpi=300)
	plt.subplot(121)
	plt.hist(random_node_degree)
	ymin,ymax = plt.ylim()
	plt.plot([source_deg_num,source_deg_num],[ymin,ymax])
	name=gene+'_node_P:'+str(p_deg)
	plt.title(name)
	plt.subplot(122)
	plt.figure(1,figsize=(6,3),dpi=300)
	plt.hist(random_edges)
	ymin,ymax = plt.ylim()
	plt.plot([source_edges,source_edges],[ymin,ymax])
	name=gene+'_edges_P:'+str(p_edges)
	plt.title(name)
	plt.show()	


