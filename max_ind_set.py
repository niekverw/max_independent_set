#!/usr/bin/env python3
#

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import random
import sys
import numpy as np
random.seed(13379)


def pick_best_candidate(H, nodes, keyfunc=None):
	"""unless a key funtion is supploed it makes a random choice from the set it is handed."""
	if keyfunc == None:
		selected_node = random.choice( nodes )
	else:
		selected_node = min(nodes, key=keyfunc)
	return( selected_node )

def findindependentsubset(H):
#### find maximum:

	cnt_maxindep=0
	tmp_maxindep=0
	tmp_max_n_mis=0
	cnt_iter=0
	independentset=[]

	if len(H) < 1:
		pass
	elif len(H) == 1:
		independentset.append( H[0] )
		return(independentset)
	# elif len(H) == 2:
	# 	picked_node = pick_best_candidate(H, H.nodes())

	# 	independentset.append( picked_node)
	# 	return(independentset)

	while True:
		cnt_iter+=1
		maxindepset=nx.maximal_independent_set(H)
		## get sum of missingness
		s_mis=(1-df_sampleqc[ df_sampleqc['n_eid'].isin(maxindepset)   ]['sample.qc.missing.rate']).sum()
		maxindep_n=len(maxindepset)
		max_n_mis=s_mis
		sys.stdout.write("#nodes in subgraph: {} , current independent: {},current n_mis {} , iterations {},max indep {} : {}, max 1-missingness {} ,cnt_maxindep: {}, \r".format(len(H.nodes()),maxindep_n,round(max_n_mis),cnt_iter,round(tmp_maxindep),len(independentset),round(tmp_max_n_mis,3),cnt_maxindep) )


		if max_n_mis>tmp_max_n_mis:
			cnt_maxindep=0 #reset
		tmp_maxindep=max(maxindep_n,tmp_maxindep)
		tmp_max_n_mis=max(max_n_mis,tmp_max_n_mis)
		if tmp_max_n_mis==max_n_mis:
			independentset=maxindepset
			cnt_maxindep+=1
			cnt_iter-=(len(H.nodes())/2)
			if cnt_maxindep==10:
				print("")
				return(independentset)

		if(cnt_iter>=len(H.nodes())*2 or cnt_iter>=3000) :
			print("")
			return(independentset)

file_neid_subset="/path/to/file/with/identifiers.txt"
suffix="subset"


file_relatedness="ukb962_rel_s488374.dat"# .headtail"
file_sampleqc="/ukb_sqc_v2.txt.gz.processed.tsv" #.headtail"
file_output=file_neid_subset + ".rel.tsv"

df_neid_subset=pd.read_csv(file_neid_subset, sep=' ')
df_relatedness=pd.read_csv(file_relatedness, sep=' ')#[:100]
df_sampleqc=pd.read_csv(file_sampleqc, sep='\t')
#print(df_relatedness.dtypes)
#print(df_sampleqc.dtypes)
print(df_neid_subset.dtypes)
############ ALL INDIVIDUALS ########################
#######
### FILTER:
condition=((df_sampleqc['het.missing.outliers']==1) | (df_sampleqc['Submitted.Gender']!=df_sampleqc['Inferred.Gender']))
n_eid_outliers=df_sampleqc.loc[ condition ]['n_eid'].values
print(df_relatedness.shape)
df_relatedness=df_relatedness[  -(df_relatedness['ID1'].isin(n_eid_outliers) | df_relatedness['ID2'].isin(n_eid_outliers)) & (df_relatedness['ID1'].isin(df_neid_subset.ix[:,0]) & df_relatedness['ID2'].isin(df_neid_subset.ix[:,0] ) ) ]
print(df_relatedness.shape)
print(df_sampleqc.shape)
df_sampleqc=df_sampleqc[ -(df_sampleqc['n_eid'].isin(n_eid_outliers)) & df_sampleqc['n_eid'].isin(df_neid_subset.ix[:,0])  ]
print(df_sampleqc.shape)

####
G=nx.from_pandas_dataframe(df_relatedness, 'ID1', 'ID2' )  # convert to graph

vct_independentsamples=[]
dct_clust={}
c=0
for s in sorted(nx.connected_components(G), key=len, reverse=True):
	c+=1 #cluster.
	H=G.subgraph(s)
	vct_independentsamples+=findindependentsubset(H) ## find independent subjects in subgraph.
	for key in H.nodes():
		dct_clust[key]=c

### COMBINE CLUSTER WITH SAMPLEQC DF.
df_dct_clust=pd.DataFrame([dct_clust]).transpose()
df_dct_clust['n_eid']=df_dct_clust.index

StrGenrelclust='Genrelclust'+suffix
StrGenIndependent='GenIndependent'+suffix

df_dct_clust.columns =[StrGenrelclust, 'n_eid']
df_sampleqc=df_sampleqc.join(df_dct_clust.set_index(['n_eid']), on=['n_eid'])
df_sampleqc.loc[ -df_sampleqc['n_eid'].isin(G.nodes()) & -condition, StrGenrelclust]=df_sampleqc['n_eid']

### COMBINE INDEPENDENT SAMPLES WITH SAMPLEQC DF.
df_sampleqc[StrGenIndependent]=1
df_sampleqc.loc[ df_sampleqc['n_eid'].isin(G.nodes()), StrGenIndependent]=0 #replace all that are a node with 0
df_sampleqc.loc[ df_sampleqc['n_eid'].isin(vct_independentsamples), StrGenIndependent ]=1 # replace all that are independent
#df_sampleqc.loc[ -df_sampleqc['n_eid'].isin(G.nodes()) & -condition, 'GenIndependent']=1 # everything not in G.nodes + condition is also independent


############################################################
df_sampleqc=df_sampleqc[['n_eid',StrGenrelclust,StrGenIndependent]]
print("writing..",file_output+".rel.tsv")
df_sampleqc.to_csv(file_output+".rel.tsv",na_rep='NA', sep='\t',index=False)


exit()