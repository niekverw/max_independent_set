#!/usr/bin/env python3


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
		s_mis=(1-df_sampleqc[ df_sampleqc['identifier'].isin(maxindepset)   ]['sample.qc.missing.rate']).sum()
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

file_subset=sys.argv[1]  
file_relatedness=sys.argv[2]  
file_sampleqc=sys.argv[3]  
file_output=file_subset + ".rel.tsv"

df_neid_subset=pd.read_csv(file_subset, sep=' ')
df_relatedness=pd.read_csv(file_relatedness, sep=' ')#[:100]
df_sampleqc=pd.read_csv(file_sampleqc, sep='\t')
#print(df_relatedness.dtypes)
#print(df_sampleqc.dtypes)
print(df_neid_subset.dtypes)
############ ALL INDIVIDUALS ########################
#######
condition=((df_sampleqc['het.missing.outliers']==1) | (df_sampleqc['Submitted.Gender']!=df_sampleqc['Inferred.Gender']))
identifier_outliers=df_sampleqc.loc[ condition ]['identifier'].values
print(df_relatedness.shape)
df_relatedness=df_relatedness[  -(df_relatedness['ID1'].isin(identifier_outliers) | df_relatedness['ID2'].isin(identifier_outliers)) & (df_relatedness['ID1'].isin(df_neid_subset.ix[:,0]) | df_relatedness['ID2'].isin(df_neid_subset.ix[:,0] ) ) ]
print(df_relatedness.shape)
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
df_dct_clust['identifier']=df_dct_clust.index
df_dct_clust.columns =['Genrelclust', 'identifier']
df_sampleqc=df_sampleqc.join(df_dct_clust.set_index(['identifier']), on=['identifier'])
df_sampleqc.loc[ -df_sampleqc['identifier'].isin(G.nodes()) & -condition, 'Genrelclust']=df_sampleqc['identifier']

### COMBINE INDEPENDENT SAMPLES WITH SAMPLEQC DF.
df_sampleqc['GenIndependent']= np.nan
df_sampleqc.loc[ df_sampleqc['identifier'].isin(G.nodes()), 'GenIndependent']=0

df_sampleqc.loc[ df_sampleqc['identifier'].isin(vct_independentsamples), 'GenIndependent']=1
df_sampleqc.loc[ -df_sampleqc['identifier'].isin(G.nodes()) & -condition, 'GenIndependent']=1 # everything not in G.nodes + condition is also independent


############################################################
df_sampleqc=df_sampleqc[['identifier','Genrelclust','GenIndependent']]
print("writing..",file_output+".rel.tsv")
df_sampleqc.to_csv(file_output+".rel.tsv",na_rep='NA', sep='\t',index=False)


exit()