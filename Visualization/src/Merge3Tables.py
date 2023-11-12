#%%
# Merge the three - Abundance, Browsing, and Peptide - tables into one
#
import pandas as pd
import os
import json
# import re
# import warnings
# import plotly.graph_objects as go
# To export static images from go, need pip install -U kaleido
# import plotly.express as px
import numpy as np
# import ExpoDecayFit as edf
# from scipy.stats import hmean # harmonic mean
# import matplotlib.pyplot as plt
# import seaborn as sns
#

#%%
filepath1 = os.path.join("..","data","ProteinTurnover","dfProteins20230808.csv") # Had to rename 'Gene' to 'gene' for consistency
tPeptides = pd.read_csv(filepath1, index_col="Gene")
tPeptides.index.name = 'gene'

filepath2 = os.path.join("..","data","ProteinT12Browsing","ProteinT12Browsing.csv")
tBrowsing = pd.read_csv(filepath2, index_col="gene")

filepath3 = os.path.join("..","data","ProteinAbundance","ProteinAbundance.csv")
tAbundance = pd.read_csv(filepath3, index_col="gene")

# %%
tPeptides.info()
tBrowsing.info()
tAbundance.info()
print(tPeptides.columns)
print(tBrowsing.columns)
print(tAbundance.columns)

#%%[markdown]
# 'gene' is the index \
# ['prtn', 'gene', 'desc','proteinT12', 'chart', 'Npeptides', 'peptides'] .
# tPeptides: ['gene', 'Protein', 'Protein_Description', 'b_CFit_all', 't12_CFit_all', 'b_CFit_pass', 't12_CFit_pass', 'chart', 'peptides', 't12_CFit_best', 'rank'] 
# Need to drop [ 'b_CFit_all',' b_CFit_pass', 'rank'] 
# and left with tPeptides: ['gene', 'Protein', 'Protein_Description', 't12_CFit_all', 't12_CFit_pass', 'chart', 'peptides', 't12_CFit_best'] 
# tBrowsing: ['gene', 'rank', 'UniProtID', 'pDesc', 't12', 'mneuron', 'cneuron', 'subLoc']
# tAbundance: ['gene', 'UniProtID', 'pDesc', 'pLength', 'pMass', 'subLoc', 'abundanceLevel', 'avgCopyN'] 
# 
# Will rename as 
# tPeptides: ['gene', 'UniProtID1', 'pDesc1', 't12_CFit_all', 't12_CFit_pass', 'chart', 'peptides', 't12_1', 'rank']
# tBrowsing: ['gene', 'rank', 'UniProtID', 'pDesc', 't12', 'mneuron', 'cneuron', 'subLoc']
# tAbundance: ['gene', 'UniProtID', 'pDesc', 'pLength', 'pMass', 'subLoc', 'abundanceLevel', 'avgCopyN'] 
# .

#%%
tPeptides.drop([ 'b_CFit_all', 'b_CFit_pass', 'rank' ] , axis=1, inplace=True)
#%%
tPeptides.rename(columns={'Protein':'UniProtID-001', 'Protein_Description':'pDesc-001', 't12_CFit_all':'t12_CFit_all-001', 't12_CFit_pass':'t12_CFit_pass-001', 'chart':'chart-001', 'peptides':'peptides-001', 't12_CFit_best':'t12-001'}, inplace=True)
tPeptides['show-001']=1
tBrowsing.rename(columns={ 'rank':'rank-010', 'UniProtID':'UniProtID-010', 'pDesc':'pDesc-010', 't12':'t12-010', 'mneuron':'mneuron-010', 'cneuron':'cneuron-010', 'subLoc':'subLoc-010'}, inplace=True)
tBrowsing['show-010']=10
tAbundance.rename(columns={ 'UniProtID':'UniProtID-100', 'pDesc':'pDesc-100', 'pLength':'pLength-100', 'pMass':'pMass-100', 'cneuron':'cneuron-100', 'subLoc':'subLoc-100', 'abundanceLevel':'abundanceLevel-100', 'avgCopyN':'avgCopyN-100'}, inplace=True)
tAbundance['show-100']=100

#%%
tcombined = tAbundance.copy().join(tBrowsing, how='outer').join(tPeptides, how='outer')
tcombined.shape #10869  from 9861 + 10372 + 10130 
# ['UniProtID-100', 'pDesc-100', 'pLength-100', 'pMass-100', 'subLoc-100', 'abundanceLevel-100', 'avgCopyN-100', 'show-100', 'rank-010', 'UniProtID-010', 'pDesc-010', 't12-010', 'mneuron-010', 'cneuron-010', 'subLoc-010', 'show-010', 'UniProtID-001', 'pDesc-001', 't12_CFit_all-001', 't12_CFit_pass-001', 'chart-001', 'peptides-001', 't12-001', 'show-001']

#%%
# not working for nan
# tcombined['show'] = np.nan_to_num(tcombined['show-001'])+np.nan_to_num(tcombined['show-010'])+np.nan_to_num(tcombined['show-100'])
tcombined['show'] = tcombined['show-001'].replace(np.nan,0)+tcombined['show-010'].replace(np.nan,0)+tcombined['show-100'].replace(np.nan,0)
tcombined['show'] = tcombined['show'].astype('int')
del(tcombined['show-001'])
del(tcombined['show-010'])
del(tcombined['show-100'])

#%% 
tcombined['UniProtID-100'].value_counts()
# found 52 gene had multiple rows in the tPeptides df, in 417 rows??!!
# In essence, there should be 365 less rows...
# for examples, gene  ESYT1 (uniProtID Q9BSJ8)  H1-4 (uniProtID P10412) 

#%%
tdups = tcombined.copy()
tdups.reset_index(inplace=True)
tdups.set_index(['gene','UniProtID-100'], inplace=True)
tdupsfiltered = tdups[tdups.index.duplicated()]
tdupsfiltered.shape # (365, 21)
tpeptide2chk = tdupsfiltered[['pDesc-001']]
tpeptide2chk.to_csv('peptideDuplicateGene.csv')
#%%