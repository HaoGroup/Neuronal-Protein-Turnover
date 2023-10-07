#%%
# import json
# import warnings
# import plotly.graph_objects as go
# To export static images from go, need pip install -U kaleido
# import plotly.express as px
# import numpy as np
# import seaborn as sns

# prepare data from given excel to web display

#%%
import pandas as pd
import os
import re

basename = "ProteinAbundance"

folderpath = os.path.join("..","data",basename)
# basename = "Neuron_Profile_Abundance_4website.xlsx"
# proteinAbundanceTb = pd.read_excel(os.path.join(folderpath, basename+".xlsx" ))
proteinAbundanceTb = pd.read_excel(os.path.join(folderpath, "Neuron_Profile_Abundance_4websiteShortHeader.xlsx" ), header=0, index_col=0) # make sure 'gene' is the first column
proteinAbundanceTb.columns = [ x.strip() for x in proteinAbundanceTb.columns ] # in case colheads need trimming
proteinAbundanceTb.index.name = proteinAbundanceTb.index.name.strip() # in case "gene" needs trimming
# proteinAbundanceTb.head()
# original column names :
# ['gene', 'UniProtID', 'Protein Description', 'Protein Length (#AA)', 'Protein Mass (Da)', 'Subcellular Location ', 'Abundance Level', 'Average Protein Copy Number']
# new column names :
# ['gene', 'UniProtID', 'pDesc', 'pLength', 'pMass', 'subLoc', 'abundanceLevel', 'avgCopyN']

#%% 
def txtCleanSubLoc(t):
  """
  Clean up subcellular location string
  Args:
      t (str): clean up from original string
  """
  if (type(t)==float) : return "" # most likely nan
  res = t
  # SUBCELLULAR LOCATION: Cell membrane {ECO:0000269|PubMed:19258317, ECO:0000269|PubMed:19556522, ECO:0000269|PubMed:24097981, ECO:0000269|PubMed:35974019}; Multi-pass membrane protein {ECO:0000255}. Endosome {ECO:0000269|PubMed:24097981}. 
  # Simplify to: Cell membrane. Endosome
  #
  # further samples
  # SUBCELLULAR LOCATION: Golgi apparatus membrane {ECO:0000250|UniProtKB:Q8CF82}; Multi-pass membrane protein {ECO:0000250|UniProtKB:Q8CF82}. Lysosome membrane {ECO:0000250|UniProtKB:Q8K448}; Multi-pass membrane protein {ECO:0000250|UniProtKB:Q8K448}. Late endosome membrane {ECO:0000250|UniProtKB:Q8K448}; Multi-pass membrane protein {ECO:0000250|UniProtKB:Q8K448}. Cell membrane {ECO:0000250|UniProtKB:Q8K448}. Note=Localized at cell membrane under high cholesterol levels. {ECO:0000250|UniProtKB:Q8K448}.
  # SUBCELLULAR LOCATION: Cytoplasm. Cytoplasm, cytoskeleton. Nucleus. Cell junction, focal adhesion. Note=Associates with the actin cytoskeleton near the adhesion plaques. Enters the nucleus in the presence of HESX1.
  # SUBCELLULAR LOCATION: Nucleus {ECO:0000255|PROSITE-ProRule:PRU00625}.  
  res = re.sub(r'\s*SUBCELLULAR LOCATION:\s*','',res) # beginning phrase 
  res = re.sub(r'\};','',res) # remove any Multi-pass membrane breaks 
  res = re.sub(r'\}[^\}]*Note=','',res) # remove any }***Notes 
  res = re.sub('Note=','{',res) # if there were still Note= here, make it a start to-be-deleted
  # some cases there is no ending '}', messing up the search.
  dummystr = '{xyxyxy}'
  res += dummystr # add dummystr just in case.
  res = re.sub(r'\s*\{[^\}]+\}\s*','',res) # remove additional info, from curly bracket-start to curly bracket-end, including multi-pass info
  res = re.sub(r'\.\s+','; ',res) # replace periods between locations as semi-colons. Not using commas to avoid cvs nuance.
  res = re.sub(dummystr,'',res) # remove dummystr if still here.
  res = re.sub(r'[\.\s]*$','',res) # remove ending periods. Sometimes two back-to-back periods, might have space between
  return res 

#%%[markdown]
# 
# # Normal rows
# 
# SUBCELLULAR LOCATION: Golgi apparatus membrane {ECO:0000250|UniProtKB:Q8CF82}; Multi-pass membrane protein {ECO:0000250|UniProtKB:Q8CF82}. Lysosome membrane {ECO:0000250|UniProtKB:Q8K448}; Multi-pass membrane protein {ECO:0000250|UniProtKB:Q8K448}. Late endosome membrane {ECO:0000250|UniProtKB:Q8K448}; Multi-pass membrane protein {ECO:0000250|UniProtKB:Q8K448}. Cell membrane {ECO:0000250|UniProtKB:Q8K448}. Note=Localized at cell membrane under high cholesterol levels. {ECO:0000250|UniProtKB:Q8K448}.
# 
# SUBCELLULAR LOCATION: Cytoplasm. Cytoplasm, cytoskeleton. Nucleus. Cell junction, focal adhesion. Note=Associates with the actin cytoskeleton near the adhesion plaques. Enters the nucleus in the presence of HESX1.
# 
# # Problem rows 
# 
# ABCD1 
# SUBCELLULAR LOCATION: Peroxisome membrane {ECO:0000269|PubMed:10777694, ECO:0000269|PubMed:16946495, ECO:0000269|PubMed:17609205, ECO:0000269|PubMed:18757502, ECO:0000269|PubMed:29397936}; Multi-pass membrane protein {ECO:0000255}. Mitochondrion membrane {ECO:0000269|PubMed:16946495}; Multi-pass membrane protein. Lysosome membrane {ECO:0000269|PubMed:16946495}; Multi-pass membrane protein. Endoplasmic reticulum membrane {ECO:0000269|PubMed:16946495}; Multi-pass membrane protein.
# 
# ABCG1
# SUBCELLULAR LOCATION: Endoplasmic reticulum membrane {ECO:0000269|PubMed:22042635}; Multi-pass membrane protein {ECO:0000269|PubMed:22042635}. Golgi apparatus membrane {ECO:0000269|PubMed:22042635}; Multi-pass membrane protein {ECO:0000269|PubMed:22042635}. Cell membrane {ECO:0000269|PubMed:16702602, ECO:0000269|PubMed:24576892}. Note=Predominantly localized in the intracellular compartments mainly associated with the endoplasmic reticulum (ER) and Golgi membranes.
# try parsing
# Endoplasmic reticulum membrane {ECO:0000269|PubMed:22042635Multi-pass membrane protein {ECO:0000269|PubMed:22042635}. Golgi apparatus membrane {ECO:0000269|PubMed:22042635Multi-pass membrane protein {ECO:0000269|PubMed:22042635}. Cell membrane {ECO:0000269|PubMed:16702602, ECO:0000269|PubMed:24576892}. Note=Predominantly localized in the intracellular compartments mainly associated with the endoplasmic reticulum (ER) and Golgi membranes.
# #
# ABI2
# SUBCELLULAR LOCATION: Cytoplasm {ECO:0000269|PubMed:11516653, ECO:0000269|PubMed:7590236, ECO:0000269|PubMed:8649853}. Nucleus {ECO:0000269|PubMed:7590236}.; SUBCELLULAR LOCATION: [Isoform 1]: Cell projection, lamellipodium {ECO:0000269|PubMed:11516653, ECO:0000269|PubMed:15572692}. Cell projection, filopodium {ECO:0000269|PubMed:11516653}. Cytoplasm, cytoskeleton {ECO:0000269|PubMed:15572692}. Cell junction, adherens junction {ECO:0000269|PubMed:15572692}. Note=Isoform 1 but not isoform 3 is localized to protruding lamellipodia and filopodia tips (PubMed:11516653, PubMed:15572692). Present at nascent adherens junctions, where it clusters adjacent to the tips of F-actin protrusions (PubMed:15572692). {ECO:0000269|PubMed:11516653, ECO:0000269|PubMed:15572692}.

# ABLIM2
# SUBCELLULAR LOCATION: Cytoplasm. Note=In skeletal muscle, sarcomeric or cosarcomeric localization. {ECO:0000250}.


#%%
proteinAbundanceTb.subLoc = proteinAbundanceTb.subLoc.map(txtCleanSubLoc)


#%%
# output to csv, json, excel
proteinAbundanceTb.to_csv( os.path.join(folderpath, basename+".csv"))
proteinAbundanceTb.to_excel( os.path.join(folderpath, basename+".xlsx"))
proteinAbundanceTb.reset_index().to_json( os.path.join(folderpath, basename+".json"), orient="records")

# basename = "ProteinT12Browsing"

# %%

# %%
