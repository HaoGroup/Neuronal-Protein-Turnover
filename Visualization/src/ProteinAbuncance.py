#%%
# import json
# import re
# import warnings
# import plotly.graph_objects as go
# To export static images from go, need pip install -U kaleido
# import plotly.express as px
# import numpy as np
# import seaborn as sns

# prepare data from given excel to web display

#%%
import numpy as np
import pandas as pd
import os
import re

basename = "ProteinAbundance"

folderpath = os.path.join("..","data",basename)
# basename = "Neuron_Profile_Abundance_4website.xlsx"
# proteinAbundanceTb = pd.read_excel(os.path.join(folderpath, basename+".xlsx" ))
proteinAbundanceTb = pd.read_excel(os.path.join(folderpath, "Neuron_Profile_Abundance_4websiteShortHeader.xlsx" ))
# proteinAbundanceTb.head()
# original column names :
# ['GeneID', 'UniProtID', 'Protein Description', 'Protein Length (#AA)', 'Protein Mass (Da)', 'Subcellular Location ', 'Abundance Level', 'Average Protein Copy Number']
# new column names :
# ['GeneID', 'UniProtID', 'pDesc', 'pLength', 'pMass', 'subLoc', 'abundanceLevel', 'avgCopyN']

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
  res = re.sub(r'\}[^\}]*Note=','',res) # remove any Notes 
  res = re.sub(r'\s*\{[^\}]+\}\s*','',res) # remove additional info, including multi-pass info
  res = re.sub(r'\.\s+','; ',res) # replace periods between locations as semi-colons. Not using commas to avoid cvs nuance.
  
  return res 

#%%
proteinAbundanceTb.subLoc = proteinAbundanceTb.subLoc.map(txtCleanSubLoc)



#%%
# output to csv, json, excel
proteinAbundanceTb.to_csv( os.path.join(folderpath, basename+".csv"))
proteinAbundanceTb.to_excel( os.path.join(folderpath, basename+".xlsx"))
proteinAbundanceTb.to_json( os.path.join(folderpath, basename+".json"), orient="records")

# basename = "ProteinT12Browsing"

# %%

# %%