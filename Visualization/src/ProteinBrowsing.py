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

basename = "ProteinT12Browsing"

folderpath = os.path.join("..","data",basename)
# basename = "ProteinT12Browsing.xlsx"
# proteinBrowseTb = pd.read_excel(os.path.join(folderpath, basename+".xlsx" ))
proteinBrowseTb = pd.read_excel(os.path.join(folderpath, "ProteinT12Browsing_shortHeader.xlsx" ), header=0, index_col=0) # make sure 'rank' is the first column
proteinBrowseTb.columns = [ x.strip() for x in proteinBrowseTb.columns ] # in case colheads need trimming
proteinBrowseTb.index.name = proteinBrowseTb.index.name.strip() # in case "rank" needs trimming
# proteinBrowseTb.head()
# original column names :
# Rank	Gene	UniProt ID	Protein Description	Average Half-Life (Days)	Motor Neuron (Days) 	Cortical Neuron (Days) 	Subcellular Location 
# ['Rank', 'Gene', 'UniProt ID', 'Protein Description', 'Average Half-Life (Days)', 'Motor Neuron (Days) ', 'Cortical Neuron (Days) ', 'Subcellular Location ']
# new column names :
# ['rank', 'gene', 'UniProtID', 'pDesc', 't12', 'mneuron ', 'cneuron', 'subLoc']

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


#%%
proteinBrowseTb.subLoc = proteinBrowseTb.subLoc.map(txtCleanSubLoc)


#%%
# output to csv, json, excel
proteinBrowseTb.to_csv( os.path.join(folderpath, basename+".csv"))
proteinBrowseTb.to_excel( os.path.join(folderpath, basename+".xlsx"))
proteinBrowseTb.reset_index().to_json( os.path.join(folderpath, basename+".json"), orient="records")

# basename = "ProteinT12Browsing"

# %%
# Create Protein Halflife rank plot
import plotly.graph_objects as go

def createRankPlot(df) -> None: # From ProteinTurnover class
  fig = go.Figure()
  # dfplot = df[["rank","t12"]]
  df.apply(singleProteinDatumPlot, fig = fig, axis=1 )

  # add horizontal double arrow
  arrhead, arrsize, arrwidth, arrcolor = (3,1,3,'blue')
  fig.add_annotation(x=50,y=65,ax=1000,ay=65,xref='x',yref='y',axref='x',ayref='y',text='faster', yanchor="bottom", xanchor="left", showarrow=True,arrowhead=arrhead,arrowsize=arrsize,arrowwidth=arrwidth,arrowcolor=arrcolor)
  fig.add_annotation(x=50,y=65,ax=9500,ay=65,xref='x',yref='y',axref='x',ayref='y',text='', showarrow=True,arrowhead=arrhead,arrowsize=arrsize,arrowwidth=arrwidth,arrowcolor=arrcolor)
  fig.add_annotation(x=10000,y=65,ax=8500,ay=65,xref='x',yref='y',axref='x',ayref='y',text='slower turnover', yanchor="bottom", xanchor="right", showarrow=True,arrowhead=arrhead,arrowsize=arrsize,arrowwidth=arrwidth,arrowcolor=arrcolor)
  # fig.add_hline(y=65, x0=1000, x1=9000, line_width=1, line_dash="dash", line_color="black") 

  fig.update_layout( 
    title={
      'text': 'Protein Half-life rank',
      'x': 0.45,
      'xanchor': 'center'
      }, 
    xaxis_title='Rank of Protein', 
    yaxis_title='Protein Half-life (d)',
    font=dict(
        family="Garamond",
        size=16,
        color="Black" # "RebeccaPurple"
    )
  )
  
  fig.show()
  incPlotlyJs = "plotly.min.js"  # False # add cdn plotly js on the html directly to minimize size of each html # symbolic link
  fig.write_html( "proteinRank.html", include_plotlyjs=incPlotlyJs )
  
  return

def singleProteinDatumPlot(proteinrow, fig) -> None: # From ProteinTurnover class
  # protein = proteinrow.name[0]
  gene = proteinrow['gene']
  specs = dict(name=gene, connectgaps=False, mode='markers',  showlegend=False, legendgroup = 't12rank')
  markeropt = dict(color="black", symbol='circle', size=4)
  add1goTrace(fig, x=[int(proteinrow.name)], y=[proteinrow["t12"]], specs=specs, markeropt=markeropt)
  return # fig

def add1goTrace(fig, x, y, specs=dict(), lineopt=dict(), markeropt=dict()): # From ProteinTurnover class
  """
  Add 1 trace to go.Figure object
  Args:
      fig (Plotly go): Plotly graph_object
      x (list): x data
      y (list): y data
      specs (dict, optional): dict(mode, name, showlegend, connectgaps). Defaults to {}. name can have html tags.
      lineopt (dict, optional): line options. {color, width}. Defaults to {}.
      markeropt (dict, optional): marker options. {symbol, size, color}. Defaults to {}.
  """
  # check lineopt standards
  # check markeropt standards
  # https://plotly.com/python/marker-style/ # hourglass, cross, x-thin, ...
  mode = specs['mode'] if (specs.__contains__('mode') and specs['mode']) else 'lines'
  name = specs['name'] if (specs.__contains__('name') and specs['name']) else '-'
  connectgaps = specs['connectgaps'] if (specs.__contains__('connectgaps') and isinstance(specs['connectgaps'], bool) ) else True

  showlegend = specs['showlegend'] if ( specs.__contains__('showlegend') and isinstance(specs['showlegend'], bool) ) else False
  legendgroup = specs['legendgroup'] if (specs.__contains__('legendgroup') and specs['legendgroup']) else specs['name']
  
  # if ( (not yerroropt['array'] is None) and yerroropt['visible']):
  #   fig.add_trace(go.Scatter( mode=mode, x=x, y=y, error_y=yerroropt, showlegend=showlegend, legendgroup=legendgroup, name=name, connectgaps=connectgaps, line=lineopt, marker=markeropt))
  # else: # add more conditions/scenarios here if needed
  #   fig.add_trace(go.Scatter( mode=mode, x=x, y=y,showlegend=showlegend, legendgroup=legendgroup, name=name, connectgaps=connectgaps, line=lineopt, marker=markeropt))
  fig.add_trace(go.Scatter( mode=mode, x=x, y=y,showlegend=showlegend, legendgroup=legendgroup, name=name, connectgaps=connectgaps, line=lineopt, marker=markeropt))
    
  # if showlegend: fig.update_layout( showlegend=showlegend, legend=legendOpts )
  return fig
    
    
# %%
createRankPlot(proteinBrowseTb)
# %%
