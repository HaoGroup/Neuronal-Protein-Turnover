#%%
import pandas as pd
import os
import json
import re
import warnings
import plotly.graph_objects as go
# To export static images from go, need pip install -U kaleido
import plotly.express as px
import numpy as np
import ExpoDecayFit as edf
# from scipy.stats import hmean # harmonic mean
# import matplotlib.pyplot as plt
# import seaborn as sns

# Create class ProteinTurnover to handle and organize the various parts of the study
# v3 - getting results table first, then (optional) create plots. Also use .apply() whenever possible
# v2 - switching to plotly from matplotlib
# added column "time" instead of using the index to allow better/easier plots 
#
class ProteinTurnover:
  """
  Neuronal Protein Turnover Study
  Ingest data
  Process various computations and data fit
  Creates plots
  """
  def __init__(self, datafiles=dict(raw=None,peptides=None, proteins=None)) -> None:
    """_summary_
    Args:
        datafiles (dict): { 'raw': file location, 'peptides': df file location, 'proteins': df file location }
    """

    self.__yAxisName = 'Relative Abundance'
    self.__yUnit = '%'
    self.__xAxisName = 'Time'
    self.__xUnit = 'day'
    self.__xvalues = ['1', '2', '4', '6']
    # Match with ExpoDecayFit module:
    self.__statsTypes = ('b','t12','r2') # keeping track of three results/statics from exponential fit
    self.__modelTypes = ("CFit",) # use only Curve-fit, 20230713
    # self.__modelTypes = ("LnLM1", "LnLM2", "CFit") # trying three different models in ExpoDecayFit: LnLM1 (sklearn LM), LnLM2 (statsmodels LM), CFit (scipy curve_fit)

    self.__compoIndexPeptide = ['Protein', 'Peptide'] # basic composite index used in Peptide and related DFs
    self.__compoIndexGene = ['Protein', 'Gene'] # basic composite index used in other dfs.
    self.__maxNAcnt = 1 # Each series only has at most 4 data points at day = 1,2,4,6.  Only allow at most 1 missing to be considered good peptide data
    self.__supportMin = len(self.__xvalues) - self.__maxNAcnt
    self.__r2cutoff = 0.8 # set lower bound for r-square to be plotted
    self.__colorPalette = px.colors.qualitative.Dark24 
    self.__colorMax = len(self.__colorPalette)
    self.__colorInd = 0 # initialize
    self.df_Peptides = None # initialize, cleaned and re-structured df for analysis ['PG.ProteinGroups', 'PG.Genes', 'PG.ProteinDescriptions', 'Peptide', 0, 1, 2, 4, 6, 'k_results', 'Protein_Turnover_from_k', 'residuals']
    self.df_Proteins = None # after all calcuations done, prepare for protein half life chart 
    self.__ingestData(datafiles=datafiles) # self.__setDfPeptides(), # self.__setDfProteins() # called inside __ingestData # set df_Peptides, including sorting, and df_Proteins
    return
  
  def __ingestData(self, datafiles) -> None:
    """
    Args:
        datafiles (dict): { 'raw': file location, 'peptides': df file location, 'proteins': df file location }
    Returns: None
    """
    datafiles = self.__setArgDatafiles(datafiles)
    raw, peptides, proteins = datafiles.values()
    
    if (raw): 
      df = pd.read_excel(raw) if (raw[-5:] == '.xlsx' or raw[-4:] == '.xls') else pd.read_csv(raw) if (raw[-4:] == '.csv') else None
      
      # clean up column names, avoid dots in names
      df.columns =  [ str(col) if str(col).isnumeric() else col.replace('.','_') for col in df.columns.values ]
      # 20230620 data file no longer has day0 columns. Re-create here:
      if (not '0' in df.columns): df['0'] = 1
      # df['iMN_Day0_Light_Relative_Abundance'] = 1
      # df['iMN_Day0_Heavy_Relative_Abundance'] = 0
      
      # Keep only the light series, drop heavy
      df.drop(list(df.filter(regex = 'Heavy_Relative_Abundance')), axis = 1, inplace = True)
      
      # import re  # rename columns as 0, 1, 2, 4, 6 for days if needed
      # df.rename(columns=lambda x: re.sub(r'iMN_Day(\d).*_Relative_Abundance$', r'\1', x), inplace = True)
      self.__setDfPeptides(df_init=df)
    elif (peptides):
      self.df_Peptides = pd.read_csv(peptides, index_col=self.__compoIndexPeptide)
      
    if not self.chkDfPeptidesIndexUnique(): warnings.warn(f"Warning: df_Peptides index ({self.__compoIndexPeptide}) not unique!!")

    # get df_Proteins table either from given file, or from setDfProteins function.
    if (proteins):
      self.df_Proteins = pd.read_csv(proteins, index_col=self.__compoIndexGene)
      self.df_Proteins['peptides']=self.df_Proteins['peptides'].apply( lambda x: json.loads(x.replace("'",'"')) ) # csv saved the peptides columns as strings, while it should be dict(). Json.loads() expect key to be in double quote, not single.
    else:
      self.__setDfProteins() # set df_Proteins when df_Peptides all set

    return None
  
  def __setDfPeptides(self, df_init) -> None:
    self.df_Peptides = df_init.set_index(self.__compoIndexPeptide)  # set composite index 

    self.df_Peptides['chart'] = 0 # keep track of whether the peptide has enough data points to make a plot. default to 0
    self.df_Peptides['support'] = len(self.__xvalues) -  self.df_Peptides[self.__xvalues].isna().sum(axis=1) # support is the number of non-NA data points for the peptide data
    newcols = [ s+'_'+m for s in self.__statsTypes for m in self.__modelTypes ] # for the modelcurve-fit statistics results
    # newcols += [ 'proteinT12', 'proteinT12est' ] # protein level half-lives, first one with only high quality data (at most one missing time data, and r^2 > cutoff of 0.8), second one with all data regardless. The t1/2 is calculated using the harmonic mean of the peptide t1/2s, or the regular mean of their decay constants.
    self.df_Peptides[newcols] = np.nan
    
    self.__setPeptidesStats() # get all peptides stats from ExpoDecayFits
    
    # sort them here once all t12, r2, and chart values are set
    # need two temporary columns to help sort
    self.df_Peptides['chart_sort']=self.df_Peptides['chart']+np.around(self.df_Peptides['support']/2) - 1 # from 0,1 becomes 0,1,2 for sorting only
    self.df_Peptides['r2_sort']=self.df_Peptides['chart_sort']*self.df_Peptides['r2_CFit'] # for sorting only
    self.df_Peptides.sort_values(by=['Gene','chart_sort','r2_sort','support','Peptide'], ascending=[True, False, False, False, True], inplace=True)
    del(self.df_Peptides['chart_sort'])
    del(self.df_Peptides['r2_sort'])
    
    return
  
  def __setPeptidesStats(self) -> None:
    self.df_Peptides.apply(self.__set1PeptideStats, axis = 1) # get 1 peptide stats from ExpoDecayFits
    return
  
  def __set1PeptideStats(self, peptiderow) -> None:
    thissupport = peptiderow['support']
    proteinname, peptidename  = peptiderow.name
    model = edf.ExpoDecayFit(peptiderow, xAxisName = self.__xAxisName, xvalues = self.__xvalues, modelTypes=self.__modelTypes, statsTypes=self.__statsTypes ) # model.startx, starty, samplexs, sampleys
    bs, t12s, r2s = [ model.modelsummary.loc[t,:] for t in self.__statsTypes ]
    stats = dict( b=bs, t12=t12s, r2=r2s )
    if ( thissupport > self.__supportMin -1 and (stats['r2']>self.__r2cutoff).any() ) : self.df_Peptides.loc[ (proteinname, peptidename) , 'chart' ] = 1 # instead of .any, consider using .all instead. When multiple models are used, then critical
    # save results in df_Peptides
    for m in self.__modelTypes:
      for s in self.__statsTypes:
        self.df_Peptides.loc[( proteinname, peptidename ), s+'_'+m] = stats[s][m]
    
    return
  
  def __setDfProteins(self):
    """
    Generate df_Proteins dataframe after all results were obtained. 
    df_Proteins will have structure that supports the web data 
    """
    df = self.df_Peptides.drop(list(self.df_Peptides.filter(regex = '^\d$')), axis = 1) # do not keep the time data, just prtn, gene, desc, bs, t12s, r2s.
    
    pepCol = 'peptides' # new column name for all peptide info per gene, for web json
    # newHeaders = {"Protein": "prtn", "Gene": "gene", "Protein_Description": "desc"}
    statsHeaders = { s: [ s+'_'+m for m in self.__modelTypes ] for s in self.__statsTypes } # the stats metrics that we look for, in particular, b-decay constant, and t12-halfLife #  'b_LnLM1', 'b_LnLM2', 'b_CFit' #  't12_LnLM1', 't12_LnLM2', 't12_CFit' # 'r2_LnLM1', 'r2_LnLM2', 'r2_CFit'
    modelMetricHeaders = sum( statsHeaders.values(), [])  #  'b_LnLM1', 'b_LnLM2', ... 't12_LnLM1' ... 'r2_CFit'
    colHeads = [ self.__compoIndexPeptide[1] ]+modelMetricHeaders + ['support'] # add 'Peptide' column (in df_Peptides) to the modelMetricHeaders, plus others

    df.reset_index(inplace=True)
    # df.rename(newHeaders, axis=1, inplace=True)
    df_gb = df.groupby( list(self.__compoIndexGene + ['Protein_Description'] ), as_index=True, dropna=False)  # groupby object
    df_res = df_gb[statsHeaders['b']].agg('mean')  # first find mean b-values decay constants
    df_res[ [ t+'_all' for t in statsHeaders['t12'] ] ] = np.log(2)/df_res[statsHeaders['b']]  # next find half lives from bs. Essentially, we are finding the harmonic mean of half lives, using all peptides
    df_res.rename( { 'b_'+m : 'b_'+m+'_all' for m in self.__modelTypes } , axis=1, inplace=True) # rename b_CFit to b_CFit_all, etc
    
    # the good peptides only
    df_res[ [ b+'_pass' for b in statsHeaders['b'] ] ]  = df_gb.apply(lambda x: x[x['chart']>0][statsHeaders['b']].mean()) # find mean of only the good/passed peptides
    df_res[ [ t+'_pass' for t in statsHeaders['t12'] ] ] = np.log(2)/df_res[ [ b+'_pass' for b in statsHeaders['b'] ]  ]  # next corresponding half lives

    df_res['chart'] = df_gb['chart'].agg('sum') # find total number of peptides have charts
    # next, get list of peptides and their metrics from different models into a dict()
    df_res[pepCol] = df_gb[colHeads].apply(lambda g: { h: list(g[h]) for h in colHeads} ) # adding the resulting pd series to df_res
    
    # set one more column for protein plot, with t12_pass value if available, or else use t12_all
    df_res[ [t+'_best' for t in statsHeaders['t12']] ] = pd.isna(df_res[[t+'_pass' for t in statsHeaders['t12']]]).values * df_res[ [t+'_all' for t in statsHeaders['t12']] ].values + pd.notna(df_res[[t+'_pass' for t in statsHeaders['t12']]]).values * np.nan_to_num(df_res[[t+'_pass' for t in statsHeaders['t12']]])
    
    df_res.sort_values(by='t12_'+self.__modelTypes[0]+'_best', inplace=True)
    df_res.reset_index(inplace=True)
    df_res['rank'] = df_res.index +1
    
    self.df_Proteins = df_res.set_index(self.__compoIndexGene)

    return

  def chkDfPeptidesIndexUnique(self): return self.df_Peptides.index.is_unique
  
  def exportResults(self, filename="", format='json'):
    """
    Export df_Peptides table for web dev
    Args:
        format (str, optional): _description_. Defaults to 'json', where only df_Proteins is exported. If csv, both df_Peptides and df_Proteins are exported.
    """
    filename = filename if filename else "result"
    if format=='csv': 
      self.df_Peptides.drop(list(self.df_Peptides.filter(regex = '^\d$')), axis = 1).to_csv(os.path.join("../data/",filename+"_peptides.csv"))
      self.df_Proteins.to_csv(os.path.join("../data/",filename+"_protein.csv"))
      return
    
    # default json, will have protein as key, values will be Gene, description, decay constants, half-lives, rsquared, and peptide list
    newHeaders = {"Protein": "prtn", "Gene": "gene", "Protein_Description": "desc"}
    # newHeaders = {"Protein": "prtn", "Gene": "gene", "Protein_Description": "desc", "proteinT12" : "proteinT12" , "proteinT12est" : "proteinT12est" }
    df = self.df_Proteins.copy()
    df.reset_index(inplace=True)
    df.rename(newHeaders, axis=1, inplace=True)
    
    print(f'resulting results has dim: {df.shape}')
    print(df.head(2))
    df.to_json( os.path.join("../data/",filename+"_protein.json"), orient="records")
    return
  
  
  def __setNextColorInd(self, n=1) -> None:
    self.__colorInd = (self.__colorInd + n) % self.__colorMax
    return
  
  def __setArgDatafiles(self, datafiles=dict(raw=None,peptides=None, proteins=None) ):
    """
    setting generic argument datafiles for __ingestData
    Args:
        datafiles (dict): { 'raw': file location, 'peptides': df file location, 'proteins': df file location }
    Returns:
        dict: { raw, peptides, proteins }
    """
    res = datafiles.copy()
    res['raw'] = os.path.join(res['raw']) if (res.__contains__('raw') and res['raw']) else None
    res['peptides'] = os.path.join(res['peptides']) if (res.__contains__('peptides') and res['peptides']) else None
    res['proteins'] = os.path.join(res['proteins']) if (res.__contains__('proteins') and res['proteins']) else None
    return res
  
  def __setArgLabels(self, labels=dict(x=None, y=None, title=None, fontfamily=None, size=None) ):
    """
    setting generic keyword argument labels into x, y, and title
    Args:
        labels (dict, optional): x-, y-lables, title, fontfamily and size of chart/text. Defaults to dict(x=None, y=None, title=None).
    Returns:
        dict: { x, y, title }
    """
    res = labels.copy()
    res['x'] = res['x'] if (res.__contains__('x') and res['x']) else self.__xAxisName.capitalize() + ' ('+ self.__xUnit + ')'
    res['y'] = res['y'] if (res.__contains__('y') and res['y']) else self.__yAxisName.capitalize() + ' ('+ self.__yUnit + ')'
    res['title'] = res['title'] if (res.__contains__('title') and res['title']) else None
    res['fontfamily'] = res['fontfamily'] if (res.__contains__('fontfamily') and res['fontfamily']) else "Garamond"
    res['size'] = res['size'] if (res.__contains__('size') and res['size']) else 16
    return res
  
  def __setArgLineopt(self, lineopt=dict(color=None, width=None) ):
    """
    setting generic keyword argument lines into show, solid, color, and width
    Args:
        lineopt (dict, optional): color (str), width (float), ..
    Returns:
        dict: { color, width, ... }
    """
    res = lineopt.copy()
    res['color'] = res['color'] if (res.__contains__('color') and res['color']) else '#000000'
    res['width'] = res['width'] if (res.__contains__('width') and res['width']) else 2
    # preserving the rest of attributes
    return res

  def __setArgMarkeropt(self, markeropt=dict(color=None, symbol=None, size=None) ):
    """
    setting generic keyword argument markers into color, symbol, size
    Args:
        markeropt (dict, optional): color (str), symbol (str), size (float), ...
    Returns:
        dict: { color, symbol, size, ... }
    """
    res = markeropt.copy()
    # https://plotly.com/python/marker-style/ # hourglass, x-thin, ...
    res['symbol'] = res['symbol'] if (res.__contains__('symbol') and res['symbol']) else 'circle'
    res['color'] = res['color'] if (res.__contains__('color') and res['color']) else '#000000'
    res['size'] = res['size'] if (res.__contains__('size') and res['size']) else 4
    return res
  
  def __setArglegendOpts(self, legendOpts=dict(title_font_family=None, font=dict(size=None)) ):
    """
    setting generic keyword argument legend with title_font_family, size
    Args:
        legendOpts (dict, optional): legendgroup (str), title_font_family (str), font (dict), ...
    Returns:
        dict: { title_font_family, font, ... }
    """
    res = legendOpts.copy()
    res['title_font_family'] = res['title_font_family'] if (res.__contains__('title_font_family') and res['title_font_family']) else 'Garamond'
    res['font'] = res['font'] if (res.__contains__('font') and res['font']) else dict(size=8)
    res['font']['size'] = res['font']['size'] if (res['font'].__contains__('size') and res['font']['size']>0) else 7
    return res
  
  # def __setArgYerroropt(self, yerroropt=dict(type=None, array=None, arrayminus=None, symmetric=None, visible=None) ):
  #   """
  #   setting generic keyword argument yerror with type, array, arrayminus, symmetric, visible, value
  #   Args:
  #       yerroropt (dict, optional): type (str), array (list), arrayminus (list), symmetric (bool), visible (bool), value (float, for percent type)...
  #   Returns:
  #       dict: { type, array, ... }
  #   """
  #   # fig = go.Figure(data=go.Scatter(
  #   #     x=[0, 1, 2],
  #   #     y=[6, 10, 2],
  #   #     error_y=dict(
  #   #         type='data', # value of error bar given in data coordinates
  #   #         array=[1, 2, 3],
  #   #         visible=True)
  #   # ))
  #   res = yerroropt.copy()
  #   res['type'] = res['type'] if (res.__contains__('type') and res['type']) else 'data' # 'data', 'percent'
  #   res['symmetric'] = res['symmetric'] if ( res.__contains__('symmetric') and isinstance(res['symmetric'], bool) ) else True
  #   res['visible'] = res['visible'] if ( res.__contains__('visible') and isinstance(res['visible'], bool) ) else True
  #   res['array'] = res['array'] if (res.__contains__('array') and (not res['array'] is None) ) else None
  #   res['arrayminus'] = res['arrayminus'] if (res.__contains__('arrayminus') and (not res['arrayminus'] is None) ) else None
  #   res['value'] = res['value'] if (res.__contains__('value') and (not res['value'] is None) and res['type']=='percent') else None
  #   return res
  
  # def __setArgTrendlines(self, trendlines=dict(show=False, solid=True, color=None, width=None) ):
  #   """
  #   setting generic keyword argument trendlines into show, solid, color, and width
  #   Args:
  #       trendlines (dict, optional): show, solid (bool), color (str), width (float) 
  #   Returns:
  #       dict: { show, solid, color, width }
  #   """
  #   show = trendlines['show'] if (trendlines.__contains__('show') and isinstance(trendlines['show'], bool) ) else False
  #   solid = trendlines['solid'] if (trendlines.__contains__('solid') and isinstance(trendlines['solid'], bool) ) else True
  #   color = trendlines['color'] if (trendlines.__contains__('color') and trendlines['color']) else '#00FE35' # trendline use '#00FE35', which is Light24[1] color, lime green
  #   width = trendlines['width'] if (trendlines.__contains__('width') and trendlines['width']) else 4
  #   return dict(show=show, solid=solid, color=color, width=width)
  
  def __setArgSaveFigOpts(self, saveFigOpts=dict(savePlot=False, showPlot=True , folder=None) ):
    """
    setting generic keyword argument saveFigOpts into savePlot and folder
    Args:
        saveFigOpts (dict, optional): savePlot (binary) , showPlot (binary) and folder (str). Defaults to dict(savePlot=False, folder=None).
    Returns:
        dict: { savePlot, folder }
    """
    res = saveFigOpts.copy()
    res['savePlot'] = res['savePlot'] if ( res.__contains__('savePlot') and isinstance(res['savePlot'], bool) ) else False
    res['showPlot'] = res['showPlot'] if ( res.__contains__('showPlot') and isinstance(res['showPlot'], bool) ) else True
    res['folder'] = res['folder'] if ( res.__contains__('folder') and res['folder'] ) else './'
    return res
  
  def proteinHalflifeChart(self, saveFigOpts = dict(savePlot=True, showPlot=True, folder='../media/plots/proteinRank/')) -> None:
    saveFigOpts=self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    d = self.df_Proteins
    df = d.drop(list(d.filter(regex = '(_all|_pass)')) + ['peptides', 'Protein_Description'] , axis = 1)
    # df.plot(x="rank", y="t12_"+self.__modelTypes[0]+"_best")
    fig = go.Figure()
    df.apply(self.__singleProteinDatumPlot, fig=fig, axis=1)
    
    fig.update_layout( 
      title={
        'text': 'Protein Half-life rank',
        'x': 0.45,
        'xanchor': 'center'
        }, 
      xaxis_title='Rank', 
      yaxis_title='Average half-lives (day)',
      font=dict(
          family="Garamond",
          size=16,
          color="Black" # "RebeccaPurple"
      )
    )

    if saveFigOpts['savePlot']:
      options = saveFigOpts.copy() 
      self.__savePlot(options, fig, "proteinHalflifeRank")
    if saveFigOpts['showPlot']: fig.show()
    
    return 
  
  def __singleProteinDatumPlot(self, proteinrow, fig ) -> None:
    model=self.__modelTypes[0]
    # protein = proteinrow.name[0]
    gene = proteinrow.name[1]
    specs = dict(name=gene, connectgaps=False, mode='markers',  showlegend=False, legendgroup = 't12rank')
    markeropt = dict(color="black", symbol='circle', size=4)
    self.__add1goTrace(fig, x=[proteinrow["rank"]], y=[proteinrow["t12_"+model+"_best"]], specs=specs, markeropt=markeropt)
    return # fig
  
  # def __add1goTrace(self, fig, x, y, specs=dict(), lineopt=dict(), markeropt=dict(), legendOpts=dict(), yerroropt=dict() ):
  def __add1goTrace(self, fig, x, y, specs=dict(), lineopt=dict(), markeropt=dict(), legendOpts=dict() ):
    """
    Add 1 trace to go.Figure object
    Args:
        fig (Plotly go): Plotly graph_object
        x (list): x data
        y (list): y data
        error_y (dict, optional): dict(type, array, arrayminus, visible, symmetric, value)
        specs (dict, optional): dict(mode, name, showlegend, connectgaps). Defaults to {}. name can have html tags.
        lineopt (dict, optional): line options. {color, width}. Defaults to {}.
        markeropt (dict, optional): marker options. {symbol, size, color}. Defaults to {}.
        error_y (dict, optional): error bars for each datapoint. Defaults to None.
    """
    # yerroropt = self.__setArgYerroropt(yerroropt=yerroropt)
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
  
  def abundancePlotPgAll(self, labels=dict(), saveFigOpts=dict()): 
    """
    Args:
        labels (dict, optional): x-, y-labels and title. Defaults to empty dictionary
        saveFigOpts (dict, optional): savePlot (binary) and folder (str). Defaults to dict(savePlot=False, folder=None).
    return: None
    """
    # assumes labels and saveFigOpts are in the right forms.
    labels = self.__setArgLabels(labels=labels)
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    self.df_Proteins.apply(self.abundancePlot1Pg , labels=labels, saveFigOpts=saveFigOpts , axis=1 )

    return
  
  def abundancePlot1Pg(self, prtnrow, labels=dict(), saveFigOpts = dict(), legendOpts=dict()  ):
    """
    Args:
        prtnGrp (tuple): the (ProteinGrp, Gene) being processed
        labels (dict, optional): x-, y-labels and title. Defaults to empty dictionary
        saveFigOpts (dict, optional): savePlot (binary) and folder (str). Defaults to dict(savePlot=False, folder=None).
        legendOpts (dict, optional): legendgroup (str), title_font_family (str), font (dict), ...
    return: None
    """
    labels = self.__setArgLabels(labels=labels)
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    legendOpts = self.__setArglegendOpts(legendOpts=legendOpts)
    t12best = prtnrow['t12_'+self.__modelTypes[0]+'_best'] # use the default model best estimate of peptide half lives harmonic mean
    # xrange = round( min( max(6.5, 1.6*t12best) , 10) ) # no more than 10, between 6.5 and 10. If t12 is close, show 1.6*t12
    xrange = 0.5 * round( min( max(13, 3.2*t12best) , 20) ) # no more than 10, between 6.5 and 10. If t12 is close, show 1.6*t12
    sampleN = 300
    samplexs = np.linspace(start=0, stop=xrange, num = sampleN)

    df1prtn = self.df_Peptides.loc[prtnrow.name[0],:] # subset only the protein info in df_Peptides
    
    if df1prtn.shape[0] < 1: return None # nothing to plot
    # prtnrow.name = ('Protein', 'Gene')
    gene = prtnrow.name[1] if type(prtnrow.name[1]) == str else 'Protein Group: '+prtnrow.name[0]
    
    # labels['title'] = f'Gene: {gene}'
    labels['title'] = gene
    # cpalette = self.__colorPalette # trendline use '#00FE35', which is Light24[1] color, lime green
    fig = go.Figure()    
    
    df1prtn.apply(self.abundancePlot1Peptide, fig=fig, prtnGrp=prtnrow.name, samplexs = samplexs, axis=1 )
    # _ = self.abundancePlotAllPeptides(df1prtn, prtnGrp=prtnrow.name, labels=labels, saveFigOpts=saveFigOpts)
    
    if len(fig.data) < 1 : return #  if nothing showing, skip

    # show half life if within range
    # if t12best < xrange: fig.add_vline(x=t12best, line_width=1, line_dash="dash", line_color="black", annotation_text="&nbsp;<b>t<sub>½</sub></b> = "+str(t12best.__round__(2)), annotation_position='bottom right')
    
    fig.update_layout( 
      title={
        'text': labels['title'],
        'x': 0.45,
        'xanchor': 'center'
        }, 
      xaxis_title=labels['x'], 
      yaxis_title=labels['y'],
      legend_title="Peptide",
      legend_tracegroupgap=1, 
      font=dict(
          family=labels['fontfamily'],
          size=labels['size'],
          color="Black" # "RebeccaPurple"
      ),
      legend=legendOpts
    )
    
    proteingenename = prtnrow.name[1] if type(prtnrow.name[1]) == str else prtnrow.name[0]
    
    if saveFigOpts['savePlot']:
      options = saveFigOpts.copy() 
      options['folder'] += 'htmls/peptideLevel/'
      self.__savePlot(options, fig, "RelAbundance_Gene-"+ proteingenename +"-peptides")
    if saveFigOpts['showPlot']: fig.show()

    # Now plot protein level average with trendline
    # df = self.df_Peptides.loc[prtnrow.name[0],:] # df most likely have changed by abundancePlotAllPeptides 
    # labels['title'] = f'Protein level: {gene}' 
    labels['title'] = gene 
    _ = self.abundancePlotProteinLevel(xrange=xrange, xs = samplexs, t12 = t12best, prtnGrp=prtnrow.name, labels=labels, saveFigOpts=saveFigOpts)

    return 

  def abundancePlot1Peptide(self, peptiderow, fig, prtnGrp, samplexs, markeropt=dict()):
    """
    Create curve fitting for each peptide, find the decay constant, half life, and r-squared; update the results table, and create plot if applicable
    Args:
        fig (plotly fig): to be added with this new line plot
        peptiderow (Pandas series): Pandas series data for one peptide
        prtnGrp (str): Protein group name/id
        peptide (str): Peptide name
        lineopt (dict, optional): show, solid (bool), color (str), width (float) 
        markeropt (dict, optional): show (bool), symbol (str), size (float) 
        legendOpts (dict, optional): legendgroup (str), title_font_family (str), font (dict),...
    Return: 
        plotly graph object GO
    """
    peptide = peptiderow.name
    lineopt = self.__setArgLineopt(lineopt=dict( color=self.__colorPalette[ self.__colorInd ], width=2))
    markeropt = self.__setArgMarkeropt(markeropt=markeropt)
    markeropt['color'] = lineopt['color']
    # legendOpts = self.__setArglegendOpts(legendOpts=legendOpts)
    # thissupport = self.df_Peptides.loc[ (prtnGrp[0],peptide) , 'support' ]
    thischart = self.df_Peptides.loc[ (prtnGrp[0],peptide) , 'chart' ]
    modelChoice = self.__modelTypes[0]
    b = peptiderow['b_'+modelChoice]
    
    # only if support > threshold, and r-square > cutoff, then show graph
    if thischart < 1: return # fig # not eligible for a plot 

    # model lines
    sampleys = 100 * np.exp(-b*samplexs)
    specs = dict(name=peptide+f' t = {peptiderow["t12_"+ modelChoice].__round__(1)}d', connectgaps=False, mode='lines',  showlegend=False, legendgroup = peptide)
    markeropt['size'] = 2
    # legendOpts['showlegend'] = True 
    specs['showlegend'] = True 
    fig = self.__add1goTrace(fig, x=tuple(samplexs), y=tuple(sampleys), specs=specs, lineopt=lineopt, markeropt=markeropt )
    # fig = self.__add1goTrace(fig, x=tuple(samplexs)+(np.nan,)+tuple(samplexs ), y=tuple(sampleys)+(np.nan,)+tuple(100-sampleys) , specs=specs, lineopt=lineopt, markeropt=markeropt )

    specs = dict(connectgaps=False, mode='lines', showlegend=False, legendgroup = peptide)
    fig = self.__add1goTrace(fig, x=tuple(samplexs), y=tuple(100-sampleys), specs=specs, lineopt=lineopt, markeropt=markeropt )
    
    # data points
    # specs = dict(name=peptide+f' t = {stats["t12"][modelChoice].__round__(1)}d', connectgaps=False, mode='markers',  showlegend=False, legendgroup = peptide)
    specs = dict(name=peptide, connectgaps=False, mode='markers',  showlegend=False, legendgroup = peptide)
    markeropt['size'] = 4
    xvals = ['0']+self.__xvalues
    fig = self.__add1goTrace(fig, x=tuple(xvals)*2, y=tuple(100*peptiderow[ xvals ])+tuple(100*(1-peptiderow[ xvals ])), specs=specs, lineopt=lineopt, markeropt=markeropt )
    
    self.__setNextColorInd()
    
    return # fig
  
  def __savePlot(self, saveFigOpts, fig, filename):
    """
    Save graph object GO as static png, dynamic html, etc

    Args:
        saveFigOpts (dict): _description_
        fig (_type_): _description_
        filename (_type_): _description_
    """
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)

    # # import os
    # folder = os.path.join( saveFigOpts['folder'], 'pngs' )
    # # save static png from plotly GO, requires plotly GO needs kaleido installed
    # if not os.path.exists(folder): os.makedirs(folder)
    # filepath = os.path.join(folder,filename+'.png')
    # fig.write_image( filepath )
    
    # save plotly GO as interactive graph in html
    # 
    # <script src="./js/plotly-1.58.5.min.js" charset="utf-8"></script>
    # <script src="./js/plotly-2.20.0.min.js" charset="utf-8"></script>
    # <script src="https://cdn.plot.ly/plotly-2.20.0.min.js" charset="utf-8"></script>
    # <script src="https://cdn.plot.ly/plotly-latest.min.js" charset="utf-8"></script>
    folder = os.path.join( saveFigOpts['folder'] )
    # folder = os.path.join( saveFigOpts['folder'], 'htmls' )
    if not os.path.exists(folder): os.makedirs(folder)
    filepath = os.path.join(folder, filename+".html")
    incPlotlyJs = "plotly.min.js"  # False # add cdn plotly js on the html directly to minimize size of each html # symbolic link
    fig.write_html( filepath, include_plotlyjs=incPlotlyJs )
    return
      
  # def abundancePlotProteinLevel(self, df, prtnGrp, labels=dict(), saveFigOpts = dict(), legendOpts=dict(), yerroropt=dict() ):
  def abundancePlotProteinLevel(self, xrange, xs, t12, prtnGrp, labels=dict(), saveFigOpts = dict(), legendOpts=dict() ):
    """
    Args:
        df (Dataframe): Pandas pivot table df
        trendlines (dict, optional): show, solid (bool), color (str), width (float) 
        lines (dict, optional): show, solid (bool), color (str), width (float) 
        markers (dict, optional): show (bool), symbol (str), size (float) 
    return: plotly graph object GO
    """
    labels = self.__setArgLabels(labels=labels)
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    legendOpts = self.__setArglegendOpts(legendOpts=legendOpts)
    # yerroropt = self.__setArgYerroropt(yerroropt=yerroropt)
    # lines = self.__setArgLines(lines=lines)
    # trendlines = self.__setArgTrendlines(trendlines=trendlines)
    # markers = self.__setArgMarkers(markers=markers)
    ys = 100 * np.exp(-np.log(2)/t12*xs)
    
    fig = go.Figure() 
    
    colors = dict(heavy='rgba(199,10,165,.9)', light='rgba(56,233,99,.9)')
    symbol = 'circle'
    
    # graph light series
    # curve fit
    markeropt = dict(color=colors['light'], symbol=symbol, size=2)
    # specs = dict(mode='lines', name='Light (degradation)', showlegend=True, connectgaps=True, legendgroup = 'light')
    specs = dict(mode='lines', name='Light (degradation)', showlegend=True, connectgaps=True)
    fig = self.__add1goTrace(fig, x=xs, y=ys, specs=specs, markeropt=markeropt )
    # data points
    # markeropt = dict(color=colors['light'], symbol=symbol, size=4)
    # specs = dict(mode='markers', name='Light (degradation)', showlegend=False, connectgaps=False, legendgroup = 'light')
    # xs = ['0'] + self.__xvalues 
    # ys = 100*df[[0]+self.__xvalues].mean()
    # fig = self.__add1goTrace(fig, x=xs, y=ys, specs=specs, markeropt=markeropt)
    # graph heavy series
    # curve fit
    markeropt = dict(color=colors['heavy'], symbol=symbol, size=2)
    # specs = dict(mode='lines', name='Heavy (synthesis)', showlegend=True, connectgaps=True, legendgroup = 'heavy')
    specs = dict(mode='lines', name='Heavy (synthesis)', showlegend=True, connectgaps=True)
    fig = self.__add1goTrace(fig, x=xs, y=100-ys, specs=specs, markeropt=markeropt )
    # data points
    # markeropt = dict(color=colors['heavy'], symbol=symbol, size=4)
    # specs = dict(mode='markers', name='Heavy (synthesis)', showlegend=False, connectgaps=False, legendgroup = 'heavy')
    # fig = self.__add1goTrace(fig, x=xs, y=100-ys, specs=specs, markeropt=markeropt)    
    
    
    if t12 < xrange: fig.add_vline(x=t12, line_width=1, line_dash="dash", line_color="black", annotation_text="&nbsp;<b>t<sub>½</sub></b> = "+str(t12.__round__(2)), annotation_position='bottom right')
    
    fig.update_layout( 
      title={
        'text': labels['title'],
        'x': 0.45,
        'xanchor': 'center'
        }, 
      xaxis_title=labels['x'], 
      yaxis_title=labels['y'],
      legend_title="Protein type",
      legend_tracegroupgap=4,
      font=dict(
          family=labels['fontfamily'],
          size=labels['size'],
          color="Black" # "RebeccaPurple"
      ), 
      legend=legendOpts
    )
    
    proteingenename = prtnGrp[1] if type(prtnGrp[1]) == str else prtnGrp[0]
    
    if saveFigOpts['savePlot']: 
      options = saveFigOpts.copy() 
      options['folder'] += 'htmls/proteinLevel/'
      self.__savePlot(options, fig, "RelAbundance_Gene-"+ proteingenename)
    if saveFigOpts['showPlot']: fig.show()
    
    return fig
  
  
  

#%%
# file = os.path.join(os.getcwd(), "../data/iMN_Peptide_Dataset.xlsx") # assuming cwd is .../Visualization/src/ folder
# file = os.path.join(os.getcwd(), "../data/20230522_dSILAC_Turnover_LightRelativeAbundances.xlsx") # assuming cwd is .../Visualization/src/ folder
# file = os.path.join(os.getcwd(), "../data/06202023_FinalReport_dSILAC_iMN_MultiTimePoint.xlsx") # assuming cwd is .../Visualization/src/ folder
# new file 06202023_FinalReport_dSILAC_iMN_MultiTimePoint.xlsx

# pto = ProteinTurnover(datafiles= dict(raw="../data/06202023_FinalReport_dSILAC_iMN_MultiTimePoint.xlsx", peptides=None, proteins=None) )
# pto = ProteinTurnover(datafiles= dict(raw=None, peptides="../data/dfPeptides20230724.csv", proteins="../data/dfProteins20230724.csv") )
# pto = ProteinTurnover(datafiles= dict(raw="../data/data20230801/1_iMN_alldata_forwebsite.xlsx", peptides=None, proteins=None) )
# pto = ProteinTurnover(datafiles= dict(raw="../data/data20230801/2_iMN_4Fraction.xlsx", peptides=None, proteins=None) )
# pto = ProteinTurnover(datafiles= dict(raw="../data/data20230801/3_iCN_DIAfractionate.xlsx", peptides=None, proteins=None) )
# pto = ProteinTurnover(datafiles= dict(raw="../data/data20230801/4_iCN_DDA.xlsx", peptides=None, proteins=None) )
# pto = ProteinTurnover(datafiles= dict(raw="../data/20230802_Peptide4web.xlsx", peptides=None, proteins=None) )
pto = ProteinTurnover(datafiles= dict(raw=None, peptides="../data/dfPeptides20230808.csv", proteins="../data/dfProteins20230808.csv") )

#%%
# saveplot, showplot = False, True
# saveplot, showplot = True, False
# saveplot, showplot = True, True
saveplot, showplot = False, False
savePath = "../media/plots/"

# pto.abundancePlotPgAll( saveFigOpts = dict(savePlot=saveplot, showPlot=showplot, folder=savePath) )
# pto.proteinHalflifeChart()


# %%
# Note, the new dataset (20230522_dSILAC_Turnover_LightRelativeAbundances.xlsx) has these three duplicate combos
# Protein	Peptide	  							
# Q6ZTK2	GLQLEGELEELRQDR	
# Q6ZSR9	QSIAGSVSITSLSSR	
# Q6ZSR9  TQNNLESDYLAR	

# Details:
# Q6ZTK2	GLQLEGELEELRQDR	
#    idx	Protein	Peptide	        Gene	Description	                                D0_L	D1_L	    D2_L	    D4_L	    D6_L	    D0_H	D1_H	    D2_H	    D4_H	    D6_H
#  75418	Q6ZTK2	GLQLEGELEELRQDR	NaN	  Putative uncharacterized protein LOC400499	1	    0.858829	0.677253	0.510069	0.260913	0	    0.141171	0.322747	0.489931	0.739087
# 147109	Q6ZTK2	GLQLEGELEELRQDR	NaN	  Putative uncharacterized protein LOC400499	1	    0.830570	0.708910	0.443543	0.158337	0	    0.169430	0.291090	0.556457	0.841663

# Q6ZSR9	QSIAGSVSITSLSSR	
#    idx	Protein	Peptide	        Gene	Description	                      D0_L	D1_L	    D2_L	D4_L	    D6_L	    D0_H	D1_H	    D2_H	D4_H	    D6_H
#  75423	Q6ZSR9	QSIAGSVSITSLSSR	NaN	  Uncharacterized protein FLJ45252	1	    NaN     	NaN	  0.444746	0.139033	0	    NaN	      NaN	  0.555254	0.860967
# 148536	Q6ZSR9	QSIAGSVSITSLSSR	NaN	  Uncharacterized protein FLJ45252	1	    0.815267	NaN	  NaN	      NaN	      0	    0.184733	NaN	  NaN	      NaN

# Q6ZSR9	TQNNLESDYLAR	
#    idx	Protein	Peptide	      Gene	Description	                      D0_L	D1_L	    D2_L	    D4_L	    D6_L	    D0_H	D1_H	    D2_H	    D4_H	    D6_H
#  75425	Q6ZSR9	TQNNLESDYLAR	NaN	  Uncharacterized protein FLJ45252	1	    0.676101	0.545013	0.636126	0.099641	0	    0.323899	0.454987	0.363874	0.900359
# 148537	Q6ZSR9	TQNNLESDYLAR	NaN	  Uncharacterized protein FLJ45252	1	    0.607315	0.413771	NaN	      0.255406	0	    0.392685	0.586229	NaN	      0.744594

# 20230620 data
# Peptide df has these extras, compared to PgGeneDescSummary df
# [('Q15147', 'PLCB4'),
#  ('Q9NRZ5', 'AGPAT4'),
#  ('P16885', 'PLCG2'),
#  ('Q9P2J9', 'PDP2'),
#  ('Q96JD6', 'AKR1E2'),
#  ('O75038', 'PLCH2'),
#  ('Q4KWH8', 'PLCH1'),
#  ('Q15119', 'PDK2'),
#  ('Q92506', 'HSD17B8'),
#  ('Q7RTP6', 'MICAL3'),
#  ('P19174', 'PLCG1'),
#  ('Q8TDZ2', 'MICAL1'),
#  ('P51178', 'PLCD1'),
#  ('Q15120', 'PDK3'),
#  ('Q99943', 'AGPAT1'),
#  ('Q9Y2I7', 'PIKFYVE'),
#  ('Q01970', 'PLCB3'),
#  ('Q9NRZ7', 'AGPAT3'),
#  ('Q8N3E9', 'PLCD3'),
#  ('Q96QU6', 'ACCS'),
#  ('Q04446', 'GBE1'),
#  ('P61604', 'HSPE1'),
#  ('Q9NUQ2', 'AGPAT5'),
#  ('Q9NQ66', 'PLCB1'),
#  ('Q15118', 'PDK1'),
#  ('Q9C0C9', 'UBE2O'),
#  ('O14874', 'BCKDK'),
#  ('Q9P0J1', 'PDP1')]
# 
# Summary df has these four extras however:
# [('Q8N1Y9', nan), ('B2RBV5', nan), ('Q6ZTK2', nan), ('Q6ZSR9', nan)].

# 20230714
# for A0A0B4J2D5(protein) GATD3B(gene) PIGLCCIAPVLAAK(peptide) , somehow it has 0.841 r2, but look terrible. 