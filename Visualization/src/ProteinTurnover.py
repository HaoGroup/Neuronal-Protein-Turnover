#%%
import pandas as pd
import os
import re
import plotly.graph_objects as go
# To export static images from go, need pip install -U kaleido
import plotly.express as px
import numpy as np
import ExpoDecayFit as edf
# from scipy.stats import hmean # harmonic mean
# import matplotlib.pyplot as plt
# import seaborn as sns

# Create class ProteinTurnover to handle and organize the various parts of the study
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
  def __init__(self, filepath) -> None:
    """_summary_
    Args:
        filepath (str): the source file location. Use os to ensure right format
        x_values (list, optional): the time marks on x-axis. Defaults to [0,1,2,4,6].
        x_unit (str, optional): the unit for x-axis. Defaults to "day".
    """

    self.__yAxisName = 'Relative Abundance'
    self.__yUnit = '%'
    self.__xAxisName = 'Time'
    self.__xUnit = 'day'
    # Match with ExpoDecayFit module:
    self.__statsTypes = ('b','t12','r2') # keeping track of three results/statics from exponential fit
    self.__modelTypes = ("LnLM1", "LnLM2", "CFit") # trying three different models in ExpoDecayFit: LnLM1 (sklearn LM), LnLM2 (statsmodels LM), CFit (scipy curve_fit)

    self.__compoIndexPeptide = ['Protein', 'Peptide'] # basic composite index used in Peptide and related DFs
    self.__compoIndexGene = ['Protein', 'Gene_x'] # basic composite index used in other dfs.
    self.__maxNAcnt = 1 # Each series only has at most 4 data points at day = 1,2,4,6.  Only allow at most 1 missing to be plotted
    self.df_Peptides = None # initialize, cleaned and re-structured df for analysis ['PG.ProteinGroups', 'PG.Genes', 'PG.ProteinDescriptions', 'Peptide', 0, 1, 2, 4, 6, 'k_results', 'Protein_Turnover_from_k', 'residuals']
    # self.PgList = [] # initialize, the list of ProteinGroups in the dataset
    self.df_PgPivot = [] # initialize, for plotting
    self.__ingestData(filepath=filepath) # set df_Peptides
    self.__setPgPivot()
    self.__modifyDfPeptides()
    return
  
  def __ingestData(self, filepath):
    """
    Args:
        filepath (_type_): the source file location. Use os to ensure right format
        x_unit (str, optional): the unit for x-axis. Defaults to None.
    Returns: None, will set df_Peptides
    """
    df = pd.read_excel(filepath) if (filepath[-5:] == '.xlsx' or filepath[-4:] == '.xls') else pd.read_csv(filepath) if (filepath[-4:] == '.csv') else None
    
    # clean up column names, avoid dots in names
    df.columns =  [ col.replace('.','_') for col in df.columns.values ]
    # 20230620 data file no longer has day0 columns. Re-create here:
    df['iMN_Day0_Light_Relative_Abundance'] = 1
    # df['iMN_Day0_Heavy_Relative_Abundance'] = 0
    
    # Keep only the light series, drop heavy
    df.drop(list(df.filter(regex = 'Heavy_Relative_Abundance')), axis = 1, inplace = True)
    
    # import re  # rename columns as 0, 1, 2, 4, 6 for days
    df.rename(columns=lambda x: re.sub(r'iMN_Day(\d).*_Relative_Abundance$', r'\1', x), inplace = True)
    df['chart'] = 0 # keep track of whether the peptide has enough data points to make a plot. default to 0

    self.df_Peptides = df.set_index(self.__compoIndexPeptide)  # set composite index and save

    if not self.chkDfPeptidesIndexUnique(): print("Protein Group in df_Peptides not unique")
    
    return None
  
  def __setPgPivot(self):
    """
    set up pivot table from df_Peptides table for analysis and plots
    Returns: None
    """
    # df_PgPivot table is indexed by the time/day value starting at 0.
    xAxisName = self.__xAxisName
    # select only rel_Abundance columns
    df = self.df_Peptides.loc[:, self.df_Peptides.columns.str.contains('^\d$')] # only columns with single digit as colnames
    df = df.stack(dropna = False) # melt method does not keep the original index keys. df is now a pd series
    df.index.names = df.index.names[:-1] + [xAxisName] # name new index column
    df.name = self.__yAxisName # set name of the series (column)
    df = pd.DataFrame(df).reset_index() # from pd series to dataframe, and reset index to use pivot functions
    # df = pd.melt()
    df[xAxisName] = df[xAxisName].astype(int) # need these x-values be taken as numeric
    self.df_PgPivot = 100 * df.pivot_table(index=xAxisName, columns = self.__compoIndexPeptide, values = self.__yAxisName, aggfunc='mean') # There are a small number of prtnGrp/peptide combo that are not unique (2 or more rows in the original excel. Take average.)
    return
  
  def __modifyDfPeptides(self): # add fitting result columns after pivot table is set.
    newcols = [ s+'_'+m for s in self.__statsTypes for m in self.__modelTypes ]
    self.df_Peptides[newcols] = np.nan
    return
  
  def peptidesFromPrtnGrp(self, prtnGrp): 
    """
    Args:
        prtnGrp (str): Protein Group name
        Returns: tuple
    """
    try: df = self.df_PgPivot[[prtnGrp]] # filter only one prtnGrp, can have multiple peptides, 
    except KeyError: return print(f"Protein Group {prtnGrp} is not found in dataset.") or ()
    except: return print(f"Unknown error when finding peptides from Protein Group {prtnGrp}.") or ()

    df.columns = df.columns.droplevel(0) # remove column index ProteinGroup
    peptides = self.__peptidesFromPivottable(df)
    return tuple(peptides)
  
  def __peptidesFromPivottable(self, df):
    """
    Obtain the list of peptides in the selected Pivot Table (unique values)
    Args:
        df (DataFrame): the pivot table in question
    Returns:
        list: list of peptide names
    """
    peptides = []
    _ = [ peptides.append(peptuple) for peptuple in df.columns.values if peptuple not in peptides ]
    
    return peptides
  
  # def __setPgGeneDescSummary(self):
  #   """
  #   setting up the df_PgGeneDescSummary as lookup table self.df_PgGeneDescSummary as well as 
  #   setting up the list of ProteinGroup values as attribute self.PgList
  #   It is also saving the main results of the protein charts.
    
  #   Returns: None
  #   """
  #   # set df_PgGeneDescSummary
  #   # compoIndexGene = ['PG.ProteinGroups', 'PG.Genes']
  #   # selectedCols = compoIndexGene + ['PG.ProteinDescriptions']
  #   compoIndexGene = self.__compoIndexGene
  #   selectedCols = compoIndexGene + ['Protein_Description']
  #   df = self.df_Peptides.reset_index()[selectedCols]
  #   df.set_index(compoIndexGene, inplace=True)
  #   # df['chart'] = 0 # keep track of whether the protein has enough data points to make a plot. default to 0
  #   # df.assign( chart = 0 ) # keep track of whether the protein has enough data points to make a plot. default to 0
  #   self.df_PgGeneDescSummary = df.drop_duplicates()
  #   if not self.chkDfPgGeneSummaryIndexUnique(): print("Protein Group - GeneId combo not unique")
    
  #   # next add list (tuple) of peptides to each row of prtnGrp
  #   # self.df_PgGeneDescSummary.loc[:,"peptides"] = self.df_PgGeneDescSummary.apply(lambda x: self.peptidesFromPrtnGrp(x.name[0]), axis = 1)
  #   # Warning: A value is trying to be set on a copy of a slice from a DataFrame. 
  #   # False positive: https://stackoverflow.com/questions/42105859/pandas-map-to-a-new-column-settingwithcopywarning.
  #   self.df_PgGeneDescSummary.assign( peptides = self.df_PgGeneDescSummary.apply(lambda x: self.peptidesFromPrtnGrp(x.name[0]), axis = 1) )
  #   return
  
  # def exportPgGeneDescSummary(self, filename="", format='json'):
  #   """
  #   Export df_PgGeneDescSummary table for web dev
  #   Args:
  #       format (str, optional): _description_. Defaults to 'json'.
  #   """
  #   filename = filename if filename else "PgGeneDescSummary"
  #   if format=='csv': 
  #     self.df_PgGeneDescSummary.to_csv(filename+".csv")
  #   else: # default json, will have protein as key, values will be gene_x, description, and peptide list
  #     df = self.df_PgGeneDescSummary.copy()
  #     df.reset_index(inplace=True)
  #     df.rename({'Protein': 'prtn', 'Gene_x': 'gene', 'Protein_Description': 'desc'}, axis=1, inplace=True)
  #     # df.set_index('Protein', inplace=True)
  #     df.to_json(filename+".json", orient="records")
  #   return
  
  def exportResults(self, filename="", format='json'):
    """
    Export df_Peptides table for web dev
    Args:
        format (str, optional): _description_. Defaults to 'json'.
    """
    filename = filename if filename else "results"
    df = self.df_Peptides.drop(list(self.df_Peptides.filter(regex = '^\d$')), axis = 1) # do not keep the time data, just prtn, gene, desc, bs, t12s, r2s.
    if format=='csv': 
      df.to_csv(filename+".csv")
      return
    
    # default json, will have protein as key, values will be gene_x, description, decay constants, half-lives, rsquared, and peptide list
    pepCol = 'peptides' # new column name for all peptide info per gene, for web json
    newHeaders = {'Protein': 'prtn', 'Gene_x': 'gene', 'Protein_Description': 'desc'}
    # [ v for i,v in newd.items() ]
    statsHeaders = { s: [ s+'_'+m for m in self.__modelTypes ] for s in self.__statsTypes } # the stats metrics that we look for, in particular, b-decay constant, and t12-halfLife #  'b_LnLM1', 'b_LnLM2', 'b_CFit' #  't12_LnLM1', 't12_LnLM2', 't12_CFit' # 'r2_LnLM1', 'r2_LnLM2', 'r2_CFit'
    modelMetricHeaders = sum( statsHeaders.values(), [])  #  'b_LnLM1', 'b_LnLM2', ... 't12_LnLM1' ... 'r2_CFit'
    colHeads = [ self.__compoIndexPeptide[1] ]+modelMetricHeaders # add 'Peptide' column (in df_Peptides) to the modelMetricHeaders

    df.reset_index(inplace=True)
    df.rename(newHeaders, axis=1, inplace=True)
    df_gb = df.groupby( list(newHeaders.values()), as_index=True, dropna=False)  # groupby object
    df_res = df_gb[statsHeaders['b']].agg('mean')  # first find mean b-values decay constants
    df_res[statsHeaders['t12']] = np.log(2)/df_res[statsHeaders['b']]  # next find half lives from bs. Essentially, we are finding the harmonic mean of half lives.
    df_res['chart'] = df_gb['chart'].agg('sum') # find total number of peptides have charts
    # next, get list of peptides and their metrics from different models into a dict()
    df_res[pepCol] = df_gb[colHeads].apply(lambda g: { h: tuple(g[h]) for h in colHeads} ) # adding the resulting pd series to df_res
    print(f'resulting results has dim: {df_res.shape}')
    print(df_res.head(2))
    df_res.reset_index(inplace=True)
    df_res.to_json(filename+".json", orient="records")
    return
  
  def chkDfPeptidesIndexUnique(self): 
    # print(f'df_Peptides composite index ({", ".join( self.df_Peptides.index.names )}) is unique? {self.df_Peptides.index.is_unique}') 
    return self.df_Peptides.index.is_unique
  
  # def chkDfPgGeneSummaryIndexUnique(self): 
  #   # print(f'df_PgGeneDescSummary composite index ({", ".join( self.df_PgGeneDescSummary.index.names )}) is unique? {self.df_PgGeneDescSummary.index.is_unique}')
  #   return self.df_PgGeneDescSummary.index.is_unique
  
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
  
  def __setArgLegendopt(self, legendopt=dict(title_font_family=None, font=dict(size=None)) ):
    """
    setting generic keyword argument legend with title_font_family, size
    Args:
        legendopt (dict, optional): legendgroup (str), title_font_family (str), font (dict), ...
    Returns:
        dict: { title_font_family, font, ... }
    """
    res = legendopt.copy()
    res['title_font_family'] = res['title_font_family'] if (res.__contains__('title_font_family') and res['title_font_family']) else 'Garamond'
    res['font'] = res['font'] if (res.__contains__('font') and res['font']) else dict(size=8)
    res['font']['size'] = res['font']['size'] if (res['font'].__contains__('size') and res['font']['size']>0) else 7
    return res
  
  def __setArgYerroropt(self, yerroropt=dict(type=None, array=None, arrayminus=None, symmetric=None, visible=None) ):
    """
    setting generic keyword argument yerror with type, array, arrayminus, symmetric, visible, value
    Args:
        yerroropt (dict, optional): type (str), array (list), arrayminus (list), symmetric (bool), visible (bool), value (float, for percent type)...
    Returns:
        dict: { type, array, ... }
    """
    # fig = go.Figure(data=go.Scatter(
    #     x=[0, 1, 2],
    #     y=[6, 10, 2],
    #     error_y=dict(
    #         type='data', # value of error bar given in data coordinates
    #         array=[1, 2, 3],
    #         visible=True)
    # ))
    res = yerroropt.copy()
    res['type'] = res['type'] if (res.__contains__('type') and res['type']) else 'data' # 'data', 'percent'
    res['symmetric'] = res['symmetric'] if ( res.__contains__('symmetric') and isinstance(res['symmetric'], bool) ) else True
    res['visible'] = res['visible'] if ( res.__contains__('visible') and isinstance(res['visible'], bool) ) else True
    res['array'] = res['array'] if (res.__contains__('array') and (not res['array'] is None) ) else None
    res['arrayminus'] = res['arrayminus'] if (res.__contains__('arrayminus') and (not res['arrayminus'] is None) ) else None
    res['value'] = res['value'] if (res.__contains__('value') and (not res['value'] is None) and res['type']=='percent') else None
    return res
  
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
  
  def __add1goTrace(self, fig, x, y, specs=dict(), lineopt=dict(), markeropt=dict(), legendopt=dict(), yerroropt=dict() ):
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
    yerroropt = self.__setArgYerroropt(yerroropt=yerroropt)
    # check lineopt standards
    # check markeropt standards
    # https://plotly.com/python/marker-style/ # hourglass, cross, x-thin, ...
    mode = specs['mode'] if (specs.__contains__('mode') and specs['mode']) else 'lines'
    name = specs['name'] if (specs.__contains__('name') and specs['name']) else '-'
    connectgaps = specs['connectgaps'] if (specs.__contains__('connectgaps') and isinstance(specs['connectgaps'], bool) ) else True

    showlegend = specs['showlegend'] if ( specs.__contains__('showlegend') and isinstance(specs['showlegend'], bool) ) else False
    legendgroup = specs['legendgroup'] if (specs.__contains__('legendgroup') and specs['legendgroup']) else specs['name']
    
    if ( (not yerroropt['array'] is None) and yerroropt['visible']):
      fig.add_trace(go.Scatter( mode=mode, x=x, y=y, error_y=yerroropt, showlegend=showlegend, legendgroup=legendgroup, name=name, connectgaps=connectgaps, line=lineopt, marker=markeropt))
    else: # add more conditions/scenarios here if needed
      fig.add_trace(go.Scatter( mode=mode, x=x, y=y,showlegend=showlegend, legendgroup=legendgroup, name=name, connectgaps=connectgaps, line=lineopt, marker=markeropt))
      
    # if showlegend: fig.update_layout( showlegend=showlegend, legend=legendopt )
    return fig
  
  def abundancePlot1Peptide(self, fig, df, prtnGrp, peptide, lineopt=dict(), markeropt=dict(), legendopt=dict()):
    """
    Args:
        fig (plotly fig): to be added with this new line plot
        df (Dataframe): Pandas pivot table df
        peptide (str): Peptide name
        lines (dict, optional): show, solid (bool), color (str), width (float) 
        markers (dict, optional): show (bool), symbol (str), size (float) 
    return: plotly plot object, colorcntadd (0 or 1)
    """
    lineopt = self.__setArgLineopt(lineopt=lineopt)
    markeropt = self.__setArgMarkeropt(markeropt=markeropt)
    markeropt['color'] = lineopt['color']
    legendopt = self.__setArgLegendopt(legendopt=legendopt)
    #
    # df.columns = df.columns.droplevel(0) # further drop peptide level header 
    df.reset_index(inplace=True) # reset index to have 'time' column, for plotting
    
    # now with ExpoDecayFit module
    model = edf.ExpoDecayFit(df) # model.startx, starty, samplexs, sampleys
    bs, t12s, r2s = [ model.modelsummary.loc[t,:] for t in self.__statsTypes ]
    stats = dict( b=bs, t12=t12s, r2=r2s )

    modelChoice = 'CFit' # pick one to plot here. ("LnLM1", "LnLM2", "CFit")
    ysamples = model.sampleys[modelChoice]
    # xsamples = model.samplexs
    # xrange = round( min( max(6, 1.8*t12s[modelChoice]) , 10) ) # no more than 10, between 6 and 10. If t12 is close, show 1.8*t12
    
    # model lines
    specs = dict(name=peptide+f' t = {stats["t12"][modelChoice].__round__(1)}d', connectgaps=False, mode='lines',  showlegend=False, legendgroup = peptide)
    markeropt['size'] = 2
    # legendopt['showlegend'] = True 
    specs['showlegend'] = True 
    fig = self.__add1goTrace(fig, x=tuple(model.samplexs), y=tuple(ysamples), specs=specs, lineopt=lineopt, markeropt=markeropt )
    # fig = self.__add1goTrace(fig, x=tuple(model.samplexs)+(np.nan,)+tuple(model.samplexs ), y=tuple(ysamples)+(np.nan,)+tuple(100-ysamples) , specs=specs, lineopt=lineopt, markeropt=markeropt )

    specs = dict(connectgaps=False, mode='lines', showlegend=False, legendgroup = peptide)
    fig = self.__add1goTrace(fig, x=tuple(model.samplexs), y=tuple(100-ysamples), specs=specs, lineopt=lineopt, markeropt=markeropt )
    

    # data points
    specs = dict(name=peptide+f' t = {stats["t12"][modelChoice].__round__(1)}d', connectgaps=False, mode='markers',  showlegend=False, legendgroup = peptide)
    markeropt['size'] = 4
    fig = self.__add1goTrace(fig, x=tuple(df[self.__xAxisName])*2, y=tuple(df[peptide])+tuple(100-df[peptide]), specs=specs, lineopt=lineopt, markeropt=markeropt )
    
    # combine both light and heavy series data points and model lines
    # xall = tuple(df[self.__xAxisName])*2 +tuple(model.samplexs)*2
    # yall = tuple(df[peptide])+tuple(100-df[peptide])+tuple(ysamples)+tuple(100-ysamples)

    
    # fig = self.__add1goTrace(fig, x=xall, y=yall, specs=specs, lineopt=lineopt, markeropt=markeropt )
    
    # update chart value for this peptide
    if len(fig.data) >0 : self.df_Peptides.loc[ (prtnGrp[0],peptide) , 'chart' ] = 1
    # save results in df_Peptides
    for m in self.__modelTypes:
      for s in self.__statsTypes:
        self.df_Peptides.loc[( prtnGrp[0], peptide ), s+'_'+m] = stats[s][m]
    
    return fig
  
  def abundancePlotAllPeptides(self, df, prtnGrp, labels=dict(), saveFigOpts=dict(), legendopt=dict() ):
    """
    Args:
        df (Dataframe): Pandas pivot table df
        prtnGrp (tuple): Protein-Gene index
        labels (dict, optional): x-, y-labels and title. Defaults to empty dictionary
        lines (dict, optional): show, solid (bool), color (str), width (float) 
        markers (dict, optional): show (bool), symbol (str), size (float) 
    return: plotly plot object
    """    
    labels = self.__setArgLabels(labels=labels)
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    legendopt = self.__setArgLegendopt(legendopt=legendopt)

    peptides = self.__peptidesFromPivottable(df)

    cpalette = px.colors.qualitative.Dark24  # trendline use '#00FE35', which is Light24[1] color, lime green
    colorcnt = 0 # initialize # can consider making this global if a random start is preferred.
    colorcntmax = len(cpalette)    
    fig = go.Figure()
    for peptide in peptides: 
      fig = self.abundancePlot1Peptide(fig=fig, df=df[[peptide]], prtnGrp=prtnGrp, peptide=peptide, lineopt = dict( color=cpalette[ colorcnt%colorcntmax ], width=2), legendopt=legendopt )
      colorcnt += 1
    
    if len(fig.data) < 1 : return #  if nothing showing, skip
    
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
      legend=legendopt
    )
    
    proteingenename = prtnGrp[1] if type(prtnGrp[1]) == str else prtnGrp[0]
    
    if saveFigOpts['savePlot']:
      options = saveFigOpts.copy() 
      options['folder'] += 'htmls/peptides/'
      self.__savePlot(options, fig, "RelAbundance_Gene-"+ proteingenename +"-peptides")
    if saveFigOpts['showPlot']: fig.show()
    
    return fig
  
  def __savePlot(self, saveFigOpts, fig, filename):
    """
    Save plot object as static png, dynamic html, etc

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
      
  def abundancePlotProteinLevel(self, df, prtnGrp, labels=dict(), saveFigOpts = dict(), legendopt=dict(), yerroropt=dict() ):
    """
    Args:
        df (Dataframe): Pandas pivot table df
        trendlines (dict, optional): show, solid (bool), color (str), width (float) 
        lines (dict, optional): show, solid (bool), color (str), width (float) 
        markers (dict, optional): show (bool), symbol (str), size (float) 
    return: plotly plot objecta
    """
    labels = self.__setArgLabels(labels=labels)
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    legendopt = self.__setArgLegendopt(legendopt=legendopt)
    yerroropt = self.__setArgYerroropt(yerroropt=yerroropt)
    dataColName = 'ProteinAverage' # use this when averging multiple peptides
    # lines = self.__setArgLines(lines=lines)
    # trendlines = self.__setArgTrendlines(trendlines=trendlines)
    # markers = self.__setArgMarkers(markers=markers)
    # 
    # Prep the data
    if (len(df.columns) > 1) : # Pg_level needing to average all peptides, also determine error bars
      stdevs = df.std(axis=1) # these are in percentages, 0-100
      confMultiplier = 1.0 # For 1 sd, 68% confidence if random sample. But Proteins are not consisting of random samples of peptide?
      yerroropt['array'] = confMultiplier * stdevs if not np.allclose(stdevs,0, atol=1E-1) else None
      df_avg = pd.DataFrame( { dataColName : df.mean(axis = 1) } ) # na will be ignored !! Re-think strategy
    else:
      dataColName = df.columns.values[-1] # only 1 peptide, use it as dataColName
      yerroropt['array'] = None
      df_avg = df.copy() # peptide level data
    df_avg.reset_index(inplace=True) # reset index to have 'time' column, for plotting
    xs = df_avg[self.__xAxisName]
    ys = dict( light=df_avg[dataColName] , heavy=100-df_avg[dataColName])
    
    cpalette = px.colors.qualitative.Dark24
    
    # call ExpoDecayFit module
    fit = edf.ExpoDecayFit(df_avg, xAxisName = self.__xAxisName, modelTypes=self.__modelTypes, statsTypes=self.__statsTypes) # model.startx, starty, samplexs, sampleys
    
    bs, t12s, r2s = [ fit.modelsummary.loc[t,:] for t in self.__statsTypes ]
    
    xsamples = fit.samplexs
    ysamples = dict( light=fit.sampleys, heavy=100-fit.sampleys )
    
    # calculate harmonic mean half-life for protein here
    bavg = sum(bs)/len(bs)
    t12avg = np.log(2)/bavg
    
    modelChoice = 'CFit'
    xrange = round( min( max(6, 1.8*t12avg) , 10) ) # no more than 10, between 6 and 10. If t12 is close, show 1.8*t12
 
    # colors = dict(heavy='rgba(199,10,165,.9)', light='rgba(56,233,99,.9)')
    # symbols = dict(heavy='x', light='cross') # try hexagon
    # symbols = dict(heavy='circle', light='circle') # try hexagon
    symbol = 'circle'
    types = ('light','heavy')
    models = self.__modelTypes # (self.__modelTypes[0], self.__modelTypes[2]) # 'LnLM1', 'LnLM2', 'CFit'
    colorCnt, incr = 0, 2
    legendnames = dict( light='Light data (degradation)', heavy='Heavy data (synthesis)')
    # legendopt['showlegend'] =  True 
    
    fig = go.Figure()
    for t in types:
      # data
      markeropt = dict(color=cpalette[colorCnt], symbol=symbol, size=4)
      specs = dict(mode='markers', name=legendnames[t], showlegend=True, connectgaps=False)
      fig = self.__add1goTrace(fig, x=xs, y=ys[t], specs=specs, markeropt=markeropt, yerroropt=yerroropt )
      colorCnt += incr
      
      # curve model fit
      for m in models:
        markeropt = dict(color=cpalette[colorCnt], symbol=symbol, size=2)
        specs = dict(mode='lines', name=f'{t.capitalize()} ({m}: {t12s[m].__round__(1)}d, {r2s[m].__round__(3)})', showlegend=True, connectgaps=True)
        fig = self.__add1goTrace(fig, x=xsamples, y=ysamples[t][m], specs=specs, markeropt=markeropt )
        colorCnt += incr
        
    if len(fig.data) < 1 : return #  if nothing showing, skip
    
    # show half life if within range
    if t12avg < xrange: fig.add_vline(x=t12avg, line_width=1, line_dash="dash", line_color="black", annotation_text="&nbsp;<b>t<sub>½</sub></b> = "+str(t12avg.__round__(2)), annotation_position='bottom right')
    
    fig.update_layout( 
      title={
        'text': labels['title'],
        'x': 0.45,
        'xanchor': 'center'
        }, 
      xaxis_title=labels['x'], 
      yaxis_title=labels['y'],
      legend_title="Data vs Model",
      legend_tracegroupgap=4,
      font=dict(
          family=labels['fontfamily'],
          size=labels['size'],
          color="Black" # "RebeccaPurple"
      ), 
      legend=legendopt
    )
    
    proteingenename = prtnGrp[1] if type(prtnGrp[1]) == str else prtnGrp[0]
    
    if saveFigOpts['savePlot']: 
      options = saveFigOpts.copy() 
      options['folder'] += 'htmls/proteins/'
      self.__savePlot(options, fig, "RelAbundance_Gene-"+ proteingenename)
    if saveFigOpts['showPlot']: fig.show()
    
    return fig
  
  def abundancePlot1Pg(self, prtnGrp, labels=dict(), saveFigOpts = dict() ):
    """
    Args:
        prtnGrp (tuple): the (ProteinGrp, Gene) being processed
        labels (dict, optional): x-, y-labels and title. Defaults to empty dictionary
        saveFigOpts (dict, optional): savePlot (binary) and folder (str). Defaults to dict(savePlot=False, folder=None).
    return: None
    """
    labels = self.__setArgLabels(labels=labels)
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    
    df = self.df_PgPivot[[prtnGrp[0]]] # filter only one prtnGrp, can have multiple peptides, 
    df.columns = df.columns.droplevel(0) # remove column index ProteinGroup
    
    # REMOVE peptides with more than __maxNAcnt na values
    df = df.drop(columns=[col for col in df if df[col].isna().sum() > self.__maxNAcnt ]) # __maxNAcnt is set globally
    if df.shape[1] < 1: return None # nothing to plot
    
    labels['title'] = 'Gene: ' + prtnGrp[1] if type(prtnGrp[1]) == str else 'Protein Group: '+prtnGrp[0]
    _ = self.abundancePlotAllPeptides(df, prtnGrp=prtnGrp, labels=labels, saveFigOpts=saveFigOpts)
    
    # Now plot protein level average with trendline
    labels['title'] = "Average for "+ labels['title']
    _ = self.abundancePlotProteinLevel(df, prtnGrp=prtnGrp, labels=labels, saveFigOpts=saveFigOpts)

    return 
  
  
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
    plotmax = 12 # in/out
    # for prtnGrp in self.PgList:
    #
    # Need to get list of (protein,gene) 
    # This method somehow is dropping the ones with NA description. Don't know why.
    # compoIndexGene = pto._ProteinTurnover__compoIndexGene
    # dfp = pto.df_Peptides.reset_index()[compoIndexGene+['Protein_Description']].set_index(compoIndexGene).drop_duplicates()
    #
    # Use this groupby method (results in 28 more rows, 9883 vs 9855)
    # 
    df = self.df_Peptides.reset_index().loc[:, pto._ProteinTurnover__compoIndexGene + ['chart']].groupby( self.__compoIndexGene , as_index=True, dropna=False)[['chart']].agg('sum')
    #
    for prtnGrp in df.index.values: # multi-level index with (protein, gene)
      if plotmax == 0 : break # in/out
      self.abundancePlot1Pg(prtnGrp=prtnGrp, labels=labels, saveFigOpts=saveFigOpts)
      plotmax -= 1 # in/out
      
    return
  

#%%
# file = os.path.join(os.getcwd(), "../data/iMN_Peptide_Dataset.xlsx") # assuming cwd is .../Visualization/src/ folder
# file = os.path.join(os.getcwd(), "../data/20230522_dSILAC_Turnover_LightRelativeAbundances.xlsx") # assuming cwd is .../Visualization/src/ folder
file = os.path.join(os.getcwd(), "../data/06202023_FinalReport_dSILAC_iMN_MultiTimePoint.xlsx") # assuming cwd is .../Visualization/src/ folder
# new file 06202023_FinalReport_dSILAC_iMN_MultiTimePoint.xlsx

pto = ProteinTurnover(filepath= file)
# pto.chkDfPgIndexUnique()
#%%
# saveplot, showplot = False, True
saveplot, showplot = True, False
# saveplot, showplot = True, True
# saveplot, showplot = False, False
savePath = "../media/plots/"

# pto.abundancePlotPgAll(savePlot=saveplot, saveFolder=savePath)
pto.abundancePlotPgAll( saveFigOpts = dict(savePlot=saveplot, showPlot=showplot, folder=savePath) )


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
