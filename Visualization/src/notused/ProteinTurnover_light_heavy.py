#%%
import pandas as pd
import os
import re
import plotly.graph_objects as go
# To export static images from go, need pip install -U kaleido
import plotly.express as px
import numpy as np
import ExpoDecayFit as edf
# from ExpoDecayFit import ExpoDecayFit # in the same folder
# from scipy.stats import hmean # harmonic mean
# from sklearn.linear_model import LinearRegression
# import matplotlib.pyplot as plt
# import seaborn as sns
# import streamlit as st

#%%
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
    # self.__compoIndex = ['PG.ProteinGroups', 'Peptide', 'wtype'] # basic composite index used in DFs
    # Match with ExpoDecayFit module:
    self.__modelTypes = ("LnLM1", "LnLM2", "CFit") # trying three different models in ExpoDecayFit: LnLM1 (sklearn LM), LnLM2 (statsmodels LM), CFit (scipy curve_fit)
    self.__statsTypes = ('b','t12','r2') # keeping track of three results/statics from exponential fit

    self.__compoIndex = ['Protein', 'Peptide', 'wtype'] # basic composite index used in DFs
    self.__maxNAcnt = 1 # Each series only has at most 4 data points at day = 1,2,4,6.  Only allow at most 1 missing to be plotted
    self.df_PgRaw = None # initialize, current structure has columns: ['PG.ProteinGroups', 'PG.Genes', 'PG.ProteinDescriptions', 'Peptide', 'iMN_Day{x}_{Light/Heavy}_Relative_Abundance', 'k_results', 'Protein_Turnover_from_k', 'residuals']
    self.df_Peptides = None # initialize, cleaned and re-structured df for analysis ['PG.ProteinGroups', 'PG.Genes', 'PG.ProteinDescriptions', 'Peptide', 'wtype', 0, 1, 2, 4, 6, 'k_results', 'Protein_Turnover_from_k', 'residuals']
    self.df_PgGeneDescSummary = None # initialize, serve as lookup table between PgProteinGroups, PgGenes, and PgProteinDescriptions, PLUS list/tuple of peptides
    # self.PgList = [] # initialize, the list of ProteinGroups in the dataset
    self.df_PgPivot = [] # initialize, for plotting
    self.__ingestData(filepath=filepath) # set df_PgRaw and df_Peptides
    self.__setPgPivot()
    self.__modifyDfPeptides()
    self.__setPgGeneDescSummary() # set df_PgGeneDescSummary as well as PgList, but need df_PgPivot first
    return
  
  def __ingestData(self, filepath):
    """
    Args:
        filepath (_type_): the source file location. Use os to ensure right format
        x_unit (str, optional): the unit for x-axis. Defaults to None.

    Returns: None, will set df_PgRaw and df_Peptides
    """
    self.df_PgRaw = self.df_PgRaw = pd.read_excel(filepath) if (filepath[-5:] == '.xlsx' or filepath[-4:] == '.xls') else pd.read_csv(filepath) if (filepath[-4:] == '.csv') else None
    
    # clean up column names, avoid dots in names
    self.df_PgRaw.columns =  [ col.replace('.','_') for col in self.df_PgRaw.columns.values ]
    # 20230620 data file no longer has day0 columns. Re-create here:
    self.df_PgRaw['iMN_Day0_Light_Relative_Abundance'] = 1
    self.df_PgRaw['iMN_Day0_Heavy_Relative_Abundance'] = 0
    
    # separate out data for light and heavy types, append them instead, with day_x as the same column
    df_light = self.df_PgRaw.copy()
    df_heavy = df_light.copy()
    # set new column to keep track of light/heavy
    df_light["wtype"] = "light"
    df_heavy["wtype"] = "heavy"
    # keep only light/heavy data for each df
    df_light.drop(list(df_light.filter(regex = 'Heavy_Relative_Abundance')), axis = 1, inplace = True)
    df_heavy.drop(list(df_heavy.filter(regex = 'Light_Relative_Abundance')), axis = 1, inplace = True)
    # df_light.drop(list(df_light.filter(regex = 'Relative Heavy Intensity')), axis = 1, inplace = True) # Not used after renamed new xlxs column names
    # df_heavy.drop(list(df_heavy.filter(regex = 'Relative Light Intensity')), axis = 1, inplace = True) # Not used after renamed new xlxs column names
    # rename columns as 0, 1, 2, 4, 6 for days
    # import re
    df_light.rename(columns=lambda x: re.sub(r'iMN_Day(\d).*_Relative_Abundance$', r'\1', x), inplace = True)
    df_heavy.rename(columns=lambda x: re.sub(r'iMN_Day(\d).*_Relative_Abundance$', r'\1', x), inplace = True)
    # df_light.rename(columns=lambda x: re.sub(r'Relative.*Intensity.Day.(\d)$', r'\1', x), inplace = True) # Not used after renamed new xlxs column names
    # df_heavy.rename(columns=lambda x: re.sub(r'Relative.*Intensity.Day.(\d)$', r'\1', x), inplace = True) # Not used after renamed new xlxs column names
    # set composite index 
    df_light.set_index(self.__compoIndex, inplace=True)  
    df_heavy.set_index(self.__compoIndex, inplace=True)  
    # stack the two df together
    self.df_Peptides = pd.concat([df_light, df_heavy])
    # self.chkDfPgIndexUnique()
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
    self.df_PgPivot = 100 * df.pivot_table(index=xAxisName, columns = self.__compoIndex, values = self.__yAxisName, aggfunc='mean') # There are a small number of prtnGrp/peptide combo that are not unique (2 or more rows in the original excel. Take average.)
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
    try: df = self.df_PgPivot[[prtnGrp]] # filter only one prtnGrp, with both light and heavy data, can have multiple peptides, 
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
    _ = [ peptides.append(peptuple[0]) for peptuple in df.columns.values if peptuple[0] not in peptides ]
    
    return peptides
  
  def __setPgGeneDescSummary(self):
    """
    setting up the df_PgGeneDescSummary as lookup table self.df_PgGeneDescSummary as well as 
    setting up the list of ProteinGroup values as attribute self.PgList
    It is also saving the main results of the protein charts.
    
    Returns: None
    """
    # set PgList
    # self.PgList = self.df_PgRaw["PG.ProteinGroups"].unique()
    # self.PgList = self.df_PgRaw["Protein"].unique()
    # set df_PgGeneDescSummary
    # compoIndex = ['PG.ProteinGroups', 'PG.Genes']
    compoIndex = ['Protein', 'Gene_x']
    # selectedCols = compoIndex + ['PG.ProteinDescriptions']
    selectedCols = compoIndex + ['Protein_Description']
    df = self.df_PgRaw[selectedCols]
    df.set_index(compoIndex, inplace=True)
    df['chart'] = 0 # keep track of whether the protein has enough data points to make a plot. default to 0
    self.df_PgGeneDescSummary = df.drop_duplicates()
    if not self.chkDfPgGeneSummaryIndexUnique(): print("Protein Group - GeneId combo not unique")
    
    # next add list (tuple) of peptides to each row of prtnGrp
    self.df_PgGeneDescSummary.loc[:,"peptides"] = self.df_PgGeneDescSummary.apply(lambda x: self.peptidesFromPrtnGrp(x.name[0]), axis = 1)
    # Warning: A value is trying to be set on a copy of a slice from a DataFrame. 
    # False positive: https://stackoverflow.com/questions/42105859/pandas-map-to-a-new-column-settingwithcopywarning.
    return
  
  def exportPgGeneDescSummary(self, filename="", format='json'):
    """
    Export df_PgGeneDescSummary table for web dev
    Args:
        format (str, optional): _description_. Defaults to 'json'.
    """
    filename = filename if filename else "PgGeneDescSummary"
    if format=='csv': 
      self.df_PgGeneDescSummary.to_csv(filename+".csv")
    else: # default json, will have protein as key, values will be gene_x, description, and peptide list
      df = self.df_PgGeneDescSummary.copy()
      df.reset_index(inplace=True)
      df.rename({'Protein': 'prtn', 'Gene_x': 'gene', 'Protein_Description': 'desc'}, axis=1, inplace=True)
      # df.set_index('Protein', inplace=True)
      df.to_json(filename+".json", orient="records")
    return
  
  def exportPgPeptides(self, filename="", format='json'):
    """
    Export df_Peptides table for web dev
    Args:
        format (str, optional): _description_. Defaults to 'json'.
    """
    filename = filename if filename else "PgPeptides"
    df = self.df_Peptides[self.df_Peptides['wtype']=='light']
    df.filter(~df.columns.isin(['0','1','2','4','6'])) # do not keep the time data, just prtn, gene, desc, bs, t12s, r2s.
    if format=='csv': 
      df.to_csv(filename+".csv")
    else: # default json, will have protein as key, values will be gene_x, description, and peptide list
      df.reset_index(inplace=True)
      df.rename({'Protein': 'prtn', 'Gene_x': 'gene', 'Protein_Description': 'desc'}, axis=1, inplace=True)
      # df.set_index('Protein', inplace=True)
      df.to_json(filename+".json", orient="records")
    return
  
  def chkDfPeptidesIndexUnique(self): 
    # print(f'df_Peptides composite index ({", ".join( self.df_Peptides.index.names )}) is unique? {self.df_Peptides.index.is_unique}') 
    return self.df_Peptides.index.is_unique
  
  def chkDfPgGeneSummaryIndexUnique(self): 
    # print(f'df_PgGeneDescSummary composite index ({", ".join( self.df_PgGeneDescSummary.index.names )}) is unique? {self.df_PgGeneDescSummary.index.is_unique}')
    return self.df_PgGeneDescSummary.index.is_unique
  
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
        legendopt (dict, optional): title_font_family (str), font (dict), ...
    Returns:
        dict: { title_font_family, font, ... }
    """
    res = legendopt.copy()
    res['title_font_family'] = res['title_font_family'] if (res.__contains__('title_font_family') and res['title_font_family']) else 'Garamond'
    res['font'] = res['font'] if (res.__contains__('font') and res['font']) else dict(size=8)
    res['font']['size'] = res['font']['size'] if (res['font'].__contains__('size') and res['font']['size']>0) else 7
    return res
  
  # legend=dict(
  #   title_font_family='Courier New',
  #   font=dict(
  #       size=8
  #   )
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
  
  def __setArgSaveFigOpts(self, saveFigOpts=dict(savePlot=False, folder=None) ):
    """
    setting generic keyword argument saveFigOpts into savePlot and folder
    Args:
        saveFigOpts (dict, optional): savePlot (binary) and folder (str). Defaults to dict(savePlot=False, folder=None).
    Returns:
        dict: { savePlot, folder }
    """
    res = saveFigOpts.copy()
    res['savePlot'] = res['savePlot'] if ( res.__contains__('savePlot') and isinstance(res['savePlot'], bool) ) else False
    res['folder'] = res['folder'] if ( res.__contains__('folder') and res['folder'] ) else './'
    return res
  
  def __add1goTrace(self, fig, x, y, specs=dict(), lineopt=dict(), markeropt=dict(), legendopt=dict() ):
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
    # https://plotly.com/python/marker-style/ # hourglass, cross, x-thin, ...
    mode = specs['mode'] if (specs.__contains__('mode') and specs['mode']) else 'lines'
    showlegend = specs['showlegend'] if ( specs.__contains__('showlegend') and isinstance(specs['showlegend'], bool) ) else False
    name = specs['name'] if (specs.__contains__('name') and specs['name']) else '-'
    connectgaps = specs['connectgaps'] if (specs.__contains__('connectgaps') and isinstance(specs['connectgaps'], bool) ) else True
    # check lineopt standards
    # check markeropt standards
    
    fig.add_trace(go.Scatter( mode=mode, x=x, y=y, showlegend=showlegend, name=name, connectgaps=connectgaps, line=lineopt, marker=markeropt ))
    # if showlegend: fig.update_layout( showlegend=showlegend, legend=legendopt )
    
    return fig
  
  def abundancePlot1Peptide(self, fig, df, prtnGrp, peptide, wtype='both', lineopt=dict(), markeropt=dict()):
    """
    Args:
        fig (plotly fig): to be added with this new line plot
        df (Dataframe): Pandas pivot table df
        peptide (str): Peptide name
        wtype (str, optional): weight_type: "heavy", "light", or "both". Defaults to 'both'.
        lines (dict, optional): show, solid (bool), color (str), width (float) 
        markers (dict, optional): show (bool), symbol (str), size (float) 
    return: plotly plot object, colorcntadd (0 or 1)
    """
    lineopt = self.__setArgLineopt(lineopt=lineopt)
    markeropt = self.__setArgMarkeropt(markeropt=markeropt)
    markeropt['color'] = lineopt['color']
    # legendopt = self.__setArgLegendopt(legendopt=legendopt)
    #
    df.columns = df.columns.droplevel(0) # further drop peptide level header 
    # if ( isinstance(df.columns, pd.core.indexes.multi.MultiIndex) and len(df.columns.levels) > 1 ) : df.columns = df.columns.droplevel(0) # was using this function for some other cases. No longer need this check.
    df.reset_index(inplace=True) # reset index to have 'time' column, for plotting
    # 
    # series = ('light') if wtype == 'light' else ('heavy') if wtype == 'heavy' else ('light','heavy')
    # line plots
    # for i, wtype in enumerate(series):
    #   showlegend = False if i>0 else True, # only show when i=0 # somehow this line is resulting showlegend as a tuple (True, )
    #   showlegend = showlegend[0] if isinstance(showlegend,tuple) else showlegend if isinstance(showlegend,bool) else True
    #   specs = dict(name=peptide, connectgaps=True, showlegend=showlegend, mode='lines+markers')
      
    #   fig = self.__add1goTrace(fig, x=df[self.__xAxisName], y=df.loc[:,wtype], specs=specs, lineopt=lineopt, markeropt=markeropt )
    
    # fitted plots
    # types=('light','heavy')
    # Plot scatter data points first
    showlegend = False
    specs = dict(name=peptide, connectgaps=False, showlegend=showlegend, mode='markers')
    fig = self.__add1goTrace(fig, x=df[self.__xAxisName], y=df.loc[:,'light'], specs=specs, lineopt=lineopt, markeropt=markeropt )
    fig = self.__add1goTrace(fig, x=df[self.__xAxisName], y=df.loc[:,'heavy'], specs=specs, lineopt=lineopt, markeropt=markeropt )
    
    # now with ExpoDecayFit module
    model = edf.ExpoDecayFit(df) # model.startx, starty, samplexs, sampleys
    bs = model.modelsummary.loc['b',:]
    t12s = model.modelsummary.loc['t12',:]
    r2s = model.modelsummary.loc['r2',:]
    stats = dict( b=bs, t12=t12s, r2=r2s )
    xsamples = model.samplexs
    ysamples = dict()
    ysamples['light'] = model.sampleys
    ysamples['heavy'] = 100 - ysamples['light']    
    modelChoice = 'CFit' # pick one to plot here. ("LnLM1", "LnLM2", "CFit")

    xrange = round( min( max(6, 1.8*t12s[modelChoice]) , 10) ) # no more than 10, between 6 and 10. If t12 is close, show 1.8*t12
    
    # heavy first, no legend
    specs = dict(name=peptide, connectgaps=True, showlegend=showlegend, mode='lines')
    fig = self.__add1goTrace(fig, x=xsamples, y=ysamples['heavy'][modelChoice], specs=specs, markeropt=markeropt )
    # now light series, with legend
    specs['showlegend'] = True 
    specs['name'] += f' t = {stats["t12"][modelChoice].__round__(3)}'
    fig = self.__add1goTrace(fig, x=xsamples, y=ysamples['light'][modelChoice], specs=specs, markeropt=markeropt )
    
    # save results in df_Peptides
    for m in self.__modelTypes:
      for s in self.__statsTypes:
        self.df_Peptides.loc[( prtnGrp[0], peptide ,'light'), s+'_'+m] = stats[s][m]
    
    return fig
  
  def abundancePlotAllPeptides(self, df, prtnGrp, wtype='both', labels=dict(), saveFigOpts=dict(), legendopt=dict() ):
    """
    Args:
        df (Dataframe): Pandas pivot table df
        prtnGrp (tuple): Protein-Gene index
        wtype (str, optional): weight_type: "heavy", "light", or "both". Defaults to 'both'.
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
      fig = self.abundancePlot1Peptide(fig=fig, df=df[[peptide]], prtnGrp=prtnGrp, peptide=peptide, lineopt = dict( color=cpalette[ colorcnt%colorcntmax ], width=2) )
      colorcnt += 1
    
    if len(fig.data) < 1 : return #  if nothing showing, skip
    
    fig.update_layout( 
      title=labels['title'], 
      xaxis_title=labels['x'], 
      yaxis_title=labels['y'],
      legend_title="Peptide",
      font=dict(
          family=labels['fontfamily'],
          size=labels['size'],
          color="Black" # "RebeccaPurple"
      ),
      legend=legendopt
    )
    
    proteingenename = prtnGrp[1] if type(prtnGrp[1]) == str else prtnGrp[0]

    
    if saveFigOpts['savePlot']:
      filename = "RelAbundance_Gene-"+ proteingenename +"-peptides" # title starts with "Gene: " or "Protein Group: "
      self.__savePlot(saveFigOpts, fig, filename)
    else: 
      fig.show()
    
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

    # import os
    folder = os.path.join( saveFigOpts['folder'], 'pngs' )
    # save static png from plotly GO, requires plotly GO needs kaleido installed
    if not os.path.exists(folder): os.mkdir(folder)
    filepath = os.path.join(folder,filename+'.png')
    fig.write_image( filepath )
    
    # save plotly GO as interactive graph in html
    # 
    # <script src="./js/plotly-1.58.5.min.js" charset="utf-8"></script>
    # <script src="./js/plotly-2.20.0.min.js" charset="utf-8"></script>
    # <script src="https://cdn.plot.ly/plotly-2.20.0.min.js" charset="utf-8"></script>
    # <script src="https://cdn.plot.ly/plotly-latest.min.js" charset="utf-8"></script>
    folder = os.path.join( saveFigOpts['folder'], 'htmls' )
    if not os.path.exists(folder): os.mkdir(folder)
    filepath = os.path.join(folder, filename+".html")
    incPlotlyJs = "plotly.min.js"  # False # add cdn plotly js on the html directly to minimize size of each html # symbolic link
    fig.write_html( filepath, include_plotlyjs=incPlotlyJs )
    return
      
  def abundancePlotProteinLevel(self, df, prtnGrp, labels=dict(), saveFigOpts = dict(), legendopt=dict() ):
    """
    Args:
        df (Dataframe): Pandas pivot table df
        trendlines (dict, optional): show, solid (bool), color (str), width (float) 
        lines (dict, optional): show, solid (bool), color (str), width (float) 
        markers (dict, optional): show (bool), symbol (str), size (float) 
    return: plotly plot objecta
    """
    # import numpy as np
    # import pandas as pd
    # import ExpoDecayFit as edf
    # from ExpoDecayFit import ExpoDecayFit 
    labels = self.__setArgLabels(labels=labels)
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    legendopt = self.__setArgLegendopt(legendopt=legendopt)
    # lines = self.__setArgLines(lines=lines)
    # trendlines = self.__setArgTrendlines(trendlines=trendlines)
    # markers = self.__setArgMarkers(markers=markers)
    #
    if ( type(df.columns) == pd.core.indexes.multi.MultiIndex ): # Pg_level needing to average all peptides
      df.columns = df.columns.reorder_levels(['wtype','Peptide'])
      types=tuple(df.columns.levels[0]) # essentially ('heavy','light')
      types=('heavy','light')
      df_avg = pd.DataFrame()
      for t in types: df_avg[t] = df[t].mean(axis = 1) # na will be ignored !! Re-think strategy
    else:
      df_avg = df.copy() # peptide level data
      
    df_avg.reset_index(inplace=True) # reset index to have 'time' column, for plotting
    
    # Obtain exponential fit params for the light series
    # expfitres = self.__expFitAvgLightSeries(df_avg) # {b, t12, r2}
    # b = expfitres['b']
    # t12 = expfitres['t12']
    # r2 = expfitres['r2']
    # # determine x-range. If t12 is too long, try to do something...
    # xrange = round( min( max(6, 1.8*t12) , 10) ) # no more than 10, between 6 and 10. If t12 is close, show 1.8*t12
    # xsamples = np.linspace(start = 0, stop = xrange, num = 300)
    # ysamples = dict()
    # ysamples['light'] = 100*np.exp(-b*xsamples)
    # ysamples['heavy'] = 100 - ysamples['light']
    # decayf =
    #
    # now do these with ExpoDecayFit module
    model = edf.ExpoDecayFit(df_avg, xAxisName = self.__xAxisName, modelTypes=self.__modelTypes, statsTypes=self.__statsTypes) # model.startx, starty, samplexs, sampleys
    bs = model.modelsummary.loc['b',:]
    t12s = model.modelsummary.loc['t12',:]
    r2s = model.modelsummary.loc['r2',:]
    
    xsamples = model.samplexs
    ysamples = dict()
    ysamples['light'] = model.sampleys
    ysamples['heavy'] = 100 - ysamples['light']
    
    xrange = round( min( max(6, 1.8*t12s['LnLM1']) , 10) ) # no more than 10, between 6 and 10. If t12 is close, show 1.8*t12
    
 
    colors = dict(heavy='rgba(199,10,165,.9)', light='rgba(56,233,99,.9)')
    # symbols = dict(heavy='x', light='cross') # try hexagon
    symbols = dict(heavy='circle', light='circle') # try hexagon

    fig = go.Figure()
    for t in types:
      markeropt = dict(color=colors[t], symbol=symbols[t])
      specs = dict(mode='markers', name=t.capitalize(), showlegend=False, connectgaps=False)
      fig = self.__add1goTrace(fig, x=df_avg[self.__xAxisName], y=df_avg[t], specs=specs, markeropt=markeropt)
      legendname = t.capitalize()+' (synthesis)' if t=='heavy' else t.capitalize()+' (degradation)' if t=='light' else '-'
      specs = dict(mode='lines', name=legendname, showlegend=True, connectgaps=True)
      fig = self.__add1goTrace(fig, x=xsamples, y=ysamples[t]['LnLM1'], specs=specs, markeropt=markeropt )
      fig = self.__add1goTrace(fig, x=xsamples, y=ysamples[t]['CFit'], specs=specs, markeropt=markeropt )
        
    if len(fig.data) < 1 : return #  if nothing showing, skip
    
    # show half life if within range
    if t12s['LnLM1'] < xrange: fig.add_vline(x=t12s['LnLM1'], line_width=1, line_dash="dash", line_color="black", annotation_text="&nbsp;<b>t<sub>½</sub></b> = "+str(t12s['LnLM1'].__round__(2)), annotation_position='bottom right')
    
    fig.update_layout( title=labels['title'], xaxis_title=labels['x'], yaxis_title=labels['y'],
      # legend_title="_",
      font=dict(
          family=labels['fontfamily'],
          size=labels['size'],
          color="Black" # "RebeccaPurple"
      ), 
      legend=legendopt
    )
    
    # update df_PgGeneDescSummary table
    self.df_PgGeneDescSummary.loc[prtnGrp,'chart'] = 1
    
    proteingenename = prtnGrp[1] if type(prtnGrp[1]) == str else prtnGrp[0]
    
    if saveFigOpts['savePlot']:
      filename = "RelAbundance_Gene-"+ proteingenename 
      self.__savePlot(saveFigOpts, fig, filename)
    else: 
      fig.show()
    
    return fig
  
  # def __addExpTrendline(self, fig, df, peptide, wtype='both', lines=dict(), markers=dict(), trendlines = dict()): # add trendline for the "light" peptide with exponential decay trendline
  #   # if light is present, fit it with exponential trendline
  #   if df.loc[:"light"] : pass
  #   return fig
  
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
    
    df = self.df_PgPivot[[prtnGrp[0]]] # filter only one prtnGrp, with both light and heavy data, can have multiple peptides, 
    df.columns = df.columns.droplevel(0) # remove column index ProteinGroup
    
    # REMOVE peptides with more than __maxNAcnt na values
    df = df.drop(columns=[col for col in df if df[col].isna().sum() > self.__maxNAcnt ]) # __maxNAcnt is set globally
    if df.shape[1] < 1: return None # nothing to plot
    
    labels['title'] = 'Gene: ' + prtnGrp[1] if type(prtnGrp[1]) == str else 'Protein Group: '+prtnGrp[0]
    _ = self.abundancePlotAllPeptides(df, prtnGrp=prtnGrp, wtype='both', labels=labels, saveFigOpts=saveFigOpts)
    
    # Now plot protein level average with trendline
    # labels['title'] = "Average for Protein: "+prtnGrp
    labels['title'] = "Average for "+ labels['title']
    _ = self.abundancePlotProteinLevel(df, prtnGrp=prtnGrp, labels=labels, saveFigOpts=saveFigOpts)

    return 
  
    # Now plot with matplotlib
    # plt.figure().set_figwidth(15)
    # ax = df.plot(kind='line', marker='.')
    # handles, labels = ax.get_legend_handles_labels()  
    # lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(0,-0.15))
    # natext = "("+str(maxNAcnt)+ '-NAs/series included)' if maxNAcnt>1 else '(1-NA/series included)' if maxNAcnt>0 else ' '
    # axtext = ax.text(-0.13, -0.13, natext, transform=ax.transAxes)
    # # ax.grid('on')
    # ax.set_title(charttitle)
    # ax.set_xlabel(xlabel)
    # ax.set_ylabel(ylabel)
    # # ax
    
    # Save png
    if (savePlot): 
      #set path
      # import os
      filepath = os.path.join(saveFolder,"relAbundance_prtnGrp_"+ prtnGrp.replace(";","_") +"_"+str(maxNAcnt)+"NAs.png")
      plt.savefig(filepath, bbox_extra_artists=(lgd,axtext,axtext), bbox_inches='tight', facecolor=(.65, .75, .75))
    
    return ax 
  
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
    plotmax = 4 # in/out
    # for prtnGrp in self.PgList:
    for prtnGrp in self.df_PgGeneDescSummary.index.values: # multi-level index with (protein, gene)
      if plotmax == 0 : return # in/out
      self.abundancePlot1Pg(prtnGrp=prtnGrp, labels=labels, saveFigOpts=saveFigOpts)
      plotmax -= 1 # in/out
    return
  
  def calcProteinDecayConstFromPeptides(self):
    # after calculating b-constants and half-lives of each peptide, find the average for the protein group.
    # Use mean value among the b-constants, which is same as using harmonic mean of all the half-lives.
    # dfpep.reset_index(inplace=True)
    # dfpep = dfpep[dfpep['wtype']=='light']
    # dfpep.shape
    return

#%%
# file = os.path.join(os.getcwd(), "../data/iMN_Peptide_Dataset.xlsx") # assuming cwd is .../Visualization/src/ folder
# file = os.path.join(os.getcwd(), "../data/20230522_dSILAC_Turnover_LightRelativeAbundances.xlsx") # assuming cwd is .../Visualization/src/ folder
file = os.path.join(os.getcwd(), "../data/06202023_FinalReport_dSILAC_iMN_MultiTimePoint.xlsx") # assuming cwd is .../Visualization/src/ folder
# new file 06202023_FinalReport_dSILAC_iMN_MultiTimePoint.xlsx

pto = ProteinTurnover(filepath= file)
# pto.chkDfPgIndexUnique()
# pto.chkDfPgGeneSummaryIndexUnique()
#%%
# savenow = True
savenow = False
savePath = "../media/plots/"

#%%
# pto.abundancePlotPgAll(savePlot=savenow, saveFolder=savePath)
pto.abundancePlotPgAll( saveFigOpts = dict(savePlot=savenow, folder=savePath) )


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

# .

#%%