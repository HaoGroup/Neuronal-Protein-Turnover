#%%
import matplotlib.pyplot as plt
import pandas as pd
import os
from sklearn.linear_model import LinearRegression
import plotly.graph_objects as go
# To export static images from go, need pip install -U kaleido
import plotly.express as px
# import seaborn as sns
# import numpy as np
import streamlit as st

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
  def __init__(self, filepath, x_values = [0,1,2,4,6]) -> None:
    """_summary_

    Args:
        filepath (_type_): the source file location. Use os to ensure right format
        x_values (list, optional): the time marks on x-axis. Defaults to [0,1,2,4,6].
        x_unit (str, optional): the unit for x-axis. Defaults to "day".
    """
    self.__yAxisName = 'Relative abundance'
    self.__yUnit = 'L%'
    self.__xAxisName = 'time'
    self.__xUnit = 'day'
    # self.__compoIndex = ['PG.ProteinGroups', 'Peptide', 'wtype'] # basic composite index used in DFs
    self.__compoIndex = ['Protein', 'Peptide', 'wtype'] # basic composite index used in DFs
    self.__maxNAcnt = 1 # Each series only has at most 4 data points at day = 1,2,4,6.  Only allow at most 1 missing to be plot
    self.df_PgRaw = None # initialize, current structure has columns: ['PG.ProteinGroups', 'PG.Genes', 'PG.ProteinDescriptions', 'Peptide', 'iMN_Day{x}_{Light/Heavy}_Relative_Abundance', 'k_results', 'Protein_Turnover_from_k', 'residuals']
    self.df_Pg = None # initialize, cleaned and re-structured df for analysis ['PG.ProteinGroups', 'PG.Genes', 'PG.ProteinDescriptions', 'Peptide', 'wtype', 0, 1, 2, 4, 6, 'k_results', 'Protein_Turnover_from_k', 'residuals']
    self.df_PgGeneDescLookup = None # initialize, serve as lookup table between PgProteinGroups, PgGenes, and PgProteinDescriptions, PLUS list/tuple of peptides
    self.PgList = [] # initialize, the list of ProteinGroups in the dataset
    self.df_PgPivot = [] # initialize, for plotting
    self.__ingestData(filepath=filepath, x_values=x_values) # set df_PgRaw and df_PG
    self.__setPgPivot()
    self.__setPgGeneDescLookup() # set df_PgGeneDescLookup as well as PgList, but need df_PgPivot first
    return
  
  def __ingestData(self, filepath, x_values = None):
    """
    Args:
        filepath (_type_): the source file location. Use os to ensure right format
        x_values (list, optional): the time marks on x-axis. Defaults to [0,1,2,4,6].
        x_unit (str, optional): the unit for x-axis. Defaults to None.

    Returns: None, will set df_PgRaw and df_Pg
    """
    self.df_PgRaw = self.df_PgRaw = pd.read_excel(filepath) if (filepath[-5:] == '.xlsx' or filepath[-4:] == '.xls') else pd.read_csv(filepath) if (filepath[-4:] == '.csv') else None
    # print(self.df_PgRaw.head())
    
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
    import re
    df_light.rename(columns=lambda x: re.sub(r'iMN_Day(\d).*_Relative_Abundance$', r'\1', x), inplace = True)
    df_heavy.rename(columns=lambda x: re.sub(r'iMN_Day(\d).*_Relative_Abundance$', r'\1', x), inplace = True)
    # df_light.rename(columns=lambda x: re.sub(r'Relative.*Intensity.Day.(\d)$', r'\1', x), inplace = True) # Not used after renamed new xlxs column names
    # df_heavy.rename(columns=lambda x: re.sub(r'Relative.*Intensity.Day.(\d)$', r'\1', x), inplace = True) # Not used after renamed new xlxs column names
    # set composite index 
    df_light.set_index(self.__compoIndex, inplace=True)  
    df_heavy.set_index(self.__compoIndex, inplace=True)  
    # stack the two df together
    self.df_Pg = pd.concat([df_light, df_heavy])
    # print(self.df_Pg.head())
    # self.chkDfPgIndexUnique()
    if not self.chkDfPgIndexUnique(): print("Protein Group in df_Pg not unique")
    
    return None
  
  def __setPgPivot(self):
    """
    set up pivot table from df_PG table for analysis and plots
    Returns: None
    """
    # df_PgPivot table is indexed by the time/day value starting at 0.
    xAxisName = self.__xAxisName
    # select only rel_Abundance columns
    df = self.df_Pg.loc[:, self.df_Pg.columns.str.contains('^\d$')] # only columns with single digit as colnames
    df = df.stack(dropna = False) # melt method does not keep the original index keys. df is now a pd series
    df.index.names = df.index.names[:-1] + [xAxisName] # name new index column
    df.name = self.__yAxisName # set name of the series (column)
    df = pd.DataFrame(df).reset_index() # from pd series to dataframe, and reset index to use pivot functions
    # df = pd.melt()
    df[xAxisName] = df[xAxisName].astype(int) # need these x-values be taken as numeric
    self.df_PgPivot = 100 * df.pivot_table(index=xAxisName, columns = self.__compoIndex, values = self.__yAxisName)
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
  
  def __setPgGeneDescLookup(self):
    """
    setting up the df_PgGeneDescLookup as lookup table self.df_PgGeneDescLookup as well as 
    setting up the list of ProteinGroup values as attribute self.PgList
    Returns: None
    """
    # set PgList
    # self.PgList = self.df_PgRaw["PG.ProteinGroups"].unique()
    self.PgList = self.df_PgRaw["Protein"].unique()
    # set df_PgGeneDescLookup
    # compoIndex = ['PG.ProteinGroups', 'PG.Genes']
    compoIndex = ['Protein', 'Gene.x']
    # selectedCols = compoIndex + ['PG.ProteinDescriptions']
    selectedCols = compoIndex + ['Protein.Description']
    df = self.df_PgRaw[selectedCols]
    df.set_index(compoIndex, inplace=True)    
    self.df_PgGeneDescLookup = df.drop_duplicates()
    if not self.chkDfPgGeneLookupIndexUnique(): print("Protein Group - GeneId combo not unique")
    
    # next add list (tuple) of peptides to each row of prtnGrp
    self.df_PgGeneDescLookup.loc[:,"peptides"] = self.df_PgGeneDescLookup.apply(lambda x: self.peptidesFromPrtnGrp(x.name[0]), axis = 1)
    # Warning: A value is trying to be set on a copy of a slice from a DataFrame. 
    # False positive: https://stackoverflow.com/questions/42105859/pandas-map-to-a-new-column-settingwithcopywarning.
    return
  
  def chkDfPgIndexUnique(self): 
    # print(f'df_Pg composite index ({", ".join( self.df_Pg.index.names )}) is unique? {self.df_Pg.index.is_unique}') 
    return self.df_Pg.index.is_unique
  
  def chkDfPgGeneLookupIndexUnique(self): 
    # print(f'df_PgGeneDescLookup composite index ({", ".join( self.df_PgGeneDescLookup.index.names )}) is unique? {self.df_PgGeneDescLookup.index.is_unique}')
    return self.df_PgGeneDescLookup.index.is_unique
  
  def __setArgLabels(self, labels=dict(x=None, y=None, title=None) ):
    """
    setting generic keyword argument labels into x, y, and title
    Args:
        labels (dict, optional): x-, y-lables, and title of chart. Defaults to dict(x=None, y=None, title=None).
    Returns:
        dict: { x, y, title }
    """
    res = labels.copy()
    res['x'] = res['x'] if (res.__contains__('x') and res['x']) else self.__xAxisName.capitalize() + ' ('+ self.__xUnit + ')'
    res['y'] = res['y'] if (res.__contains__('y') and res['y']) else self.__yAxisName.capitalize() + ' ('+ self.__yUnit + ')'
    res['title'] = res['title'] if (res.__contains__('title') and res['title']) else None
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
  
  def __setArgSaveFigs(self, saveFigs=dict(savePng=False, folder=None) ):
    """
    setting generic keyword argument saveFigs into savePng and folder
    Args:
        saveFigs (dict, optional): savePng (binary) and folder (str). Defaults to dict(savePng=False, folder=None).
    Returns:
        dict: { savePng, folder }
    """
    res = saveFigs.copy()
    res['savePng'] = res['savePng'] if ( res.__contains__('savePng') and isinstance(res['savePng'], bool) ) else False
    res['folder'] = res['folder'] if ( res.__contains__('folder') and res['folder'] ) else './'
    return res
  
  def __set1PgChartTitle(self, prtnGrp):
    rows = self.df_PgGeneDescLookup.loc[prtnGrp,:] # normally should have exactly one row
    n = len(rows.index)
    if n == 1 : 
      title = 'Gene: ' + rows.index[0]
    elif n > 1 :
      allgenes = [ rows.index[i] for i in range(n) ]
      title = 'Genes (multiple): ' + ','.join(allgenes)
    else:
      title = 'Gene cannot be determined.'
    # title.replace(";","_") # can allow semi-colons in titles. Not in filenames (png)
    
    return title
  
  def __add1goTrace(self, fig, x, y, specs=dict(), lineopt=dict(), markeropt=dict() ):
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
    
    return fig
  
  def abundancePlot1Peptide(self, fig, df, peptide, wtype='both', lineopt=dict(), markeropt=dict() ):
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
    #
    df.columns = df.columns.droplevel(0) # further drop peptide level header 
    # if ( isinstance(df.columns, pd.core.indexes.multi.MultiIndex) and len(df.columns.levels) > 1 ) : df.columns = df.columns.droplevel(0) # was using this function for some other cases. No longer need this check.
    df.reset_index(inplace=True) # reset index to have 'time' column, for plotting
    # 
    series = ('light') if wtype == 'light' else ('heavy') if wtype == 'heavy' else ('heavy','light')
    for i, wtype in enumerate(series):
      legend = False if i>0 else True, # only show when i=0 # somehow this line is resulting legend as a tuple (True, )
      legend = legend[0] if isinstance(legend,tuple) else legend if isinstance(legend,bool) else True
      specs = dict(name=peptide, connectgaps=True, showlegend=legend, mode='lines+markers')
      
      fig = self.__add1goTrace(fig, x=df[self.__xAxisName], y=df.loc[:,wtype], specs=specs, lineopt=lineopt, markeropt=markeropt )
    
    return fig
  
  def abundancePlotAllPeptides(self, df, wtype='both', labels=dict(), saveFigOpts=dict() ):
    """
    Args:
        df (Dataframe): Pandas pivot table df
        wtype (str, optional): weight_type: "heavy", "light", or "both". Defaults to 'both'.
        labels (dict, optional): x-, y-labels and title. Defaults to empty dictionary
        lines (dict, optional): show, solid (bool), color (str), width (float) 
        markers (dict, optional): show (bool), symbol (str), size (float) 
    return: plotly plot object
    """    
    labels = self.__setArgLabels(labels=labels)
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)

    peptides = self.__peptidesFromPivottable(df)

    cpalette = px.colors.qualitative.Dark24  # trendline use '#00FE35', which is Light24[1] color, lime green
    colorcnt = 0 # initialize # can consider making this global if a random start is preferred.
    colorcntmax = len(cpalette)    
    fig = go.Figure()
    for peptide in peptides: 
      fig = self.abundancePlot1Peptide(fig=fig, df=df[[peptide]], peptide=peptide, lineopt = dict( color=cpalette[ colorcnt%colorcntmax ], width=2) )
      colorcnt += 1
    
    if len(fig.data) < 1 : return #  if nothing showing, skip
    
    fig.update_layout( title=labels['title'], xaxis_title=labels['x'], yaxis_title=labels['y'],
      legend_title="Peptide",
      font=dict(
          # family="Courier New, monospace",
          # size=18,
          color="Black" # "RebeccaPurple"
      )
    )
    if saveFigOpts['savePng']:
      # import os
      folder = saveFigOpts['folder']
      # save static png from plotly GO, requires plotly GO needs kaleido installed
      if not os.path.exists(folder): os.mkdir(folder)
      filename = "RelAbundance_prtnGrp-"+ labels['title'][6:].replace(";","-") +"-peptides" # title starts with "Gene: "
      filepath = os.path.join(folder,filename+'.png')
      fig.write_image( filepath )
      
      # save plotly GO as interactive graph in html
      folder = os.path.join( saveFigOpts['folder'], 'htmls' )
      if not os.path.exists(folder): os.mkdir(folder)
      filepath = os.path.join(folder, filename+".html")
      fig.write_html( filepath )
    else: 
      fig.show()
    
    return fig
  
  def __expFitLightSeries(self, df):
    """
    Create an exponential fit from the light series
    Args:
        df (pd DataFrame): df with x and y values to be fitted
    return: dict(b, t12, r2)
    """
    import numpy as np
    from sklearn.linear_model import LinearRegression
    # from scipy.optimize import curve_fit
    
    # df is pivot table with only 2 columns of averages. In case there is any NaN
    df0 = df.dropna()
    
    res = dict() # with b, t12, r2 values
    
    # Try scipy curve_fit:
    # from scipy.optimize import curve_fit
    # import numpy as np
    # x = df0['time']
    # y = np.log(0.01*df0['light'])
    # popt, pcov = curve_fit(lambda t, b: np.exp(-b * t), x, y) # popt-optimized parameters
    
    x=df0[['time']]
    y=np.log(0.01*df0['light'])
    lm = LinearRegression(fit_intercept=False)
    lm.fit(x, y)
    res['b'] = 0 - lm.coef_[0]
    res['t12'] = np.log(2)/res['b']
    res['r2'] = lm.score(x, y)
    
    return res
  
  def abundancePlotProteinLevel(self, df, labels=dict(), saveFigOpts = dict() ):
    """
    Args:
        df (Dataframe): Pandas pivot table df
        trendlines (dict, optional): show, solid (bool), color (str), width (float) 
        lines (dict, optional): show, solid (bool), color (str), width (float) 
        markers (dict, optional): show (bool), symbol (str), size (float) 
    return: plotly plot objecta
    """
    import numpy as np
    labels = self.__setArgLabels(labels=labels)
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    # lines = self.__setArgLines(lines=lines)
    # trendlines = self.__setArgTrendlines(trendlines=trendlines)
    # markers = self.__setArgMarkers(markers=markers)
    #
    df.columns = df.columns.reorder_levels(['wtype','Peptide'])
    # types=tuple(df.columns.levels[0]) # essentially ('heavy','light')
    types=('heavy','light')
    df_avg = pd.DataFrame()
    for t in types: df_avg[t] = df[t].mean(axis = 1) # na will be ignored !! Re-think strategy
    df_avg.reset_index(inplace=True) # reset index to have 'time' column, for plotting
    
    # Obtain exponential fit params for the light series
    expfitres = self.__expFitLightSeries(df_avg) # {b, t12, r2}
    b = expfitres['b']
    t12 = expfitres['t12']
    r2 = expfitres['r2']
    # determine x-range. If t12 is too long, try to do something...
    xrange = round( min( max(6, 1.8*t12) , 10) ) # no more than 10, between 6 and 10. If t12 is close, show 1.8*t12
    xsamples = np.linspace(start = 0, stop = xrange, num = 300)
    ysamples = dict()
    ysamples['light'] = 100*np.exp(-b*xsamples)
    ysamples['heavy'] = 100 - ysamples['light']
    # decayf =
    
    colors = dict(heavy='rgba(199,10,165,.9)', light='rgba(56,233,99,.9)')
    symbols = dict(heavy='x', light='cross') # try hexagon

    fig = go.Figure()
    for t in types:
      markeropt = dict(color=colors[t], symbol=symbols[t], size=7)
      specs = dict(mode='markers', name=t.capitalize(), showlegend=False, connectgaps=False)
      fig = self.__add1goTrace(fig, x=df_avg[self.__xAxisName], y=df_avg[t], specs=specs, markeropt=markeropt )
      legendname = t.capitalize()+' (synthesis)' if t=='heavy' else t.capitalize()+' (degradation)' if t=='light' else '-'
      specs = dict(mode='lines', name=legendname, showlegend=True, connectgaps=True)
      fig = self.__add1goTrace(fig, x=xsamples, y=ysamples[t], specs=specs, markeropt=markeropt )
        
    if len(fig.data) < 1 : return #  if nothing showing, skip
    
    # show half life if within range
    if t12 < xrange: fig.add_vline(x=t12, line_width=1, line_dash="dash", line_color="black", annotation_text="&nbsp;<b>t<sub>½</sub></b> = "+str(t12.__round__(2)), annotation_position='bottom right' )
    
    fig.update_layout( title=labels['title'], xaxis_title=labels['x'], yaxis_title=labels['y'],
      # legend_title="_",
      font=dict(
          # family="Courier New, monospace",
          # size=18,
          color="Black" # "RebeccaPurple"
      )
    )
    if saveFigOpts['savePng']:
      # import os
      folder = saveFigOpts['folder']
      # save static png from plotly GO, requires plotly GO needs kaleido installed
      if not os.path.exists(folder): os.mkdir(folder)
      filename = "RelAbundance_prtnGrp-"+ labels['title'][21:].replace(";","-") # title starts with "Average for Protein: "
      filepath = os.path.join(folder, filename+".png")
      fig.write_image( filepath )
      
      # save plotly GO as interactive graph in html
      folder = os.path.join( saveFigOpts['folder'], 'htmls' )
      if not os.path.exists(folder): os.mkdir(folder)
      filepath = os.path.join(folder, filename +".html")
      fig.write_html( filepath )
    else: 
      fig.show()
    
    return fig
  
    # only keep rows with at most maxNAcnt nan values
    df = df.drop(columns=[col for col in df if df[col].isna().sum() > self.__maxNAcnt ]) # __maxNAcnt is set globally
    df.reset_index(inplace=True) # reset index to have 'time' column, for plotting
    if df.shape[1] < 2: return fig # return fig # other than the time column, nothing to plot
    # 
    series = ('light') if wtype == 'light' else ('heavy') if wtype == 'heavy' else ('heavy','light')
    for i, wtype in enumerate(series):
      fig.add_trace(go.Scatter(
      x=df[self.__xAxisName],
      y=df.loc[:,wtype],
      line=dict(color=lines['color'], width=lines['width']), 
      marker_symbol = markers['symbol'], 
      marker_size = markers['size'],
      showlegend = False if i else True, # only show when i=0
      name = peptide, # Style name/legend entry with html tags
      connectgaps=True # override default to connect the gaps
      ))
    
    # if trendlines['show']: fig = self.__addExpTrendline(fig, df, peptide, wtype=wtype, lines=lines, markers=markers, trendlines = trendlines) # df has one time column with light and/or heavy
    
    return fig

  def __addExpTrendline(self, fig, df, peptide, wtype='both', lines=dict(), markers=dict(), trendlines = dict()): # add trendline for the "light" peptide with exponential decay trendline
    # if light is present, fit it with exponential trendline
    if df.loc[:"light"] : pass
    return fig
  
  def abundancePlot1Pg(self, prtnGrp, labels=dict(), saveFigs = dict() ):
    """
    Args:
        prtnGrp (str): the ProteinGrp being processed
        labels (dict, optional): x-, y-labels and title. Defaults to empty dictionary
        saveFigs (dict, optional): savePng (binary) and folder (str). Defaults to dict(savePng=False, folder=None).
    return: None
    """
    labels = self.__setArgLabels(labels=labels)
    saveFigs = self.__setArgSaveFigs(saveFigs=saveFigs)
    
    df = self.df_PgPivot[[prtnGrp]] # filter only one prtnGrp, with both light and heavy data, can have multiple peptides, 
    df.columns = df.columns.droplevel(0) # remove column index ProteinGroup
    
    # REMOVE peptides with more than __maxNAcnt na values
    df = df.drop(columns=[col for col in df if df[col].isna().sum() > self.__maxNAcnt ]) # __maxNAcnt is set globally
    if df.shape[1] < 1: return None # nothing to plot
    
    labels['title'] = labels['title'] if labels['title'] else self.__set1PgChartTitle(prtnGrp) # Title for such plots
    _ = self.abundancePlotAllPeptides(df, wtype='both', labels=labels, saveFigOpts=saveFigOpts)
    
    # Now plot protein level average with trendline
    labels['title'] = "Average for Protein: "+prtnGrp
    _ = self.abundancePlotProteinLevel(df, labels=labels, saveFigOpts=saveFigOpts)

    return 
  
    # Now plot with matplotlib
    # plt.figure().set_figwidth(15)
    ax = df.plot(kind='line', marker='.')
    handles, labels = ax.get_legend_handles_labels()  
    lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(0,-0.15))
    natext = "("+str(maxNAcnt)+ '-NAs/series included)' if maxNAcnt>1 else '(1-NA/series included)' if maxNAcnt>0 else ' '
    axtext = ax.text(-0.13, -0.13, natext, transform=ax.transAxes)
    # ax.grid('on')
    ax.set_title(charttitle)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax
    
    # Save png
    if (savepng): 
      #set path
      import os
      filepath = os.path.join(saveFolder,"relAbundance_prtnGrp_"+ prtnGrp.replace(";","_") +"_"+str(maxNAcnt)+"NAs.png")
      plt.savefig(filepath, bbox_extra_artists=(lgd,axtext,axtext), bbox_inches='tight', facecolor=(.65, .75, .75))
    
    return ax 
  
  def abundancePlotPgAll(self, labels=dict(), saveFigs=dict()): 
    """
    Args:
        labels (dict, optional): x-, y-labels and title. Defaults to empty dictionary
        saveFigs (dict, optional): savePng (binary) and folder (str). Defaults to dict(savePng=False, folder=None).
    return: None
    """
    # assumes labels and saveFigs are in the right forms.
    labels = self.__setArgLabels(labels=labels)
<<<<<<< Updated upstream
    saveFigs = self.__setArgSaveFigs(saveFigs=saveFigs)
    # savePng = saveFigs['savePng']; saveFolder = saveFigs['folder']; 
    plotmax = 10
    for prtnGrp in self.PgList:
      if plotmax == 0 : return
      self.abundancePlot1Pg(prtnGrp=prtnGrp, labels=labels, saveFigs=saveFigs)
      plotmax -= 1
=======
    saveFigOpts = self.__setArgSaveFigOpts(saveFigOpts=saveFigOpts)
    # savePng = saveFigOpts['savePng']; saveFolder = saveFigOpts['folder']; 
    # plotmax = 10
    for prtnGrp in self.PgList:
      # if plotmax == 0 : return
      self.abundancePlot1Pg(prtnGrp=prtnGrp, labels=labels, saveFigOpts=saveFigOpts)
      # plotmax -= 1
>>>>>>> Stashed changes
    return

#%%
# file = os.path.join(os.getcwd(), "../data/iMN_Peptide_Dataset.xlsx") # assuming cwd is .../Visualization/src/ folder
file = os.path.join(os.getcwd(), "../data/20230522_dSILAC_Turnover_LightRelativeAbundances.xlsx") # assuming cwd is .../Visualization/src/ folder

pto = ProteinTurnover(filepath= file)
# pto.chkDfPgIndexUnique()
# pto.chkDfPgGeneLookupIndexUnique()
#%%
savenow = True
# savenow = False
savePath = "../media/plots/"

#%%
# pto.abundancePlotPgAll(savepng=savenow, saveFolder=savePath)
pto.abundancePlotPgAll( saveFigs = dict(savePng=savenow, folder=savePath) )


# %%
# Note, the new dataset (20230522_dSILAC_Turnover_LightRelativeAbundances.xlsx) has these six duplicate combos
# Protein	Peptide	  							
# Q6ZTK2	GLQLEGELEELRQDR	
# Q6ZSR9	QSIAGSVSITSLSSR	
# Q6ZSR9  TQNNLESDYLAR	
# Q6ZTK2	GLQLEGELEELRQDR	
# Q6ZSR9	QSIAGSVSITSLSSR	
# Q6ZSR9	TQNNLESDYLAR	
# .

#%%