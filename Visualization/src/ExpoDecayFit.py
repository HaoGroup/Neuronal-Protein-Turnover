#%%
import pandas as pd
import numpy as np
import os
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit

# import plotly.graph_objects as go
# import plotly.express as px
# import seaborn as sns
# import streamlit as st

#%%
# Create class ExpoDecayFit to create a consistent model fit
#
class ExpoDecayFit:
  """
  Exponential Decay Fit 
  """
  def __init__(self, df) -> None:
    """_summary_

    Args:
        df (pandas Dataframe): Expected to be from pivot table, with time as index, and Multi-index column either (wtype,peptide), or (wtype,average) for protein level
    """
    # self.__indexName = df.index.name # 'time'
    # self.__indexVals = tuple( df.index.values )
    # self.__indexX = pd.Series( self.__indexVals , name=self.__indexName) # 0,1,2,4,6
    self.__wtypes = tuple(df.columns.levels[0]) # should be ('light','heavy')
    
    self.fitxweights = np.ones( df.index.shape ) # wights for the different x values # dafaults all ones.
    self.fityweights = np.array( () ) # wights for the different y SERIES. For example, a series with 1 NA out of the four values (time = 1,2,4,6) should have weights 3/4...
    self.fitxs = df.index.values # np.array() # x-data, 1d-array or pandas series
    self.fitys = np.array( () ) # y-data, typically 2d-array or pandas dataframe. Can have just average light and heavy two columns, or more if multiple peptides are included.
    
    self.__modelSkLR = None # sklearn.linearRegression
    self.__modelSmLR = None # statsmodel.linearRegression
    self.__modelSpCF = None # scipy curve_fit
    self.__modelNpPF = None # numpy polyfit (polynomial fit, like linear regression) 
    self.__modelTypes = ("SkLR", "SmLR", "SpCF", "NpPF")
    self.modelType = "SkLR" # current selected // default
    self.__modelStatsTypes = ('b','t12','r2')
    self.modelsummary = pd.DataFrame(index = pd.Series( self.__modelStatsTypes , name="stats"), columns=pd.MultiIndex.from_product([self.__wtypes, self.__modelTypes], names=['wtype','modeltype'])) # Dataframe, with two-level column headers ['light'/'heavy', modeltype]
    
    self.startx = 0
    self.starty = 0
    self.maxx = 10
    self.maxy = 100 # default shows 100%
    self.sampleN = 300 
    self.samplexs = self.__setSampleXs()
    self.sampleys = pd.DataFrame(index = range(self.sampleN), columns=pd.MultiIndex.from_product([self.__wtypes, self.__modelTypes], names=['wtype','modeltype'])) # Dataframe, with two-level column headers ['light'/'heavy', modeltype].
        
    return
  
  def __setSampleXs(self): return np.linspace(start=self.startx, stop=self.maxx, num = self.sampleN)
  

  def __setArgLabels(self, labels=dict(x=None, y=None, title=None) ):
    """
    setting generic keyword argument labels into x, y, and title
    Args:
        labels (dict, optional): x-, y-lables, and title of chart. Defaults to dict(x=None, y=None, title=None).
    Returns:
        dict: { x, y, title }
    """
    res = labels.copy()
    res['x'] = res['x'] if (res.__contains__('x') and res['x']) else self.__xAxisName + ' ('+ self.__xUnit + ')'
    res['y'] = res['y'] if (res.__contains__('y') and res['y']) else self.__yAxisName + ' ('+ self.__yUnit + ')'
    res['title'] = res['title'] if (res.__contains__('title') and res['title']) else None
    return res
  
  # def __setArgLineopt(self, lineopt=dict(color=None, width=None) ):
  #   """
  #   setting generic keyword argument lines into show, solid, color, and width
  #   Args:
  #       lineopt (dict, optional): color (str), width (float), ..
  #   Returns:
  #       dict: { color, width, ... }
  #   """
  #   res = lineopt.copy()
  #   res['color'] = res['color'] if (res.__contains__('color') and res['color']) else '#000000'
  #   res['width'] = res['width'] if (res.__contains__('width') and res['width']) else 2
  #   # preserving the rest of attributes
  #   return res

  # def __setArgMarkeropt(self, markeropt=dict(color=None, symbol=None, size=None) ):
  #   """
  #   setting generic keyword argument markers into color, symbol, size
  #   Args:
  #       markeropt (dict, optional): color (str), symbol (str), size (float), ...
  #   Returns:
  #       dict: { color, symbol, size, ... }
  #   """
  #   res = markeropt.copy()
  #   # https://plotly.com/python/marker-style/ # hourglass, x-thin, ...
  #   res['symbol'] = res['symbol'] if (res.__contains__('symbol') and res['symbol']) else 'hourglass'
  #   res['color'] = res['color'] if (res.__contains__('color') and res['color']) else '#000000'
  #   res['size'] = res['size'] if (res.__contains__('size') and res['size']) else 10
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


  def __expFitAvgLightSeries(self, df):
    """
    Create an exponential fit from the light series
    Args:
        df (pd DataFrame): df with x and y values to be fitted. y is expected to have 3 columns: time, light and heavy
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
  
  def abundancePlotProteinLevel(self, df, labels=dict() ):
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
    expfitres = self.__expFitAvgLightSeries(df_avg) # {b, t12, r2}
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
    symbols = dict(heavy='hourglass', light='cross') # try hexagon

    fig = go.Figure()
    for t in types:
      markeropt = dict(color=colors[t], symbol=symbols[t], size=10)
      specs = dict(mode='markers', name=t, showlegend=False, connectgaps=False)
      fig = self.__add1goTrace(fig, x=df_avg[self.__xAxisName], y=df_avg[t], specs=specs, markeropt=markeropt )
      specs = dict(mode='lines', name=t, showlegend=False, connectgaps=True)
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
    fig.show()
    
    return fig
  

#%%
# file = os.path.join(os.getcwd(), "../data/iMN_Peptide_Dataset.xlsx") # assuming cwd is .../Visualization/src/ folder
# pto = ProteinTurnover(filepath= file)
# pto.chkDfPgIndexUnique()
# pto.chkDfPgGeneLookupIndexUnique()
# savenow = True
# savenow = False
# savePath = "../media/plots/"

#%%
# pto.abundancePlotPgAll(savepng=savenow, saveFolder=savePath)
# pto.abundancePlotPgAll( saveFigOpts = dict(savePng=savenow, folder=savePath) )


# %%
