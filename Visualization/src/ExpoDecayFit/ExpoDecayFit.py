#%%
import pandas as pd
import numpy as np
import statsmodels.api as sm
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
  # def __init__(self, row, xAxisName = 'Time', modelTypes=('LnLM1', 'LnLM2', 'CFit'), statsTypes=('b','t12','r2')) -> None:
  def __init__(self, row, xAxisName = 'Time', xvalues = ['1','2','4','6','8'] , modelTypes=('CFit',), statsTypes=('b','t12','r2'), filter = dict()) -> None:
  # def __init__(self, row, xAxisName = 'Time', xvalues = ['1','2','4','6'] , modelTypes=('CFit',), statsTypes=('b','t12','r2'), filter = dict() ) -> None:
    """_summary_

    Args:
        # df (Pandas Dataframe): Expected to be from pivot table, with Time as index, and Multi-index column either (wtype,peptide), or (wtype,average) for protein level
        row (Pandas Series): a row of data, from peptide typically. 
        xAxisName (str): x-axis name
        modelTypes (list/tuple): list of model types to be deployed
        statsTypes (list/tuple): list of statistical parameters to be recorded. Defaults are b: decay constant, t12: half-life, and r2: r-squared value of fit
        filter (dict): dictionary of filter property, ratio-tolerence, absolute-tolerence, etc. See __setFilter() function
    """
    # self.__indexVals = tuple( df.index.values )
    # self.__indexX = pd.Series( self.__indexVals , name=self.__indexName) # 0,1,2,4,6
    
    # self.fitxweights = np.ones( df.index.shape ) # wights for the different x values # dafaults all ones.
    # self.fityweights = np.array( () ) # wights for the different y SERIES. For example, a series with 1 NA out of the four values (Time = 1,2,4,6) should have weights 3/4...
    # self.fitxs = df.index.values # np.array() # x-data, 1d-array or pandas series
    # self.fitys = np.array( () ) # y-data, typically 2d-array or pandas dataframe. Can have just average light and heavy two columns, or more if multiple peptides are included.
    
    self.__xAxisName = xAxisName # match with def in ProteinTurnover class
    self.__xvalues = xvalues # match with def in ProteinTurnover class
    # self.__modelPolyFit = None # numpy polyfit (polynomial fit, like linear regression) 
    # self.__modelTypes = ('LnLM1', 'LnLM2', 'CFit', 'PolyFit')
    self.__modelTypes = modelTypes
    # self.modelType = 'CFit' # current selected // default
    self.__statsTypes = statsTypes
    # self.modelsummary = pd.DataFrame(index = pd.Series( self.__statsTypes , name='stats'), columns=pd.MultiIndex.from_product([self.__wtypes, self.__modelTypes], names=['wtype','modeltype'])) # Dataframe, with two-level column headers ['light'/'heavy', modeltype]
    self.modelsummary = pd.DataFrame(index = pd.Series( self.__statsTypes , name='stats'), columns=self.__modelTypes) # Dataframe, column headers modeltype. For light series only
    self.__filter = self.__setFilter(filter)
    
    # moved y_predicted outside of module now as it's quite straight foward
    # self.startx = 0
    # self.starty = 0
    # self.maxx = 10
    # self.maxy = 100 # default shows 100%
    # self.sampleN = 300 
    # self.samplexs = self.__setSampleXs()
    # self.sampleys = pd.DataFrame(index = range(self.sampleN), columns=self.__modelTypes) # Dataframe, column headers modeltype. For light series only
    # was
    # self.sampleys = pd.DataFrame(index = range(self.sampleN), columns=pd.MultiIndex.from_product([self.__wtypes, self.__modelTypes], names=['wtype','modeltype'])) # Dataframe, with two-level column headers ['light'/'heavy', modeltype].
    self.__expFitAvgLightSeries(row, model=modelTypes[0] )
        
    return
  
  # def __setSampleXs(self): return np.linspace(start=self.startx, stop=self.maxx, num = self.sampleN)
  
  # Curve fitting function # for Scipy curve_fitting
  def __expDecayFcn(self, x, lambd): return np.exp(-lambd * x)
  
  def __setFilter(self, filter):
    """
    require (1) last value * rtol > next value
    and-or
    require (2) last value + atol > next value 
    type: monotone
    Args:
        filter (dict): {'type': 'monotone', 'rtol':1.2, 'atol':0.1, 'rel-abs-and-or' : 'or' }
    """
    res = {'type': 'monotone', 'rtol':1.2, 'atol':0.15, 'rel-abs-and-or' : 'or' } # default
    for key in res.keys():
      res[key] = filter(key) if (filter.__contains__(key) and filter(key)) else res[key]
    return res
  
  def __filterCheck(self, y):
    """
    check if y (list) satisfy the filter condition
    Args:
        y (list): y values in the light series, without NAs
    """
    res = True
    for i in range(len(y)-2): # ignore first one (always start from y=1), and no need to check last one
      # assume type = monotone
      cond1 = self.__filter['rtol']*y[i+1] > y[i+2] 
      cond2 = self.__filter['atol']+y[i+1] > y[i+2] 
      res = res and cond1 and cond2 if self.__filter['rel-abs-and-or']=='and' else res and (cond1 or cond2) # default 'or'
      if not res: return False
    return res

  # def __setArgLabels(self, labels=dict(x=None, y=None, title=None) ):
  #   """
  #   setting generic keyword argument labels into x, y, and title
  #   Args:
  #       labels (dict, optional): x-, y-lables, and title of chart. Defaults to dict(x=None, y=None, title=None).
  #   Returns:
  #       dict: { x, y, title }
  #   """
  #   res = labels.copy()
  #   res['x'] = res['x'] if (res.__contains__('x') and res['x']) else self.__xAxisName + ' ('+ self.__xUnit + ')'
  #   res['y'] = res['y'] if (res.__contains__('y') and res['y']) else self.__yAxisName + ' ('+ self.__yUnit + ')'
  #   res['title'] = res['title'] if (res.__contains__('title') and res['title']) else None
  #   return res

  def __expFitAvgLightSeries(self, row, model = 'all'):
    """
    Create an exponential fit from the light series, with x=0, y=1 (100%)
    Args:
        # df (pd DataFrame): df with x and y values to be fitted. y is expected to have 3 columns: Time, light and heavy
        row (Pandas Series): a row of data, from peptide typically. 
    return: None
    """
    
    # Try scipy curve_fit:
    # from scipy.optimize import curve_fit
    # import numpy as np
    # x = df0[self.__xAxisName]
    # y = np.log(0.01*df0['light'])
    
    xvals = ['0']+self.__xvalues
    data = row[xvals].dropna()
    data.index = data.index.astype(int)
    # data.index = data.index.astype(float)
    x = data.index.values
    y = data.values
    logy=np.log(y.tolist())
    
    if (model=='LnLM1' or model=='all'):  # ScikitLearn LR
      thismodel = 'LnLM1'
      mdl = LinearRegression(fit_intercept=False)
      mdl.fit(x.reshape(-1,1), logy) # make sure x.shape is (n,1), instead of (n,)
      res_b = -mdl.coef_[0]
      # save result statistics
      self.modelsummary.loc[self.__statsTypes[0],thismodel] = res_b
      self.modelsummary.loc[self.__statsTypes[1],thismodel] = np.log(2)/res_b
      self.modelsummary.loc[self.__statsTypes[2],thismodel] = mdl.score(x.reshape(-1,1), logy)
      # save model curve fit data points
      # moved y_predicted outside of module now as it's quite straight foward
      # self.sampleys[thismodel] = np.exp(-res_b*self.samplexs) # express in percentages
      
    if (model=='LnLM2' or model=='all'):  # statsmodels LR
      thismodel = 'LnLM2'
      # X = sm.add_constant(X)
      # model = sm.glm(y, X, family=log)
      mdl = sm.OLS(logy, x)
      results = mdl.fit()
      res_b = -results.params[0]
      # save result statistics
      self.modelsummary.loc[self.__statsTypes[0],thismodel] = res_b
      self.modelsummary.loc[self.__statsTypes[1],thismodel] = np.log(2)/res_b
      self.modelsummary.loc[self.__statsTypes[2],thismodel] = results.rsquared
      # save model curve fit data points
      # moved y_predicted outside of module now as it's quite straight foward
      # self.sampleys[thismodel] = np.exp(-res_b*self.samplexs) # express in percentages

    if (model=='CFit' or model=='all'):  # Scipy CurveFit
      thismodel = 'CFit'
      # use filter to avoid max iter limitation. Can also apply to other models other than 'CFit' if needed.
      if not self.__filterCheck(y): return
      # yvals = y.tolist() # or list(y) # y is a pd.core.series.Series
      # xvals = x[self.__xAxisName].tolist() # change from pd.Dataframe to pd.Series to list
      # try:
      #   popt, pcov = curve_fit(self.__expDecayFcn, x.tolist(), y.tolist())
      # except Exception as e:
      #   print(f"Caught an exception: {e}")
      #   raise  
      popt, pcov = curve_fit(self.__expDecayFcn, x.tolist(), y.tolist())
      # popt, pcov = curve_fit(lambda t, b: np.exp(-b * t), x, y) # popt-optimized parameters
      res_b = popt[0]
      # save result statistics
      self.modelsummary.loc[self.__statsTypes[0],thismodel] = res_b
      self.modelsummary.loc[self.__statsTypes[1],thismodel] = np.log(2)/res_b
      ## calculate r^2 ourselves
      # residuals = y - self.__expDecayFcn(x, res_b)[self.__xAxisName] # x needs to be np array here, not just list, and results needs to be a series, not Dataframe since y is a series
      residuals = y - self.__expDecayFcn(x, res_b) # x needs to be np array here, not just list, and results needs to be a series, not Dataframe since y is a series
      ss_res = np.sum(residuals**2)
      ss_tot = np.sum((y-np.mean(y))**2) # y needs to be np array here, not just list
      r_squared = 1 - (ss_res / ss_tot)
      # from sklearn.metrics import r2_score
      # r_squared = r2_score(y, self.__expDecayFcn(xvals, res_b))
      self.modelsummary.loc[self.__statsTypes[2],thismodel] = r_squared
      # save model curve fit data points
      # moved y_predicted outside of module now as it's quite straight foward
      # self.sampleys[thismodel] = np.exp(-res_b*self.samplexs) # express in percentages
      
    # if (model=='PolyFit' or model=='all'):  # Numpy Polyfit
      # pass
    
    return 
  

    

# %%
