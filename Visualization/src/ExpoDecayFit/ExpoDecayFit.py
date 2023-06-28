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
  def __init__(self, df, xAxisName = 'Time', modelTypes=('LnLM1', 'LnLM2', 'CFit'), statsTypes=('b','t12','r2')) -> None:
    """_summary_

    Args:
        df (pandas Dataframe): Expected to be from pivot table, with Time as index, and Multi-index column either (wtype,peptide), or (wtype,average) for protein level
    """
    # self.__indexVals = tuple( df.index.values )
    # self.__indexX = pd.Series( self.__indexVals , name=self.__indexName) # 0,1,2,4,6
    
    # self.fitxweights = np.ones( df.index.shape ) # wights for the different x values # dafaults all ones.
    # self.fityweights = np.array( () ) # wights for the different y SERIES. For example, a series with 1 NA out of the four values (Time = 1,2,4,6) should have weights 3/4...
    # self.fitxs = df.index.values # np.array() # x-data, 1d-array or pandas series
    # self.fitys = np.array( () ) # y-data, typically 2d-array or pandas dataframe. Can have just average light and heavy two columns, or more if multiple peptides are included.
    
    self.__xAxisName = xAxisName # match with def in ProteinTurnover class
    # self.__modelPolyFit = None # numpy polyfit (polynomial fit, like linear regression) 
    # self.__modelTypes = ('LnLM1', 'LnLM2', 'CFit', 'PolyFit')
    self.__modelTypes = modelTypes
    # self.modelType = 'LnLM1' # current selected // default
    self.__statsTypes = statsTypes
    # self.modelsummary = pd.DataFrame(index = pd.Series( self.__statsTypes , name='stats'), columns=pd.MultiIndex.from_product([self.__wtypes, self.__modelTypes], names=['wtype','modeltype'])) # Dataframe, with two-level column headers ['light'/'heavy', modeltype]
    self.modelsummary = pd.DataFrame(index = pd.Series( self.__statsTypes , name='stats'), columns=self.__modelTypes) # Dataframe, column headers modeltype. For light series only
    
    self.startx = 0
    self.starty = 0
    self.maxx = 10
    self.maxy = 100 # default shows 100%
    self.sampleN = 300 
    self.samplexs = self.__setSampleXs()
    # self.sampleys = pd.DataFrame(index = range(self.sampleN), columns=pd.MultiIndex.from_product([self.__wtypes, self.__modelTypes], names=['wtype','modeltype'])) # Dataframe, with two-level column headers ['light'/'heavy', modeltype].
    self.sampleys = pd.DataFrame(index = range(self.sampleN), columns=self.__modelTypes) # Dataframe, column headers modeltype. For light series only
    self.__expFitAvgLightSeries(df)
        
    return
  
  def __setSampleXs(self): return np.linspace(start=self.startx, stop=self.maxx, num = self.sampleN)

  # Curve fitting function # for Scipy curve_fitting
  def __expDecayFcn(self, x, lambd): return np.exp(-lambd * x)

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

  def __expFitAvgLightSeries(self, df, model = 'all'):
    """
    Create an exponential fit from the light series
    Args:
        df (pd DataFrame): df with x and y values to be fitted. y is expected to have 3 columns: Time, light and heavy
    return: None
    """
    # import numpy as np
    # from sklearn.linear_model import LinearRegression
    # import statsmodels.api as sm
    # from scipy.optimize import curve_fit
    
    # df is pivot table with only 2 columns of averages. In case there is any NaN
    if ( len(df.columns.values) != 2): print(f'exp fit df length = {len(df.columns.values)} with columns = {df.columns.values}')
    df0 = df.dropna()
    
    # res = dict() # with b, t12, r2 values
    
    # Try scipy curve_fit:
    # from scipy.optimize import curve_fit
    # import numpy as np
    # x = df0[self.__xAxisName]
    # y = np.log(0.01*df0['light'])
    
    x=df0[[self.__xAxisName]]
    # y=0.01*df0['light']
    # x=df0.loc[:,self.__xAxisName] # this give 1-D array, not working in sklearn
    y=0.01*df0.loc[:, df.columns.values[-1] ] # convert back from percentages into proportions
    logy=np.log(y)
    
    if (model=='LnLM1' or model=='all'):  # ScikitLearn LR
      thismodel = 'LnLM1'
      mdl = LinearRegression(fit_intercept=False)
      mdl.fit(x, logy)
      res_b = -mdl.coef_[0]
      # save result statistics
      self.modelsummary.loc[self.__statsTypes[0],thismodel] = res_b
      self.modelsummary.loc[self.__statsTypes[1],thismodel] = np.log(2)/res_b
      self.modelsummary.loc[self.__statsTypes[2],thismodel] = mdl.score(x, logy)
      # save model curve fit data points
      self.sampleys[thismodel] = 100*np.exp(-res_b*self.samplexs) # express in percentages
      
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
      self.sampleys[thismodel] = 100*np.exp(-res_b*self.samplexs) # express in percentages

    if (model=='CFit' or model=='all'):  # Scipy CurveFit
      thismodel = 'CFit'
      yvals = y.tolist() # or list(y) # y is a pd.core.series.Series
      xvals = x[self.__xAxisName].tolist() # change from pd.Dataframe to pd.Series to list
      popt, pcov = curve_fit(self.__expDecayFcn, xvals, yvals)
      # popt, pcov = curve_fit(lambda t, b: np.exp(-b * t), x, y) # popt-optimized parameters
      res_b = popt[0]
      # save result statistics
      self.modelsummary.loc[self.__statsTypes[0],thismodel] = res_b
      self.modelsummary.loc[self.__statsTypes[1],thismodel] = np.log(2)/res_b
      ## calculate r^2 ourselves
      residuals = y - self.__expDecayFcn(x, res_b)[self.__xAxisName] # x needs to be np array here, not just list, and results needs to be a series, not Dataframe since y is a series
      ss_res = np.sum(residuals**2)
      ss_tot = np.sum((y-np.mean(y))**2) # y needs to be np array here, not just list
      r_squared = 1 - (ss_res / ss_tot)
      # from sklearn.metrics import r2_score
      # r_squared = r2_score(y, self.__expDecayFcn(xvals, res_b))
      self.modelsummary.loc[self.__statsTypes[2],thismodel] = r_squared
      # save model curve fit data points
      self.sampleys[thismodel] = 100*np.exp(-res_b*self.samplexs) # express in percentages
      
    # if (model=='PolyFit' or model=='all'):  # Numpy Polyfit
      # pass
    
    return 
  

    

# %%
