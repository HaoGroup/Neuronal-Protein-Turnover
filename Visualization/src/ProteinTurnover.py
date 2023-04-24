#%%
import matplotlib.pyplot as plt
import pandas as pd
import os
# import seaborn as sns
# import numpy as np

#%%
# Create class ProteinTurnover to handle and organize the various parts of the study
#
class ProteinTurnover:
  """
  Neuronal Protein Turnover Study
  Ingest data
  Process various computations and data fit
  Creates plots
  """
  def __init__(self, filepath, x_values = [0,1,2,4,6], x_unit = "day") -> None:
    """_summary_

    Args:
        filepath (_type_): the source file location. Use os to ensure right format
        x_values (list, optional): the time marks on x-axis. Defaults to [0,1,2,4,6].
        x_unit (str, optional): the unit for x-axis. Defaults to "day".
    """
    self.__yValueName = 'Rel Abundance'
    self.__compoIndex = ['PG.ProteinGroups', 'Peptide', 'wtype'] # basic composite index used in DFs
    self.df_PgRaw = None # initialize, current structure has columns: ['PG.ProteinGroups', 'PG.Genes', 'PG.ProteinDescriptions', 'Peptide', 'iMN_Day{x}_{Light/Heavy}_Relative_Abundance', 'k_results', 'Protein_Turnover_from_k', 'residuals']
    self.df_Pg = None # initialize, cleaned and re-structured df for analysis ['PG.ProteinGroups', 'PG.Genes', 'PG.ProteinDescriptions', 'Peptide', 'wtype', 0, 1, 2, 4, 6, 'k_results', 'Protein_Turnover_from_k', 'residuals']
    self.df_PgGeneDescLookup = None # initialize, serve as lookup table between PgProteinGroups, PgGenes, and PgProteinDescriptions
    self.PgList = [] # initialize, the list of ProteinGroups in the dataset
    self.df_PgPivot = [] # initialize, for plotting
    self.__ingestData(filepath=filepath, x_values=x_values, x_unit=x_unit) # set df_PgRaw and df_PG
    self.__setPgGeneDescLookup() # set df_PgGeneDescLookup as well as PgList
    self.__setPgPivot()
    return
  
  def __ingestData(self, filepath, x_values = [], x_unit = "Day"):
    """
    Args:
        filepath (_type_): the source file location. Use os to ensure right format
        x_values (list, optional): the time marks on x-axis. Defaults to [0,1,2,4,6].
        x_unit (str, optional): the unit for x-axis. Defaults to "Day".

    Returns: None, will set df_PgRaw and df_Pg
    """
    if (filepath[-5:] == '.xlsx' or filepath[-4:] == '.xls'): 
      self.df_PgRaw = pd.read_excel(filepath)
    elif (filepath[-4:] == '.csv'):
      self.df_PgRaw = pd.read_csv(filepath)
    else:
      pass
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
    # rename columns as 0, 1, 2, 4, 6 for days
    import re
    df_light.rename(columns=lambda x: re.sub(r'iMN_Day(\d).*_Relative_Abundance$', r'\1', x), inplace = True)
    df_heavy.rename(columns=lambda x: re.sub(r'iMN_Day(\d).*_Relative_Abundance$', r'\1', x), inplace = True)
    # set composite index 
    df_light.set_index(self.__compoIndex, inplace=True)  
    df_heavy.set_index(self.__compoIndex, inplace=True)  
    # stack the two df together
    self.df_Pg = pd.concat([df_light, df_heavy])
    # print(self.df_Pg.head())
    # self.chkDfPgIndexUnique()
    
    return None
  
  def __setPgGeneDescLookup(self):
    """
    setting up the df_PgGeneDescLookup as lookup table self.df_PgGeneDescLookup as well as 
    setting up the list of ProteinGroup values as attribute self.PgList
    Returns: None
    """
    # set PgList
    self.PgList = self.df_PgRaw["PG.ProteinGroups"].unique()
    # set df_PgGeneDescLookup
    compoIndex = ['PG.ProteinGroups', 'PG.Genes']
    selectedCols = compoIndex + ['PG.ProteinDescriptions']
    df = self.df_PgRaw[selectedCols]
    df.set_index(compoIndex, inplace=True)    
    self.df_PgGeneDescLookup = df.drop_duplicates()
    # self.chkDfPgGeneLookupIndexUnique()
    
    return
  
  def __setPgPivot(self):
    """
    set up pivot table from df_PG table for analysis and plots
    Returns: None
    """
    xAxisName = 'Day'
    # select only rel_Abundance columns
    df = self.df_Pg.loc[:, self.df_Pg.columns.str.contains('^\d$')] # only columns with single digit as colnames
    df = df.stack(dropna = False) # melt method does not keep the original index keys. df is now a pd series
    df.index.names = df.index.names[:-1] + [xAxisName] # name new index column
    df.name = self.__yValueName # set name of the series (column)
    df = pd.DataFrame(df).reset_index() # from pd series to dataframe, and reset index to use pivot functions
    # df = pd.melt()
    self.df_PgPivot = df.pivot(index=xAxisName, columns = self.__compoIndex, values = self.__yValueName)
    return
  
  def chkDfPgIndexUnique(self): return print(f'df_Pg composite index ({", ".join( self.df_Pg.index.names )}) is unique? {self.df_Pg.index.is_unique}') 
  
  def chkDfPgGeneLookupIndexUnique(self): return print(f'df_PgGeneDescLookup composite index ({", ".join( self.df_PgGeneDescLookup.index.names )}) is unique? {self.df_PgGeneDescLookup.index.is_unique}') 
  
  def __set1PgChartTitle(self, prtnGrp, maxNAcnt = 0):
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
    
    # Add if maxNAcnt >0    
    # if maxNAcnt>0: title += " ("+str(maxNAcnt)+ ' NAs allowed)'
    return title
  
  def abundancePlot1Pg(self, prtnGrp, maxNAcnt = 0, xlabel=None, ylabel=None, charttitle=None, savepng=False, saveFolder = ''):
    """
    Args:
        prtnGrp (str): Which ProteinGrp is being processed
        maxNAcnt (int, optional): maximum missing/NA data points allowed to include in plot. Defaults to 0.
        ylabel (str, optional): y-axis label
        xlabel (str, optional): x-axis label. Defaults to pivot table index name
        charttitle (str, optional): chart title. Defaults to the name of PG.ProteinGroup 
        savepng (boolean, optioal): save chart to png
        saveFolder (str, optional): (relative) path to save png
    return: plot object
    """
    df = self.df_PgPivot[[prtnGrp]] # filter only one prtnGrp, with both light and heavy data, can have multiple peptides
    xlabel = xlabel if xlabel else self.df_PgPivot.index.name
    ylabel = ylabel if ylabel else self.__yValueName
    charttitle = charttitle if charttitle else self.__set1PgChartTitle(prtnGrp, maxNAcnt)

    # only keep rows with at most maxNAcnt nan values
    df = df.drop(columns=[col for col in df if df[col].isna().sum() >maxNAcnt])
    
    if len(df.columns) == 0 : return # nothing to plot
    
    # Now plot
    # plt.figure().set_figwidth(15)
    ax = df.plot()
    handles, labels = ax.get_legend_handles_labels()  
    lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(0,-0.15))
    natext = "("+str(maxNAcnt)+ '-NAs/series)' if maxNAcnt>1 else '(1-NA/series)' if maxNAcnt>0 else ' '
    axtext = ax.text(-0.15, 1.1, natext, transform=ax.transAxes)
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
  
  def abundancePlotPgAll(self, xlabel=None, ylabel=None, charttitle=None, savepng=False, saveFolder = ''):
    """
    Args:
        ylabel (str, optional): y-axis label
        xlabel (str, optional): x-axis label. Defaults to pivot table index name
        charttitle (str, optional): chart title. Defaults to the name of PG.ProteinGroup 
        savepng (boolean, optioal): save chart to png
        saveFolder (str, optional): (relative) path to save png
    return: None
    """
    for prtnGrp in self.PgList:
      for m in range( len(self.df_PgPivot.index) - 1 ):
        self.abundancePlot1Pg(prtnGrp=prtnGrp, maxNAcnt=m, savepng=savepng, saveFolder=saveFolder)
    return

#%%
file = os.path.join(os.getcwd(), "../data/iMN_Peptide_Dataset.xlsx") # assuming cwd is .../Visualization/src/ folder
pto = ProteinTurnover(filepath= file)
# pto.chkDfPgIndexUnique()
# pto.chkDfPgGeneLookupIndexUnique()
# savenow = True
savenow = False
savePath = "../media/plots/"
pto.abundancePlotPgAll(savepng=savenow, saveFolder=savePath)


# %%
