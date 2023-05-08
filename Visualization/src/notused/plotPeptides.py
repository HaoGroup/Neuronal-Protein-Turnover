# %%

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

# %%

df_PG_raw = pd.read_excel("../data/iMN_Peptide_Dataset.xlsx") # use relative paths, non-user-specific
# print(df.head().to_string())
df_PG = df_PG_raw[["PG.ProteinGroups","PG.Genes","Peptide","iMN_Day0_Light_Relative_Abundance","iMN_Day1_Light_Relative_Abundance","iMN_Day2_Light_Relative_Abundance","iMN_Day4_Light_Relative_Abundance","iMN_Day6_Light_Relative_Abundance"]]
# rename columns 
df_PG.columns = ["PG.ProteinGroups","PG.Genes","Peptide",0,1,2,4,6]
# print(df_PG.head().to_string())

PgList = df_PG_raw["PG.ProteinGroups"].unique() # ndarray (97,)


#%%
# Create function to plot peptide relative abundance
# E Lo 20230418
#

def abundancePlot1PgGroup(df, prtnGrp, maxNAcnt = 0, ylabel="Relative Abundance", xlabel='Day', charttitle='', savepng=False, saveFolder = ''):
    """
    Args:
        df (pandas dataframe): the raw dataframe 
        prtnGrp (str): Which ProteinGrp is being processed
        maxNAcnt (int, optional): maximum missing/NA data points allowed to include in plot. Defaults to 0.
        ylabel (str, optional): y-axis label
        xlabel (str, optional): x-axis label. Defaults to 'day'
        charttitle (str, optional): chart title. Defaults to the name of PG.ProteinGroup 
        savepng (boolean, optioal): save chart to png
        saveFolder (str, optional): (relative) path to save png
    return: plot object
    """
    df_1prtnGrp = df[ df["PG.ProteinGroups"] == prtnGrp ]  # there is still the PG.ProteinGroup column with single value prtnGrp
    # peptides = df_1prtnGrp["Peptide"].unique() # for prtnGrpId = -1, peptides.shape = (19,)
    
    # set the plot title
    if charttitle=='':
      # But first try find PG.Gene value for this prtnGrp:
      PgGene = df_1prtnGrp["PG.Genes"].unique() # should be length 1, or possibly zero
      # Now set chart title:
      charttitle = 'ProteinGroup = '+prtnGrp
      if len(PgGene) : charttitle += ' (Gene = '+ ','.join(PgGene.tolist()) +')'
      charttitle.replace(";","_") # replaces semi-colons, and potentially other special charcters with underscores.
      if maxNAcnt>0: charttitle += " "+str(maxNAcnt)+ ' NAs allowed'
    
    # only keep rows with at most maxNAcnt nan values
    rows2plot_bool = df_1prtnGrp.iloc[:,-5:].isna().sum(axis=1) <= maxNAcnt
    df_1prtnGrp = df_1prtnGrp[rows2plot_bool]
    
    if len(df_1prtnGrp) < 1 : return # stop if nothing to plot

    # If there are rows to plot, we now need column names 0,1,2,4,6 as index, the peptide name as new column head.
    # Complete in 2 steps.
    
    # Step 1: get column heads 0, 1, 2, 4, 6 as day values
    df_1prtnGrp_melted = pd.melt(df_1prtnGrp, id_vars=['PG.ProteinGroups', 'Peptide'], value_vars=df_1prtnGrp.columns[2:], var_name='day', value_name='rel_abundance')
    # print(df_1prtnGrp_melted.head())
    
    # Step 2: pivot back to have peptides as column names, values still rel_abundance
    df_1prtnGrp_pivot = df_1prtnGrp_melted.pivot(index='day', columns='Peptide' , values='rel_abundance')
    # print(df_1prtnGrp_pivot.head())

    # Plot all the peptides for this prtnGrp as line plot
    ax = df_1prtnGrp_pivot.plot()
    ax.set_title(charttitle)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if (savepng): 
      #set path
      import os
      filepath = os.path.join(saveFolder,"relAbundance_prtnGrp-"+ prtnGrp.replace(";","-") +"_"+str(maxNAcnt)+"NAs.png")
      plt.savefig(filepath, facecolor=(.65, .75, .75))
    
    return ax
#%%
groupId = -1
# Example: plot series without any nan
# abundancePlot1PgGroup(df_PG, PgList[groupId], savepng=True, saveFolder='media/plots')
abundancePlot1PgGroup(df_PG, PgList[groupId])
# Example: plot series at most one nan
abundancePlot1PgGroup(df_PG, PgList[groupId], maxNAcnt = 1)

#%%
# Try looping through PG.Genes
# Track all results with different maxNAcnt values in range(4) - 0,1,2,4,6 
# there can be 0, 1, 2, or 3 values missing for day 1, 2, 4, and 6.
# 
# list comprehension inside list comprehension
savepngnow = True
# os.getcwd() # in src folder
peptidePlots = [ [abundancePlot1PgGroup(df_PG,PgList[groupid], maxNAcnt=mxNA, savepng=savepngnow, saveFolder='../media/plots') for groupid in range(len(PgList))] for mxNA in range(4) ]



# %%
