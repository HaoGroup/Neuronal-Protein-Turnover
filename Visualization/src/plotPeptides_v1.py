# v1, group by PG.Gene
# %%

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

# %%

df_PGgenes_raw = pd.read_excel("../data/iMN_Peptide_Dataset.xlsx") # use relative paths, non-user-specific
# print(df.head().to_string())
df_PGgenes = df_PGgenes_raw[["PG.Genes","Peptide","iMN_Day0_Light_Relative_Abundance","iMN_Day1_Light_Relative_Abundance","iMN_Day2_Light_Relative_Abundance","iMN_Day4_Light_Relative_Abundance","iMN_Day6_Light_Relative_Abundance"]]
# rename columns 
df_PGgenes.columns = ["PG.Genes","Peptide",0,1,2,4,6]
# print(df_PGgenes.head().to_string())

PgGenes = df_PGgenes_raw["PG.Genes"].unique() # ndarray (97,)


#%%
# Create function to plot peptide relative abundance
# E Lo 20230418
#

def abundancePlot1gene(df, gene, maxNAcnt = 0, ylabel="Relative Abundance", xlabel='Day', charttitle='', savepng=False, saveFolder = ''):
    """
    Args:
        df (pandas dataframe): the raw dataframe 
        gene (str): Which gene is being processed
        maxNAcnt (int, optional): maximum missing/NA data points allowed to include in plot. Defaults to 0.
        ylabel (str, optional): y-axis label
        xlabel (str, optional): x-axis label. Defaults to 'day'
        charttitle (str, optional): chart title. Defaults to the name of PG.Gene 
        savepng (boolean, optioal): save chart to png
        saveFolder (str, optional): (relative) path to save png
    return: plot object
    """
    # set the plot title
    if charttitle=='':
      charttitle = 'For PG.Gene = '+gene 
      if maxNAcnt>0: charttitle += " ("+str(maxNAcnt)+ ' NAs allowed)'
    
    df_1gene = df[ df["PG.Genes"] == gene ]  # there is still the PG.Gene column with single value gene
    # peptides = df_1gene["Peptide"].unique() # for geneEd = -1, peptides.shape = (19,)
    
    # only keep rows with at most maxNAcnt nan values
    rows2plot_bool = df_1gene.iloc[:,-5:].isna().sum(axis=1) <= maxNAcnt
    df_1gene = df_1gene[rows2plot_bool]
    
    if len(df_1gene) < 1 : return # stop if nothing to plot

    # If there are rows to plot, we now need column names 0,1,2,4,6 as index, the peptide name as new column head.
    # Complete in 2 steps.
    
    # Step 1: get column heads 0, 1, 2, 4, 6 as day values
    df_1gene_melted = pd.melt(df_1gene, id_vars=['PG.Genes', 'Peptide'], value_vars=df_1gene.columns[1:], var_name='day', value_name='rel_abundance')
    # print(df_1gene_melted.head())
    
    # Step 2: pivot back to have peptides as column names, values still rel_abundance
    df_1gene_pivot = df_1gene_melted.pivot(index='day', columns='Peptide' , values='rel_abundance')
    # print(df_1gene_pivot.head())

    # Plot all the peptides for this gene as line plot
    ax = df_1gene_pivot.plot()
    ax.set_title(charttitle)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if (savepng): 
      #set path
      import os
      filepath = os.path.join(saveFolder,"relAbundance_gene-"+gene+"_"+str(maxNAcnt)+"NAs.png")
      plt.savefig(filepath, facecolor=(.65, .75, .75))
    
    return ax

geneId = -1
# Example: plot series without any nan
# abundancePlot1gene(df_PGgenes, PgGenes[geneId], savepng=True, saveFolder='media/plots')
abundancePlot1gene(df_PGgenes, PgGenes[geneId])
# Example: plot series at most one nan
abundancePlot1gene(df_PGgenes, PgGenes[geneId], maxNAcnt = 1)

#%%
# Try looping through PG.Genes
# Track all results with different maxNAcnt values in range(4) - 0,1,2,4,6 
# there can be 0, 1, 2, or 3 values missing for day 1, 2, 4, and 6.
# 
# list comprehension inside list comprehension
savepngnow = False
peptidePlots = [ [abundancePlot1gene(df_PGgenes,PgGenes[geneid], maxNAcnt=mxNA, savepng=savepngnow, saveFolder='media/plots') for geneid in range(len(PgGenes))] for mxNA in range(4) ]



# %%
