# %%

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

# %%

# df = pd.read_excel("/Users/medhaswetasen/Documents/GitHub/Neuronal-Protein-Turnover/Visualization/data/iMN_Peptide_Dataset.xlsx")
df_PGgenes_raw = pd.read_excel("./data/iMN_Peptide_Dataset.xlsx") # use relative paths, non-user-specific
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

def abundancePlot1gene(df, gene, maxNAcnt = 0, ylabel="Relative Abundance", xlabel='Day', charttitle=''):
    """
    Args:
        df (pandas dataframe): the raw dataframe 
        gene (str): Which gene is being processed
        maxNAcnt (int, optional): maximum missing/NA data points allowed to include in plot. Defaults to 0.
        ylabel (str, optional): y-axis label
        xlabel (str, optional): x-axis label. Defaults to 'day'
        charttitle (str, optional): chart title. Defaults to the name of PG.Gene 
    return: plot object
    """
    charttitle = 'For PG.Gene = '+gene if charttitle=='' else charttitle
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
    
    return ax

geneId = -1
# Example: plot series without any nan
abundancePlot1gene(df_PGgenes, PgGenes[geneId])
# Example: plot series at most one nan
abundancePlot1gene(df_PGgenes, PgGenes[geneId], maxNAcnt = 1)





# %%
