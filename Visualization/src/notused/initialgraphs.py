# %%

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

# %%

# df = pd.read_excel("/Users/medhaswetasen/Documents/GitHub/Neuronal-Protein-Turnover/Visualization/iMN_Peptide_Dataset.xlsx")
# df = pd.read_excel("/Users/medhaswetasen/Documents/GitHub/Neuronal-Protein-Turnover/Visualization/data/iMN_Peptide_Dataset.xlsx")
df_PGgenes = pd.read_excel("./data/iMN_Peptide_Dataset.xlsx") # use relative paths, non-user-specific
# print(df.head().to_string())
# %%
df_1 = df_PGgenes[["PG.Genes","Peptide","iMN_Day0_Light_Relative_Abundance","iMN_Day1_Light_Relative_Abundance","iMN_Day2_Light_Relative_Abundance","iMN_Day4_Light_Relative_Abundance","iMN_Day6_Light_Relative_Abundance"]]
# df_2 = df[["PG.Genes","Peptide","iMN_Day0_Heavy_Relative_Abundance","iMN_Day1_Heavy_Relative_Abundance","iMN_Day2_Heavy_Relative_Abundance","iMN_Day4_Heavy_Relative_Abundance","iMN_Day6_Heavy_Relative_Abundance"]]
df_2 = df_1.copy() 
#  ???????????????????? What is df_2 (and df_4, df_44 used for????)

df_1.columns = ["PG.Genes","Peptide",0,1,2,4,6]
# df_2.columns = ["PG.Genes","Peptide",0,1,2,4,6]
df_2.columns = df_1.columns

# print(df_1.head().to_string())

# gk_1 = pd.pivot_table(df_1, index=["PG.Genes","Peptide"])
# gk_2 = pd.pivot_table(df_2, index=["PG.Genes","Peptide"])

# Group the DataFrame by Region
# gk_1 = df_1.groupby(["PG.Genes","Peptide"])
# gk_2 = df_2.groupby(["PG.Genes","Peptide"])

#%%
PgGenes = df_PGgenes["PG.Genes"].unique() # ndarray (97,)

# PgGenes_1 = df_1["PG.Genes"].unique()
# PgGenes_2 = df_2["PG.Genes"].unique()

#%%
# ELO Testing Plot mult-series
# Test just one PgGene
# Need dataframe with different series for each column
geneId = -1 # int
PgGene = PgGenes[geneId] # the gene code in the DB
df_3 = df_1[df_1["PG.Genes"] == PgGene ]
# df_3.drop(['PG.Genes'], axis=1, inplace=True)
peptides = df_3["Peptide"].unique() # for geneEd = -1, peptides.shape = (19,)

# print(df_3.head())
# need column names 0,1,2,4,6 as index, the peptide name as new column head.
# In 2 steps.
# Step 1: get column heads 0, 1, 2, 4, 6 as time values
df_melted_3 = pd.melt(df_3, id_vars=['PG.Genes', 'Peptide'], value_vars=df_3.columns[1:], var_name='time', value_name='rel_abundance')
# df_melted_3.head()
# Step 2: pivot back to have peptides as column names, values still rel_abundance
df_pivot_3 = df_melted_3.pivot(index='time', columns='Peptide' , values='rel_abundance')
# print(df_pivot_3.head())

# Plot all the peptides for this gene as line plot
df_pivot_3.plot()

#%%
# Create function to plot

def abundancePlot1gene(gene, minDataPts = 5):
    """
    Args:
        gene (str): Which gene is being processed
        minDataPts (int, optional): minimum data points required to create plot. Defaults to 5.
    return: plot object
    """
    
    return




#%%
# TO MODIFY THIS GENERALIZATION

# for unique in PgGenes:
#     df_3 = df_1[df_1["PG.Genes"] == unique ]
#     df_4 = df_2[df_2["PG.Genes"] == unique ]
#     peptides = df_3["Peptide"].unique()
#     df_melted_4 = pd.melt(df_4, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Score')
#     f_melted_sorted_4 = df_melted_4.sort_values(by=['Peptide', 'Values'])
#     df_melted_3 = pd.melt(df_3, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Score')
#     f_melted_sorted_3 = df_melted_3.sort_values(by=['Peptide', 'Values'])
#     for p in pep:
#         df_33 = f_melted_sorted_3[f_melted_sorted_3["Peptide"] == p]
#         df_44 = f_melted_sorted_4[f_melted_sorted_4["Peptide"] == p]
#         df_33 = df_33 [["Peptide","Score"]]
#         df_44 = df_44 [["Peptide", "Score"]]
#         df_33.plot()
#     plt.show()




# #%%
# # TRYING FOR SINGLE PROTEIN:

# df_3 = df_1[df_1["PG.Genes"] == PgGenes[-1] ]
# df_4 = df_2[df_2["PG.Genes"] == PgGenes[-1] ]
# peptides = df_3["Peptide"].unique()
# df_melted_4 = pd.melt(df_4, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Score')
# f_melted_sorted_4 = df_melted_4.sort_values(by=['Peptide', 'Values'])
# f_melted_sorted_4 = f_melted_sorted_4.reset_index()
# df_melted_3 = pd.melt(df_3, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Score')
# f_melted_sorted_3 = df_melted_3.sort_values(by=['Peptide', 'Values'])
# f_melted_sorted_3 = f_melted_sorted_3.reset_index()
# for p in pep:
#     df_33 = f_melted_sorted_3[f_melted_sorted_3["Peptide"] == p]
#     df_44 = f_melted_sorted_4[f_melted_sorted_4["Peptide"] == p]
#     df_33 = df_33 [["Score"]]
#     df_44 = df_44 [[ "Score"]]
#     df_33 = df_33.reset_index()
#     df_44 = df_44.reset_index()
#     df_33.plot()
#     df_44.plot()
# plt.show()


# #%%
# # TRYING FOR THE PROTEIN IN THE PPT (USING JUST ONE PROTEIN AS BOTH DON'T EXIST): (TRY 1: USING PREVIOUS LOGIC)

# # df[df["PG.Genes"]=="SYPL1"]
# # df[df["PG.Genes"]=='CTSD']

# # unique = 'SYPL1'
# # unique = 'CTSD'

# # if unique in PgGenes:
# #    print("Element found!")
# # else:
# #    print("Element not found.")

# # df[df["PG.Genes"]=='FBLL1'] ~~ HAS NA

# #%%
# # RUN THIS FIRST
# nonna = []
# for unique in PgGenes:
#     df_3 = df_1[df_1["PG.Genes"] == unique]
#     df_4 = df_2[df_2["PG.Genes"] == unique]
#     if np.all(df_3.notnull()) and np.all(df_4.notnull()):
#         nonna.append(unique)


# for unique in nonna:
#     df_3 = df_1[df_1["PG.Genes"] == unique ]
#     df_4 = df_2[df_2["PG.Genes"] == unique ]
#     peptides = df_3["Peptide"].unique()
#     df_melted_4 = pd.melt(df_4, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Heavy Medium')
#     f_melted_sorted_4 = df_melted_4.sort_values(by=['Peptide', 'Values'])
#     f_melted_sorted_4 = f_melted_sorted_4.set_index("Values")
#     f_melted_sorted_4 = f_melted_sorted_4[["Heavy Medium"]]
#     df_melted_3 = pd.melt(df_3, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Light Medium')
#     f_melted_sorted_3 = df_melted_3.sort_values(by=['Peptide', 'Values'])
#     f_melted_sorted_3 = f_melted_sorted_3.set_index("Values")
#     merged_df = pd.merge(f_melted_sorted_3, f_melted_sorted_4, left_index=True, right_index=True)
#     for p in pep:
#         df_33 = merged_df[merged_df["Peptide"] == p]
#         df_33.plot(y=['Heavy Medium', 'Light Medium'])
#     plt.show()



# #%%
# # RUN THIS SECOND

# for unique in PgGenes:
#     df_3 = df_1[df_1["PG.Genes"] == unique ]
#     df_4 = df_2[df_2["PG.Genes"] == unique ]
#     peptides = df_3["Peptide"].unique()
#     df_melted_4 = pd.melt(df_4, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Heavy Medium')
#     f_melted_sorted_4 = df_melted_4.sort_values(by=['Peptide', 'Values'])
#     f_melted_sorted_4 = f_melted_sorted_4.set_index("Values")
#     f_melted_sorted_4 = f_melted_sorted_4[["Heavy Medium"]]
#     df_melted_3 = pd.melt(df_3, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Light Medium')
#     f_melted_sorted_3 = df_melted_3.sort_values(by=['Peptide', 'Values'])
#     f_melted_sorted_3 = f_melted_sorted_3.set_index("Values")
#     merged_df = pd.merge(f_melted_sorted_3, f_melted_sorted_4, left_index=True, right_index=True)
#     for p in pep:
#         df_33 = merged_df[merged_df["Peptide"] == p]
#         df_33.plot(y=['Heavy Medium', 'Light Medium'])
#     plt.show()


# df_3 = df_1[df_1["PG.Genes"] == PgGenes[-1] ]
# df_4 = df_2[df_2["PG.Genes"] == PgGenes[-1] ]
# peptides = df_3["Peptide"].unique()
# df_melted_4 = pd.melt(df_4, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Heavy Medium')
# f_melted_sorted_4 = df_melted_4.sort_values(by=['PG.Genes', 'Peptide', 'Values'])
# # f_melted_sorted_4 = f_melted_sorted_4.set_index("Values")
# f_melted_sorted_4 = f_melted_sorted_4.reset_index()
# f_melted_sorted_4 = f_melted_sorted_4[["Heavy Medium"]]
# df_melted_3 = pd.melt(df_3, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Light Medium')
# f_melted_sorted_3 = df_melted_3.sort_values(by=['PG.Genes', 'Peptide', 'Values'])
# # f_melted_sorted_3 = f_melted_sorted_3.set_index("Values")
# f_melted_sorted_3 = f_melted_sorted_3.reset_index()
# merged_df = pd.merge(f_melted_sorted_3, f_melted_sorted_4, left_index=True, right_index=True)
# merged_df = merged_df.sort_values(by=['Peptide', 'Values'])
# for p in pep:
#     df_33 = merged_df[merged_df["Peptide"] == p]
#     df_33.plot(x = 'Values', y=['Heavy Medium', 'Light Medium'])
# plt.show()










