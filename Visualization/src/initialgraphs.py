# %%

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

# %%

# df = pd.read_excel("/Users/medhaswetasen/Documents/GitHub/Neuronal-Protein-Turnover/Visualization/iMN_Peptide_Dataset.xlsx")
# df = pd.read_excel("/Users/medhaswetasen/Documents/GitHub/Neuronal-Protein-Turnover/Visualization/data/iMN_Peptide_Dataset.xlsx")
df = pd.read_excel("../data/iMN_Peptide_Dataset.xlsx") # use relative paths, non-user-specific
# print(df.head().to_string())
# %%
df_1 = df[["PG.Genes","Peptide","iMN_Day0_Light_Relative_Abundance","iMN_Day1_Light_Relative_Abundance","iMN_Day2_Light_Relative_Abundance","iMN_Day4_Light_Relative_Abundance","iMN_Day6_Light_Relative_Abundance"]]
df_2 = df[["PG.Genes","Peptide","iMN_Day0_Heavy_Relative_Abundance","iMN_Day1_Heavy_Relative_Abundance","iMN_Day2_Heavy_Relative_Abundance","iMN_Day4_Heavy_Relative_Abundance","iMN_Day6_Heavy_Relative_Abundance"]]

df_1.columns = ["PG.Genes","Peptide",0,1,2,4,6]
df_2.columns = ["PG.Genes","Peptide",0,1,2,4,6]

# print(df_1.head().to_string())

# gk_1 = pd.pivot_table(df_1, index=["PG.Genes","Peptide"])
# gk_2 = pd.pivot_table(df_2, index=["PG.Genes","Peptide"])

# Group the DataFrame by Region
# gk_1 = df_1.groupby(["PG.Genes","Peptide"])
# gk_2 = df_2.groupby(["PG.Genes","Peptide"])

unique_entries = df["PG.Genes"].unique()

# unique_entries_1 = df_1["PG.Genes"].unique()
# unique_entries_2 = df_2["PG.Genes"].unique()

# TO MODIFY THIS GENERALIZATION

for unique in unique_entries:
    df_3 = df_1[df_1["PG.Genes"] == unique ]
    df_4 = df_2[df_2["PG.Genes"] == unique ]
    pep = df_3["Peptide"].unique()
    df_melted_4 = pd.melt(df_4, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Score')
    f_melted_sorted_4 = df_melted_4.sort_values(by=['Peptide', 'Values'])
    df_melted_3 = pd.melt(df_3, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Score')
    f_melted_sorted_3 = df_melted_3.sort_values(by=['Peptide', 'Values'])
    for p in pep:
        df_33 = f_melted_sorted_3[f_melted_sorted_3["Peptide"] == p]
        df_44 = f_melted_sorted_4[f_melted_sorted_4["Peptide"] == p]
        df_33 = df_33 [["Peptide","Score"]]
        df_44 = df_44 [["Peptide", "Score"]]
        df_33.plot()
    plt.show()


# TRYING FOR SINGLE PROTEIN:

df_3 = df_1[df_1["PG.Genes"] == unique_entries[-1] ]
df_4 = df_2[df_2["PG.Genes"] == unique_entries[-1] ]
pep = df_3["Peptide"].unique()
df_melted_4 = pd.melt(df_4, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Score')
f_melted_sorted_4 = df_melted_4.sort_values(by=['Peptide', 'Values'])
f_melted_sorted_4 = f_melted_sorted_4.reset_index()
df_melted_3 = pd.melt(df_3, id_vars=['PG.Genes', 'Peptide'], value_vars=[0, 1, 2, 4, 6], var_name='Values', value_name='Score')
f_melted_sorted_3 = df_melted_3.sort_values(by=['Peptide', 'Values'])
f_melted_sorted_3 = f_melted_sorted_3.reset_index()
for p in pep:
    df_33 = f_melted_sorted_3[f_melted_sorted_3["Peptide"] == p]
    df_44 = f_melted_sorted_4[f_melted_sorted_4["Peptide"] == p]
    df_33 = df_33 [["Score"]]
    df_44 = df_44 [[ "Score"]]
    df_33 = df_33.reset_index()
    df_44 = df_44.reset_index()
    df_33.plot()
    df_44.plot()
plt.show()






