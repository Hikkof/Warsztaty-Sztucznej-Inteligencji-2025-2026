import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


### wiem że mogę użyć df_activity[[cecy_które_chcę]], ale jest tutaj wytłumaczenie dlaczego biorę konkretne cechy
df = pd.read_csv('initial cleaning/activity_203.csv', index_col=False)

# wyrzucam kolumny które (prawie) na pewno nie będą przydatne
# zamienianie IC50 na pIC50
df = df[df['standard_value'] > 0]
df.loc[(df['standard_type'] == 'IC50') & (df['standard_units'] == 'uM'), 'standard_value'] = -np.log10(df['standard_value'] / 1e6)
df.loc[(df['standard_type'] == 'IC50') & (df['standard_units'] == 'nM'), 'standard_value'] = -np.log10(df['standard_value'] / 1e9)
# todo: df = df[(df['standard_value'] < 12) & (df['standard_value'] > 3)]
#df = df[(df['standard_value'] < 12) & (df['standard_value'] > 3)]
#df['pIC50 '] = np.log(df['IC50'] * 0.001) * (-1)
#df.loc[df['standard_type'] == 'IC50', 'standard_value'] = 10**(-1*df['standard_value'])
# te informacje już nie są potrzebne
'''
df.drop('standard_units', axis=1, inplace=True)
df.drop('units', axis=1, inplace=True)
df.drop('uo_units', axis=1, inplace=True)
df.drop('qudt_units', axis=1, inplace=True)

df.drop('target_chembl_id', axis=1, inplace=True)
df.drop('target_organism', axis=1, inplace=True)
df.drop('target_pref_name', axis=1, inplace=True)
df.drop('target_tax_id', axis=1, inplace=True)

df.drop('standard_type', axis=1, inplace=True)
df.drop('type', axis=1, inplace=True)

# używam 'standard_value', 'value' nie jest potrzebne
df.drop('value', axis=1, inplace=True)
'''
df = df[['canonical_smiles', 'standard_value']]
df.drop_duplicates(inplace=True)
# todo: df['canonical_smiles'].duplicated()
df.drop_duplicates(subset=["canonical_smiles"], inplace=True)

plt.boxplot(df['standard_value'])
plt.show()

print(df.info())

df.to_csv('activity.csv', index=False)

'''
 0   canonical_smiles  17627 non-null  object 
 1   standard_value    17627 non-null  float64
 '''
