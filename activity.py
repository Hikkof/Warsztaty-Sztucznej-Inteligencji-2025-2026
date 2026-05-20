# jeśli ostatecznie będziemy przyjmować tylko 'SMILES' i z niego wydobywać dodatkowe informacje to ten plik jest zbędny

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


### wiem że mogę użyć df_activity[[cecy_które_chcę]], ale jest tutaj wytłumaczenie dlaczego biorę konkretne cechy
df = pd.read_csv('initial cleaning/activity.csv', index_col=False)

# wyrzucam kolumny które (prawie) na pewno nie będą przydatne
# zamienianie IC50 na pIC50
df = df[df['standard_value'] > 0]
df.loc[(df['standard_type'] == 'IC50') & (df['standard_units'] == 'uM'), 'standard_value'] = -np.log10(df['standard_value'] / 1e6)
df.loc[(df['standard_type'] == 'IC50') & (df['standard_units'] == 'nM'), 'standard_value'] = -np.log10(df['standard_value'] / 1e9)

#df['pIC50 '] = np.log(df['IC50'] * 0.001) * (-1)
#df.loc[df['standard_type'] == 'IC50', 'standard_value'] = 10**(-1*df['standard_value'])
# te informacje już nie są potrzebne
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

# nadal nie wygląda prawidłow, ale nie wiem jakie wartości 'standard_value' powinno przyjmować
plt.boxplot(df['standard_value'])
plt.show()

print(df.info())

df.to_csv('activity.csv', index=False)

'''
 0   activity_id                10000 non-null  int64  
 1   assay_chembl_id            10000 non-null  int64  
 2   assay_type                 10000 non-null  object 
 3   bao_endpoint               10000 non-null  int64  
 4   bao_format                 10000 non-null  int64  
 5   bao_label                  10000 non-null  object 
 6   canonical_smiles           10000 non-null  object 
 7   molecule_chembl_id         10000 non-null  int64  
 8   parent_molecule_chembl_id  10000 non-null  int64  
 9   standard_value             10000 non-null  float64
 10  target_chembl_id           10000 non-null  int64  
 '''
