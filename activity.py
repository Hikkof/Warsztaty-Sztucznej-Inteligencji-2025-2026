import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


df = pd.read_csv('initial cleaning/activity.csv', index_col=False)

# zamienianie IC50 na pIC50
df = df[df['standard_value'] > 0]
df.loc[(df['standard_type'] == 'IC50') & (df['standard_units'] == 'uM'), 'standard_value'] = -np.log10(df['standard_value'] / 1e6)
df.loc[(df['standard_type'] == 'IC50') & (df['standard_units'] == 'nM'), 'standard_value'] = -np.log10(df['standard_value'] / 1e9)

df = df[['canonical_smiles', 'standard_value']]

df.drop_duplicates(subset=["canonical_smiles"], inplace=True)

plt.boxplot(df['standard_value'])
plt.show()

print(df.info())

df.to_csv('activity.csv', index=False)

