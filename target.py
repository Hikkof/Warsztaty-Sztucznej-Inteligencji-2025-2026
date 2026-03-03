import pandas as pd


df = pd.read_csv('initial cleaning/target.csv', index_col=False)

# 'cross_references' nie wygląda jak ważna informacja
df.drop('cross_references', axis=1, inplace=True)
# używam 'target_chembl_id', 'organism' i 'pref_name' nie jest potrzebne
df.drop('organism', axis=1, inplace=True)
df.drop('pref_name', axis=1, inplace=True)
# 'species_group_flag' nie wygląda jak ważna informacja
df.drop('species_group_flag', axis=1, inplace=True)
# kolumna "target_components" możliwe że zawiera przydatne informacje, ale zwykły
# pd.json_normalize(df['target_components']) nie potrafi rozbić tej kolumny
df.drop('target_components', axis=1, inplace=True)
# używam 'target_chembl_id', 'tax_id' nie jest potrzebne
df.drop('tax_id', axis=1, inplace=True)

print(df.info())

df.to_csv('target.csv', index=False)

'''
 0   target_chembl_id  10000 non-null  int64 
 1   target_type       10000 non-null  object
 '''