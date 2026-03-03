import pandas as pd


df = pd.read_csv('initial cleaning/activity.csv', index_col=False)

# wyrzucam kolumny które (prawie) na pewno nie będą przydatne
df['potential_duplicate'] = df['potential_duplicate'].astype('bool')
# nie sądzę że te informacje są przydatne
df.drop('document_chembl_id', axis=1, inplace=True)
df.drop('document_year', axis=1, inplace=True)
# jedyna wartość to 1
df.drop('src_id', axis=1, inplace=True)
# 'standard_flag' nie wygląda jak ważna informacja
df.drop('standard_flag', axis=1, inplace=True)
# używam 'relation', 'standard_relatio' nie jest potrzebne
df.drop('standard_relation', axis=1, inplace=True)
# 'record_id' nie wygląda jak ważna informacja
df.drop('record_id', axis=1, inplace=True)
# używam 'unints', 'qudt_units' i 'standard_units' nie jest potrzebne
df.drop('qudt_units', axis=1, inplace=True)
df.drop('standard_units', axis=1, inplace=True)
# 'target_organism' i 'target_tax_id' są wybrakowane, 'target_chembl_id' jest wystarczające
df.drop('target_organism', axis=1, inplace=True)
df.drop('target_tax_id', axis=1, inplace=True)
# używam 'standard_type', 'type' nie jest potrzebne
df.drop('type', axis=1, inplace=True)
# 'uo_units' nie wygląda jak ważna informacja
df.drop('uo_units', axis=1, inplace=True)
# używam 'standard_value', 'value' nie jest potrzebne
df.drop('value', axis=1, inplace=True)

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
 9   potential_duplicate        10000 non-null  bool   
 10  relation                   10000 non-null  object 
 11  standard_type              10000 non-null  object 
 12  standard_value             10000 non-null  float64
 13  target_chembl_id           10000 non-null  int64  
 14  units                      7911 non-null   object 
 '''
