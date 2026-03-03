import matplotlib.pyplot as plt
import pandas as pd

from chembl_webresource_client.new_client import new_client

activity = new_client.activity

records_num = len(activity)
field_names = list(activity.all()[0])

# pola które z jakiegoś powodu nie pozwalają na filtrowanie
fields_n = ['bao_format', 'bao_label', 'data_validity_description', 'src_id', 'toid']
for field in fields_n:
    field_names.remove(field)

null_counts = [len(activity.filter(**{name + '__isnull': True})) / records_num for name in field_names]

print(null_counts)
plt.barh(field_names, null_counts)
plt.xlabel('Percentage of nulls')
plt.ylabel('Field names')
plt.show()

activity = activity.filter(standard_relation__iexact="=").filter(canonical_smiles__isnull=False).filter(value__isnull=False)
activity = activity.filter(standard_type__iexact="IC50")

df_activity = pd.DataFrame(activity[:10000])

df_activity.drop([field_names[field_num] for field_num in range(len(field_names)) if null_counts[field_num] >= 0.50],
                 axis=1, inplace=True)

to_drop = ['activity_comment', 'assay_description', 'data_validity_description', 'target_pref_name', 'toid']
df_activity.drop(to_drop, axis=1, inplace=True)

ids = ['assay_chembl_id', 'bao_endpoint', 'bao_format', 'document_chembl_id', 'molecule_chembl_id',
       'parent_molecule_chembl_id', 'target_chembl_id', 'uo_units']
for i in ids:
    df_activity[i] = df_activity[i].str.extract('(\d+)', expand=False).astype('Int32')

df_activity.to_csv('activity.csv', index=False)

print(df_activity.info())

'''
 0   activity_id                10000 non-null  int64 
 1   assay_chembl_id            10000 non-null  Int32 
 2   assay_type                 10000 non-null  object
 3   bao_endpoint               10000 non-null  Int32 
 4   bao_format                 10000 non-null  Int32 
 5   bao_label                  10000 non-null  object
 6   canonical_smiles           10000 non-null  object
 7   document_chembl_id         10000 non-null  Int32 
 8   document_year              10000 non-null  int64 
 9   molecule_chembl_id         10000 non-null  Int32 
 10  parent_molecule_chembl_id  10000 non-null  Int32 
 11  potential_duplicate        10000 non-null  int64 
 12  qudt_units                 7395 non-null   object
 13  record_id                  10000 non-null  int64 
 14  relation                   10000 non-null  object
 15  src_id                     10000 non-null  int64 
 16  standard_flag              10000 non-null  int64 
 17  standard_relation          10000 non-null  object
 18  standard_type              10000 non-null  object
 19  standard_units             8376 non-null   object
 20  standard_value             10000 non-null  object
 21  target_chembl_id           10000 non-null  Int32 
 22  target_organism            8783 non-null   object
 23  target_tax_id              8783 non-null   object
 24  type                       10000 non-null  object
 25  units                      7911 non-null   object
 26  uo_units                   7809 non-null   Int32 
 27  value                      10000 non-null  object
'''
