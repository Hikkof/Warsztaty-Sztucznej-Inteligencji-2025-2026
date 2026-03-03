import matplotlib.pyplot as plt
import pandas as pd

from chembl_webresource_client.new_client import new_client


target = new_client.target

records_num = len(target)
field_names = list(target.all()[0])

# pola które z jakiegoś powodu nie pozwalają na filtrowanie
fields_n = ['cross_references']
for field in fields_n:
    field_names.remove(field)

null_counts = [len(target.filter(**{name + '__isnull': True})) / records_num for name in field_names]

print(null_counts)
plt.barh(field_names, null_counts)
plt.xlabel('Percentage of nulls')
plt.ylabel('Field names')
plt.show()

target = target.filter(target_components__isnull=False)

df_target = pd.DataFrame(target[:10000])

df_target.drop([field_names[field_num] for field_num in range(len(field_names)) if null_counts[field_num] >= 0.50],
                 axis=1, inplace=True)

ids = ['target_chembl_id']
for i in ids:
    df_target[i] = df_target[i].str.extract('(\d+)', expand=False).astype('Int32')

df_target.to_csv('target.csv', index=False)

print(df_target.info())

'''
 0   cross_references    10000 non-null  object 
 1   organism            9999 non-null   object 
 2   pref_name           10000 non-null  object 
 3   species_group_flag  10000 non-null  bool   
 4   target_chembl_id    10000 non-null  Int32  
 5   target_components   10000 non-null  object 
 6   target_type         10000 non-null  object 
 7   tax_id              9999 non-null   float64
'''
