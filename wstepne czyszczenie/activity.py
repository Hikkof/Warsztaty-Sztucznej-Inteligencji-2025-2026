import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from chembl_webresource_client.new_client import new_client


activity = new_client.activity

records_num = len(activity)
field_names = list(activity.all()[0])

# pola które z jakiegoś powodu nie pozwalają na filtrowanie używając __isnull
fields_n = ['bao_format', 'bao_label', 'data_validity_description', 'src_id', 'toid', 'parent_molecule_chembl_id']
for field in fields_n:
    field_names.remove(field)

null_counts = [len(activity.filter(**{name + '__isnull': True})) / records_num for name in field_names]

print(null_counts)
plt.barh(field_names, null_counts)
plt.xlabel('Percentage of nulls')
plt.ylabel('Field names')
plt.show()

activity = activity.filter(canonical_smiles__isnull=False,
                           standard_relation='=',
                           standard_type__in=['IC50', 'pIC50'],
                           standard_units__in=['nM', 'uM'],
                           standard_value__isnull=False,
                           target_chembl_id='CHEMBL203') \
    .only(['activity_id', 'assay_chembl_id', 'assay_type', 'bao_endpoint', 'bao_format', 'bao_label',
           'canonical_smiles', 'molecule_chembl_id', 'qudt_units', 'standard_type', 'standard_units', 'standard_value',
           'target_chembl_id', 'target_organism', 'target_pref_name', 'target_tax_id', 'type', 'units', 'uo_units',
           'value'])

df_activity = pd.DataFrame(activity)

ids = ['assay_chembl_id', 'bao_endpoint', 'bao_format', 'molecule_chembl_id', 'target_chembl_id', 'uo_units']
for i in ids:
    df_activity[i] = df_activity[i].str.extract('(\d+)', expand=False).astype('Int32')

df_activity.to_csv('activity.csv', index=False)

print(df_activity.info())
