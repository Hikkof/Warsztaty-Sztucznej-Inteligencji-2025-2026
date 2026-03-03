import matplotlib.pyplot as plt
import pandas as pd

from chembl_webresource_client.new_client import new_client


molecule = new_client.molecule

records_num = len(molecule)

field_names = list(molecule.all()[0])
null_counts = [len(molecule.filter(**{name + '__isnull': True}))/records_num for name in field_names]

print(null_counts)
plt.barh(field_names, null_counts)
plt.xlabel('Percentage of nulls')
plt.ylabel('Field names')
plt.show()

molecule = molecule.filter(molecule_properties__isnull=False)
df_molecule = pd.DataFrame(molecule[:10000])

df_molecule.drop([field_names[field_num] for field_num in range(len(field_names)) if null_counts[field_num] >= 0.50], axis=1, inplace=True)

df_molecule = pd.concat([df_molecule, pd.json_normalize(df_molecule['molecule_properties'])], axis=1)
df_molecule.drop('molecule_properties', axis=1, inplace=True)

ids = ['molecule_chembl_id']
for i in ids:
    df_molecule[i] = df_molecule[i].str.extract('(\d+)', expand=False).astype('Int32')

# mam wrażenie że 'inorganic_flag' powinno mieć wartości True lub False, ale występują trzy wartości
df_molecule['polymer_flag'] = df_molecule['polymer_flag'].astype('bool')
df_molecule['aromatic_rings'] = df_molecule['aromatic_rings'].astype('Int32')
df_molecule['hba'] = df_molecule['hba'].astype('Int32')
df_molecule['hbd'] = df_molecule['hbd'].astype('Int32')
df_molecule['heavy_atoms'] = df_molecule['heavy_atoms'].astype('Int32')
df_molecule['rtb'] = df_molecule['rtb'].astype('Int32')

print(df_molecule.info())

df_molecule.to_csv('molecule.csv', index=False)


'''
 0   availability_type    9842 non-null   float64
 1   black_box_warning    10000 non-null  int64  
 2   chemical_probe       10000 non-null  int64  
 3   chirality            10000 non-null  int64  
 4   dosed_ingredient     10000 non-null  bool   
 5   first_in_class       10000 non-null  int64  
 6   inorganic_flag       10000 non-null  int64  
 7   molecule_chembl_id   10000 non-null  Int32  
 8   molecule_hierarchy   9860 non-null   object 
 9   molecule_structures  9959 non-null   object 
 10  molecule_type        10000 non-null  object 
 11  natural_product      10000 non-null  int64  
 12  oral                 10000 non-null  bool   
 13  orphan               10000 non-null  int64  
 14  parenteral           10000 non-null  bool   
 15  polymer_flag         10000 non-null  bool   
 16  prodrug              10000 non-null  int64  
 17  structure_type       10000 non-null  object 
 18  therapeutic_flag     10000 non-null  bool   
 19  topical              10000 non-null  bool   
 20  veterinary           10000 non-null  int64  
 21  withdrawn_flag       10000 non-null  bool   
 22  alogp                9820 non-null   object 
 23  aromatic_rings       9820 non-null   Int32  
 24  full_molformula      10000 non-null  object 
 25  full_mwt             10000 non-null  object 
 26  hba                  9820 non-null   Int32  
 27  hbd                  9820 non-null   Int32  
 28  heavy_atoms          9820 non-null   Int32  
 29  mw_freebase          10000 non-null  object 
 30  np_likeness_score    9820 non-null   object 
 31  num_ro5_violations   9820 non-null   float64
 32  psa                  9820 non-null   object 
 33  qed_weighted         9820 non-null   object 
 34  ro3_pass             9820 non-null   object 
 35  rtb                  9820 non-null   Int32
'''
