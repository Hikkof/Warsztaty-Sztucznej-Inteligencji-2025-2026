import pandas as pd


df = pd.read_csv('initial cleaning/molecule.csv', index_col=False)

# 'availability_type', black_box_warning, 'dosed_ingredient', 'first_in_class', 'molecule_hierarchy'
# nie wyglądają jak ważne informacje
df.drop('availability_type', axis=1, inplace=True)
df.drop('black_box_warning', axis=1, inplace=True)
df.drop('dosed_ingredient', axis=1, inplace=True)
df.drop('first_in_class', axis=1, inplace=True)
df.drop('molecule_hierarchy', axis=1, inplace=True)
# używam 'canonical_smiles', 'molecule_structures' nie jest potrzebne
df.drop('molecule_structures', axis=1, inplace=True)
# nie wyglądają jak ważne informacje
df.drop('oral', axis=1, inplace=True)
df.drop('orphan', axis=1, inplace=True)
df.drop('parenteral', axis=1, inplace=True)
df.drop('therapeutic_flag', axis=1, inplace=True)
df.drop('topical', axis=1, inplace=True)
df.drop('veterinary', axis=1, inplace=True)
df.drop('withdrawn_flag', axis=1, inplace=True)
# używam 'canonical_smiles', 'full_molformula' nie jest potrzebna
df.drop('full_molformula', axis=1, inplace=True)
# chyba nie są to ważne informacje
df.drop('mw_freebase', axis=1, inplace=True)
df.drop('num_ro5_violations', axis=1, inplace=True)

df.dropna(subset=['aromatic_rings'], inplace=True)

print(df.info())

df.to_csv('molecule.csv', index=False)

'''
 0   chemical_probe      9820 non-null   int64  
 1   chirality           9820 non-null   int64  
 2   inorganic_flag      9820 non-null   int64  
 3   molecule_chembl_id  9820 non-null   int64  
 4   molecule_type       9820 non-null   object 
 5   natural_product     9820 non-null   int64  
 6   polymer_flag        9820 non-null   bool   
 7   prodrug             9820 non-null   int64  
 8   structure_type      9820 non-null   object 
 9   alogp               9820 non-null   float64
 10  aromatic_rings      9820 non-null   float64
 11  full_mwt            9820 non-null   float64
 12  hba                 9820 non-null   float64
 13  hbd                 9820 non-null   float64
 14  heavy_atoms         9820 non-null   float64
 15  np_likeness_score   9820 non-null   float64
 16  psa                 9820 non-null   float64
 17  qed_weighted        9820 non-null   float64
 18  ro3_pass            9820 non-null   object 
 19  rtb                 9820 non-null   float64
'''