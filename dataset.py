import pandas as pd

from rdkit import Chem


# activity
df_activity = pd.read_csv('activity.csv', index_col=False)
# biorę kolumny które mnie interesują
df_activity = df_activity[['assay_type', 'canonical_smiles', 'molecule_chembl_id', 'standard_type', 'standard_value',
                           'target_chembl_id', 'units']]
# zamieniam 'canonical_smiles' na graf
df_activity['canonical_smiles'] = df_activity['canonical_smiles'].apply(Chem.MolFromSmiles)

# molecule
df_molecule = pd.read_csv('molecule.csv', index_col=False)
# biorę kolumny które mnie interesują
df_molecule = df_molecule[['chemical_probe', 'chirality', 'inorganic_flag', 'molecule_chembl_id', 'molecule_type',
                           'polymer_flag', 'alogp', 'aromatic_rings', 'full_mwt', 'hba', 'hbd', 'heavy_atoms',
                           'np_likeness_score', 'psa', 'qed_weighted', 'ro3_pass', 'rtb']]

# target
# tutaj trudno cokolwiek wyrzucić
df_target = pd.read_csv('target.csv', index_col=False)
df = df_activity.join(df_target.set_index("target_chembl_id"), on='target_chembl_id', how='inner')
df = df.join(df_molecule.set_index("molecule_chembl_id"), on='target_chembl_id', how='inner')

df.drop('target_chembl_id', axis=1, inplace=True)
df.drop('molecule_chembl_id', axis=1, inplace=True)

print(df.info())
