import torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator
from sklearn.model_selection import train_test_split
from collections import defaultdict
import random


BATCH_SIZE = 32
FP_SIZE = 2048
RNG_SEED = 1234


def morgan_fp(mol):
    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=FP_SIZE)
    return morgan.GetFingerprint(mol)


def get_edge_indices(mol):
    edge_indices = []

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        edge_indices.append([i, j])
        edge_indices.append([j, i])

    return torch.tensor(edge_indices, dtype=torch.long).t()

'''
def random_split(dataset, seed=RNG_SEED):
    X = dataset.drop(['standard_value'], axis=1)
    y = dataset['standard_value']

    generator = torch.Generator().manual_seed(seed)
    X_train, X_test, X_validate = \
        torch.utils.data.random_split(X, [0.8, 0.1, 0.1], generator=generator)

    generator = torch.Generator().manual_seed(seed)
    y_train, y_test, y_validate = \
        torch.utils.data.random_split(y, [0.8, 0.1, 0.1], generator=generator)

    return X_train, X_test, X_validate, y_train, y_test, y_validate
'''

'''
def random_split(dataset, frac_train=0.8, frac_valid=0.1, frac_test=0.1, seed=RNG_SEED):
    generator = torch.Generator().manual_seed(seed)
    train, test, valid = torch.utils.data.random_split(dataset, [frac_train, frac_valid, frac_test], generator=generator)

    return train, test, valid
'''


def random_split(dataset, frac_train=0.8, frac_valid=0.1, frac_test=0.1, seed=RNG_SEED):
    test, train = train_test_split(dataset, test_size=frac_train, random_state=seed)
    test, valid = train_test_split(test, test_size=frac_test/(frac_valid+frac_test), random_state=seed)

    return train, test, valid


def get_scaffold(mol, include_chirality=False):
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold, isomericSmiles=include_chirality)


def compute_scaffold(mol):
    if mol is None:
        return None
    scaf = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaf)


def scaffold_split_df(df, smiles_col="canonical_smiles", frac_train=0.8, frac_valid=0.1, frac_test=0.1):
    df_copy = df.copy()
    df_copy["scaffold"] = df[smiles_col].apply(compute_scaffold)

    scaffold_groups = (
        df_copy.groupby("scaffold")
        .indices.values()
    )

    scaffold_groups = sorted(scaffold_groups, key=len, reverse=True)

    train_idx, valid_idx, test_idx = [], [], []
    n_total = len(df_copy)

    for group in scaffold_groups:
        if len(train_idx) / n_total < frac_train:
            train_idx.extend(group)
        elif len(valid_idx) / n_total < frac_valid:
            valid_idx.extend(group)
        else:
            test_idx.extend(group)

    return df.loc[train_idx], df.loc[test_idx], df.loc[valid_idx]
