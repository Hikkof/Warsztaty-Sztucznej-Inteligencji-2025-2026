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

ATOMIC_SYMBOLS = ["H", "B", "C", "N", "O", "F", "Mg", "Si", "P", "S", "Cl", "Cu", "Zn", "Se", "Br", "Sn", "I"]
ATOMIC_VOCAB = {a: i for i, a in enumerate(ATOMIC_SYMBOLS)}
DEGREE_VOCAB = range(7)
HYBRIDIZATION_VOCAB = [
    Chem.rdchem.HybridizationType.SP,
    Chem.rdchem.HybridizationType.SP2,
    Chem.rdchem.HybridizationType.SP3,
    Chem.rdchem.HybridizationType.SP3D,
    Chem.rdchem.HybridizationType.SP3D2,
    Chem.rdchem.HybridizationType.UNSPECIFIED
]
BOND_TYPE_VOCAB = [
    Chem.rdchem.BondType.SINGLE,
    Chem.rdchem.BondType.DOUBLE,
    Chem.rdchem.BondType.TRIPLE,
    Chem.rdchem.BondType.AROMATIC
]
STEREO_VOCAB = [
    Chem.rdchem.BondStereo.STEREONONE,
    Chem.rdchem.BondStereo.STEREOANY,
    Chem.rdchem.BondStereo.STEREOZ,
    Chem.rdchem.BondStereo.STEREOE,
    Chem.rdchem.BondStereo.STEREOCIS,
    Chem.rdchem.BondStereo.STEREOTRANS
]

def one_hot(x, vocab):
    return [1 if x == val else 0 for val in vocab]


def morgan_fp(mol):
    morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=FP_SIZE)
    return morgan.GetFingerprint(mol)


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


def scaffold_split_df(df, mol_col='mol', frac_train=0.8, frac_valid=0.1):
    df_copy = df.copy()
    df_copy["scaffold"] = df[mol_col].apply(compute_scaffold)

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
