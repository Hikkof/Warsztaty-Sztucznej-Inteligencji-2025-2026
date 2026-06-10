from abc import ABC

import matplotlib.pyplot as plt
import pandas as pd
import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool
from torch_geometric.data import Data, InMemoryDataset
from torch_geometric.loader import DataLoader
from torchmetrics.regression import R2Score
from rdkit import Chem

import funkcje


class GNN(nn.Module):
    def __init__(self, in_channels):
        super(GNN, self).__init__()

        self.conv1 = GCNConv(in_channels, 256)
        self.conv2 = GCNConv(256, 128)
        self.conv3 = GCNConv(128, 64)
        self.conv4 = GCNConv(64, 64)
        self.lin1 = nn.Linear(64, 32)
        self.lin2 = nn.Linear(32, 1)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index).relu()
        x = self.conv3(x, edge_index).relu()
        x = self.conv4(x, edge_index).relu()
        x = global_mean_pool(x, batch)
        x = self.lin1(x).relu()
        x = self.lin2(x)

        return x.squeeze(1)


def mol_to_graph(mol, target):
    atom_features = []
    for atom in mol.GetAtoms():
        atom_features.append(
            funkcje.one_hot(atom.GetSymbol(), funkcje.ATOMIC_VOCAB) +
            funkcje.one_hot(atom.GetChiralTag(), funkcje.DEGREE_VOCAB) +
            funkcje.one_hot(atom.GetTotalDegree(), funkcje.DEGREE_VOCAB) +
            funkcje.one_hot(atom.GetHybridization(), funkcje.HYBRIDIZATION_VOCAB) + [
                atom.GetFormalCharge(),
                atom.GetTotalNumHs(),
                atom.GetMass(),
                atom.GetDegree(),
                atom.GetTotalNumHs(),
                atom.GetValence(which=Chem.ValenceType.EXPLICIT),
                atom.GetValence(which=Chem.ValenceType.IMPLICIT),
                atom.GetNumRadicalElectrons(),
                atom.GetIsotope(),
                atom.HasProp('_CIPCode'),
                atom.GetIsAromatic(),
                atom.IsInRing(),
            ])
    x = torch.tensor(atom_features, dtype=torch.float)

    edge_index = []
    edge_attr = []

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        edge_index.append([i, j])
        edge_index.append([j, i])

        bond_features = \
            funkcje.one_hot(bond.GetBondType(), funkcje.BOND_TYPE_VOCAB) + \
            funkcje.one_hot(bond.GetStereo(), funkcje.STEREO_VOCAB) + [
                bond.GetIsConjugated(),
                bond.IsInRing(),
                bond.GetIsAromatic()
            ]

        edge_attr.append(bond_features)
        edge_attr.append(bond_features)

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    y = torch.tensor([target], dtype=torch.float)

    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y)


class CustomDataset(InMemoryDataset, ABC):
    def __init__(self, data_frame):
        super().__init__()
        data_list = []

        for mol, target in zip(data_frame["mol"], data_frame["standard_value"]):
            graph = mol_to_graph(mol, target)
            data_list.append(graph)

        self.data, self.slices = self.collate(data_list)


df_activity = pd.read_csv('activity.csv', index_col=False)
df_activity = df_activity[['canonical_smiles', 'standard_value']]

df_activity["mol"] = df_activity["canonical_smiles"].apply(Chem.MolFromSmiles)
# rdkit wpisuje 'None' w niemożliwych cząsteczkach
df_activity = df_activity.dropna(subset=["mol"]).reset_index(drop=True)

epochs = 100
loss_function = nn.MSELoss()
r2 = R2Score()
split_types = {'Random': funkcje.random_split(df_activity),
               'Scaffold': funkcje.scaffold_split_df(df_activity, mol_col='mol')}

for split in split_types:
    print(f"Rodzaj podziału: {split}")
    train, test, valid = split_types[split]

    train_dataset = CustomDataset(train)
    train_loader = DataLoader(train_dataset, batch_size=funkcje.BATCH_SIZE, shuffle=True)

    valid_dataset = CustomDataset(valid)
    valid_loader = DataLoader(valid_dataset, batch_size=funkcje.BATCH_SIZE, shuffle=True)

    test_dataset = CustomDataset(test)
    test_loader = DataLoader(test_dataset, batch_size=funkcje.BATCH_SIZE, shuffle=True)

    model = GNN(in_channels=test_dataset[0].x.shape[1])
    optimizer = torch.optim.Adam(model.parameters(), lr=0.005)

    train_losses_MSE = []
    valid_losses_MSE = []
    valid_losses_R2 = []
    test_mse = 0

    for epoch in range(epochs):
        model.train()
        train_mse = []

        for batch in train_loader:
            optimizer.zero_grad()
            outputs = model(batch)

            loss = loss_function(outputs, batch.y)
            train_mse.append(loss.item())

            loss.backward()
            optimizer.step()

        model.eval()
        r2.reset()
        valid_mse = []

        with torch.no_grad():
            for batch in valid_loader:
                outputs = model(batch)
                valid_mse.append(loss_function(outputs, batch.y).item())
                r2.update(outputs, batch.y)

        train_losses_MSE.append(sum(train_mse) / len(train_mse))
        valid_losses_MSE.append(sum(valid_mse) / len(valid_mse))
        valid_losses_R2.append(r2.compute().item())

    with torch.no_grad():
        for batch in test_loader:
            outputs = model(batch)
            test_mse = (loss_function(outputs, batch.y).item())

    plt.plot(train_losses_MSE)
    plt.plot(valid_losses_MSE)
    plt.show()

    print(f"Najniższy MSE treningowy: {min(train_losses_MSE)}")
    print(f"Najniższy MSE validacyjny: {min(valid_losses_MSE)}")
    print(f"MSE testowy: {test_mse}")

    print(f"R2: {max(valid_losses_R2)}")
    print(valid_losses_R2)
