from abc import ABC

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
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

        self.conv1 = GCNConv(in_channels, 64)
        self.conv2 = GCNConv(64, 64)

        self.lin1 = nn.Linear(64, 32)
        self.lin2 = nn.Linear(32, 1)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index).relu()
        x = self.conv2(x, edge_index).relu()

        x = global_mean_pool(x, batch)
        x = self.lin1(x).relu()
        x = self.lin2(x)

        return x.squeeze(1)


def mol_to_graph(mol, target):
    atom_features = []
    for atom in mol.GetAtoms():
        atom_features.append([
            atom.GetAtomicNum(),
            atom.GetTotalDegree(),
            atom.GetFormalCharge(),
            atom.GetTotalNumHs(),
            atom.GetHybridization().real,
            atom.GetIsAromatic()
        ])
    x = torch.tensor(atom_features, dtype=torch.float)

    edge_index = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        edge_index.append([i, j])
        edge_index.append([j, i])

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()

    y = torch.tensor([target], dtype=torch.float)

    return Data(x=x, edge_index=edge_index, y=y)


class CustomDataset(InMemoryDataset, ABC):
    def __init__(self, data_frame):
        super().__init__()
        data_list = []

        for smiles, target in zip(data_frame["canonical_smiles"], data_frame["standard_value"]):
            graph = mol_to_graph(smiles, target)
            data_list.append(graph)

        self.data, self.slices = self.collate(data_list)


df_activity = pd.read_csv('initial cleaning/activity.csv.csv', index_col=False)
df_activity = df_activity[['canonical_smiles', 'standard_value']]

df_activity['canonical_smiles'] = df_activity['canonical_smiles'].apply(Chem.MolFromSmiles)
df_activity['canonical_smiles'].dropna(inplace=True)

epochs = 250
loss_functions = {'MSELoss': nn.MSELoss(), 'R2Score': R2Score()}
split_types = {'Random': funkcje.random_split(df_activity), 'Scaffold': funkcje.scaffold_split_df(df_activity)}

for split in split_types:
    print(f"Rodzaj podziału: {split}")
    train, test, valid = split_types[split]

    train_dataset = CustomDataset(train)
    train_loader = DataLoader(train_dataset, batch_size=funkcje.BATCH_SIZE, shuffle=True)

    test_dataset = CustomDataset(test)
    test_loader = DataLoader(test_dataset, batch_size=funkcje.BATCH_SIZE, shuffle=True)

    for loss_function in loss_functions:
        print(f"Funkcja straty: {loss_function}")

        model = GNN(in_channels=test_dataset[0].x.shape[1])
        optimizer = torch.optim.Adam(model.parameters(), lr=0.005)

        train_losses = []
        test_losses = []

        for epoch in range(epochs):
            model.train()
            train_losses.append([])

            for batch in train_loader:
                optimizer.zero_grad()
                outputs = model(batch)

                loss = loss_functions[loss_function](outputs, batch.y)
                train_losses[epoch].append(loss.item())

                loss.backward()
                optimizer.step()

            model.eval()
            test_losses.append([])

            with torch.no_grad():
                for batch in test_loader:
                    outputs = model(batch)
                    loss = loss_functions[loss_function](outputs, batch.y)
                    test_losses[epoch].append(loss.item())

        train_average = [sum(e) / len(e) for e in train_losses]
        test_average = [sum(e) / len(e) for e in test_losses]
        '''
        plt.plot(train_average)
        plt.plot(test_average)
        plt.show()
        '''

        print(f"Najniższy wynik treningowy: {min(train_average)}")
        print(f"Najniższy wynik testowy: {min(test_average)}")

print('?')

# training loop, early stop zmiana loss mniej niż 1% względem ostatniej epoki\
# porównywanie loss do ostatniej epoki i zmniejszanie lr\
# jeśli average z kilku poprzednich epok wystarczająco blisko 0 też należy zmniejszyć lr

# sprawdzanie neuronów żeby nie wygasały\
# sprawdzanie gradientów

# zapisanie modelu
