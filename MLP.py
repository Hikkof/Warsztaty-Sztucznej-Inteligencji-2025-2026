import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from torchmetrics.regression import R2Score
from rdkit import Chem

import funkcje


class MLP(nn.Module):
    def __init__(self):
        super(MLP, self).__init__()

        self.model = nn.Sequential(
            nn.Linear(funkcje.FP_SIZE, 1024),
            nn.ReLU(),
            nn.Linear(1024, 512),
            nn.ReLU(),
            nn.Linear(512, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
            # nn.Dropout(0.1),
        )

    def forward(self, x):
        return self.model(x)


df_activity = pd.read_csv('activity.csv', index_col=False)
df_activity = df_activity[['canonical_smiles', 'standard_value']]

df_activity['canonical_smiles'] = df_activity['canonical_smiles'].apply(Chem.MolFromSmiles)
df_activity['canonical_smiles'].dropna(inplace=True)  # rdkit wpisuje 'None' w niemożliwych cząsteczkach

epochs = 10
loss_functions = {'MSELoss': nn.MSELoss(), 'R2Score': R2Score()}
split_types = {'Random': funkcje.random_split(df_activity), 'Scaffold': funkcje.scaffold_split_df(df_activity)}

for split in split_types:
    print(f"Rodzaj podziału: {split}")
    train, test, valid = split_types[split]

    datasets = [train, test, valid]
    for dataset in datasets:
        dataset['canonical_smiles'] = dataset['canonical_smiles'].apply(funkcje.morgan_fp)
        dataset['canonical_smiles'] = dataset['canonical_smiles'].apply(np.array)

    X_train, X_test, X_validate = train.drop('standard_value', axis=1), test.drop('standard_value', axis=1), valid.drop(
        'standard_value', axis=1)
    y_train, y_test, y_validate = train['standard_value'], test['standard_value'], valid['standard_value']

    train_dataset = TensorDataset(torch.FloatTensor(np.array(X_train.values.tolist())), torch.FloatTensor(y_train.values))
    train_loader = DataLoader(train_dataset, batch_size=funkcje.BATCH_SIZE, shuffle=True)

    test_dataset = TensorDataset(torch.FloatTensor(np.array(X_test.values.tolist())), torch.FloatTensor(y_test.values))
    test_loader = DataLoader(test_dataset, batch_size=funkcje.BATCH_SIZE, shuffle=True)

    for loss_function in loss_functions:
        print(f"Funkcja straty: {loss_function}")

        model = MLP()
        optimizer = torch.optim.Adam(model.parameters(), lr=0.005)

        train_losses = []
        test_losses = []

        for epoch in range(epochs):
            model.train()
            train_losses.append([])

            for batch_id, (inputs, targets) in enumerate(train_loader):
                optimizer.zero_grad()
                outputs = model(inputs).squeeze()

                loss = loss_functions[loss_function](outputs, targets)
                train_losses[epoch].append(loss.item())

                loss.backward()
                optimizer.step()

            model.eval()
            test_losses.append([])

            with torch.no_grad():
                for batch_id, (inputs, targets) in enumerate(train_loader):
                    outputs = model(inputs).squeeze()
                    loss = loss_functions[loss_function](outputs, targets)
                    test_losses[epoch].append(loss.item())

        train_average = [sum(e) / len(e) for e in train_losses]
        test_average = [sum(e) / len(e) for e in test_losses]
        '''
        plt.plot(train_average)
        plt.plot(test_average)
        plt.show()
        '''
        print(f"Best train score: {min(train_average)}")
        print(f"Best test score: {min(test_average)}")
        
print('?')
