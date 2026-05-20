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
            nn.Dropout(0.1),
        )

    def forward(self, x):
        return self.model(x)


df_activity = pd.read_csv('activity.csv', index_col=False)
df_activity = df_activity[['canonical_smiles', 'standard_value']]

train, test, valid = funkcje.scaffold_split_df(df_activity)
# train, test, valid = funkcje.random_split(df_activity)

datasets = [train, test, valid]

for dataset in datasets:
    dataset['canonical_smiles'] = dataset['canonical_smiles'].apply(Chem.MolFromSmiles)
    dataset['canonical_smiles'].dropna(inplace=True)  # rdkit wpisuje 'None' w niemożliwych cząsteczkach
    dataset['canonical_smiles'] = dataset['canonical_smiles'].apply(funkcje.morgan_fp)
    dataset['canonical_smiles'] = dataset['canonical_smiles'].apply(np.array)

    print(dataset.info())

X_train, X_test, X_validate = train.drop('standard_value', axis=1), test.drop('standard_value', axis=1), valid.drop('standard_value', axis=1)
y_train, y_test, y_validate = train['standard_value'], test['standard_value'], valid['standard_value']

train_dataset = TensorDataset(torch.FloatTensor(np.array(X_train.values.tolist())), torch.FloatTensor(y_train.values))
train_loader = DataLoader(train_dataset, batch_size=funkcje.BATCH_SIZE, shuffle=True)

test_dataset = TensorDataset(torch.FloatTensor(np.array(X_test.values.tolist())), torch.FloatTensor(y_test.values))
test_loader = DataLoader(test_dataset, batch_size=funkcje.BATCH_SIZE, shuffle=True)

model = MLP()
optimizer = torch.optim.Adam(model.parameters(), lr=0.005)
#loss_function = nn.MSELoss()
loss_function = R2Score()

train_losses = []
test_losses = []
mean_valid_losses = []
epochs = 10  # 100  # zwykle im więcej parametrow tym wiecej epok
# minimum average distace, mówi ile max epok
# virtual node może pomóc żeby mocno oddalone neurony się o sobie dowiedziały
# learnning rate inny na poczatek i potem

for epoch in range(epochs):
    model.train()
    train_losses.append([])
    test_losses.append([])

    for batch_id, (inputs, targets) in enumerate(train_loader):
        optimizer.zero_grad()
        outputs = model(inputs).squeeze()

        loss = loss_function(outputs, targets)
        train_losses[epoch].append(loss.item())

        loss.backward()
        optimizer.step()

    for batch_id, (inputs, targets) in enumerate(test_loader):
        outputs = model(inputs).squeeze()
        loss = loss_function(outputs, targets)
        test_losses[epoch].append(loss.item())

train_average = [sum(e) / len(e) for e in train_losses]
test_average = [sum(e) / len(e) for e in test_losses]

plt.plot(train_average)
plt.plot(test_average)
plt.show()

print(f"Best train score: {min(train_average)}")
print(f"Best test score: {min(test_average)}")
print('?')
# jedak PIC50
# ogranizenie aktywosci po bialku
# sprawdzanie wartosci neuronow, sprawdzanie czy nie umieraja
# sprawdzanie gradientow czy nie wybuchaja ani nie wygasaja
# early stop w zalezninosci od loss fuction
# overfitting point? raczej nie wsytąpi
# R2 ~0.68, RMSE ~0.8
# hugging face: gemma?
