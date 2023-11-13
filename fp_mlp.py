import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.data as data
import pandas as pd
import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import numpy as np

device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")

class Classifier(nn.Module):
    def __init__(self, num_inputs, num_hidden, num_outputs):
        super().__init__()
        # 네트워크를 구성하는 위한 기반이 되는 layer들을 정의한다.
        self.linear1 = nn.Linear(num_inputs, num_hidden)
        self.linear2 = nn.Linear(num_hidden, num_hidden)
        self.linear3 = nn.Linear(num_hidden, num_hidden)
        self.linear4 = nn.Linear(num_hidden, num_outputs)
        self.act_fn = nn.Tanh()

    def forward(self, x):
        # Perform the calculation of the model to determine the prediction
        x = self.linear1(x)
        x = self.act_fn(x)

        x = self.linear2(x)
        x = self.act_fn(x)

        x = self.linear3(x)
        x = self.act_fn(x)

        x = self.linear4(x)
        # x = self.act_fn(x)

        return x
		

class DockingDataset(data.Dataset):
    def __init__(self):

        super().__init__()
        self.raw_df=pd.read_excel('/content/drive/MyDrive/Company/laidd/pocket_2_binding_affinity.xlsx')
        self.df_sub=self.raw_df[["SMILES_CODE", "BINDING_AFFINITY"]].copy() #SMILES_CODE, BINDING_AFFINITY 열만 추출
        self.df_sub=self.df_sub.dropna(axis='index', subset='BINDING_AFFINITY') # BINDING_AFFINITY column에 NaN이 있는 raw를 제거

        fp_list = []
        for smi in tqdm(self.df_sub["SMILES_CODE"]):
          m = Chem.MolFromSmiles(smi)
          fp = AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=1024)
          fp_list.append(fp.ToList())

        self.data=torch.tensor(fp_list, dtype=torch.float32)  #input data to tensor
        print("self.data:")
        print(self.data)
        print("self.data.shape:")
        print(self.data.shape)

        self.label=torch.tensor(self.df_sub["BINDING_AFFINITY"].values, dtype=torch.float32)  # target_label to tensor
        print("self.label:")
        print(self.label)
        print(self.label.shape)

    def __len__(self):
        # 전체 데이터의 개수를 return하는 함수
        return len(self.data)

    def __getitem__(self, idx):
        # idx 번째 데이터와 레이블을 리턴하는 함수
        data_point = self.data[idx]
        data_label = self.label[idx]
        return data_point, data_label

pocket2 = DockingDataset()

dataset_size = len(pocket2)
train_size = int(dataset_size * 0.8)
validation_size = int(dataset_size * 0.1)
test_size = dataset_size - train_size - validation_size

train_dataset, validation_dataset, test_dataset = data.random_split(pocket2, [train_size, validation_size, test_size])

print(f"Training Data Size : {len(train_dataset)}")
print(f"Validation Data Size : {len(validation_dataset)}")
print(f"Testing Data Size : {len(test_dataset)}")

loss_module = nn.MSELoss()
model = Classifier(num_inputs=1024, num_hidden=512, num_outputs=1)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

train_loader = data.DataLoader(train_dataset, batch_size=128, shuffle=True)
validation_loader = data.DataLoader(validation_dataset, batch_size=128, shuffle=True)

# Push model to device. Has to be only done once
model.to(device)

def train_model(model, optimizer, data_loader, loss_module, num_epochs=100):
    # Set model to train mode
    model.train()
    train_loss_list = []
    validation_loss_list = []

    # Training loop
    for epoch in tqdm(range(num_epochs)):
        train_loss = 0.0
        for data_inputs, data_labels in train_loader:

            ## Step 1: Move input data to device (only strictly necessary if we use GPU)
            data_inputs = data_inputs.to(device)
            data_labels = data_labels.to(device)

            ## Step 2: Run the model on the input data
            preds = model(data_inputs)
            preds = preds.squeeze(dim=1) # Output is [Batch size, 1], but we want [Batch size]

            ## Step 3: Calculate the loss
            loss = loss_module(preds, data_labels.float())
            train_loss += loss.item() # sum over batches

            ## Step 4: Perform backpropagation
            # Before calculating the gradients, we need to ensure that they are all zero.
            # The gradients would not be overwritten, but actually added to the existing ones.
            optimizer.zero_grad()
            # Perform backpropagation
            loss.backward()

            ## Step 5: Update the parameters
            optimizer.step()

        # save a training loss at a given epoch
        train_loss_list.append(train_loss)
        ### end of training

        ### evaluation using the validation set
        model.eval()
        val_loss = 0.0
        for data_inputs, data_labels in validation_loader:

            ## Step 1: Move input data to device (only strictly necessary if we use GPU)
            data_inputs = data_inputs.to(device)
            data_labels = data_labels.to(device)

            ## Step 2: Run the model on the input data
            preds = model(data_inputs)
            preds = preds.squeeze(dim=1) # Output is [Batch size, 1], but we want [Batch size]

            ## Step 3: Calculate the loss
            loss = loss_module(preds, data_labels.float())
            val_loss += loss.item() # sum over batches

        validation_loss_list.append(val_loss)
        #Print information out
        if epoch % 5 == 0:
          print(f'Epoch: {epoch:03d}, Training Loss: {train_loss:.4f}, Validation Loss: {val_loss:.4f}')

        model.train()
    return train_loss_list, validation_loss_list
	
	
train_loss, validation_loss = train_model(model, optimizer, train_loader, loss_module)


### Evaluation
test_data_loader = data.DataLoader(test_dataset, batch_size=128, shuffle=False, drop_last=False)

pred_values = []
docking_values = []

for data_inputs, data_labels in test_loader:
    data_inputs = data_inputs.to(device)
    data_labels = data_labels.to(device)

    docking_values.append(data_labels.squeeze(dim=-1))
    pred_values.append(model(data_inputs).squeeze(dim=1))
    #preds1 = preds.squeeze(dim=1)

docking_values = [value.tolist() for value in docking_values]
pred_values = [value.tolist()[0] for value in pred_values]
plt.scatter(docking_values, pred_values)


np.corrcoef(docking_values, pred_values)
