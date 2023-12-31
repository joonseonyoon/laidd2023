import torch
import numpy as np
import pandas as pd

torch.manual_seed(1234)
torch.cuda.manual_seed(1234)
torch.cuda.manual_seed_all(1234)
np.random.seed(1234)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


import torch_geometric
from torch_geometric.data import Data # basic graph class
from torch_geometric.loader import DataLoader
from torch.optim import Adam


from torch.nn import Linear
import torch.nn.functional as F
from torch.nn.functional import gelu
from torch_geometric.nn import GCNConv, BatchNorm
from torch_geometric.nn import global_mean_pool
from torch.optim.lr_scheduler import StepLR, ReduceLROnPlateau

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from tqdm import tqdm
import os

torch.manual_seed(820823)

pocket2 = pd.read_excel('C:/Users/USER/Desktop/laidd/Pocket_2_binding_affinity.xlsx')

def check_atoms(mol):
  """
  This function checks whether all atoms are valid
  """
  valid_atoms = ('H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I')
  flag = True
  for atm in mol.GetAtoms(): # 분자 안에 존재하는 모든 atom에 대해서 loop을 수행.
    if atm.GetSymbol() not in valid_atoms:
      flag = False
      break
  return flag
  

ligand_data = []
for smi, score in tqdm(zip(pocket2["SMILES_CODE"], pocket2["BINDING_AFFINITY"])):
  mol = Chem.MolFromSmiles(smi)
  score = float(score)
  if check_atoms(mol): # 모든 atom들이 원하는 atom들로만 구성되어 있을 때.
    ligand_data.append((mol, score))
  else:
    print(f"{smi} has non-targeted atoms")
  continue
  

def convert_mol_to_graph(mol, use_pos = False):

    #mol2 = Chem.AddHs(mol)
    mol2 = Chem.RemoveHs(mol)

    n_bonds = len(mol2.GetBonds()) # 분자의 공유 결합 개수
    n_atoms = len(mol2.GetAtoms()) # 분자의 원자 개수

    node_attr = []
    #### node 속성 계산 시작 ####
    # RDKit으로 계산할 수 있는 Atom의 속성은 아래 링크에서 확인할 수 있다.
    # https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Atom
    #
    #.        0.   1.   2.   3.   4.   5.   6.   7.   8.    9.  10.
    #        'H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl' 'Br'  'I'
    valid_atoms = {'H': 0, 'B':1, 'C':2, 'N':3, 'O':4, 'F':5, 'P':6, 'S':7, 'Cl':8, 'Br':9, 'I':10}

    for atm_id in range(n_atoms):
        # Select an atom.
        atm = mol2.GetAtomWithIdx(atm_id)

        # Atom symbol check (9-dim)
        sym = atm.GetSymbol()
        atm_one_hot = [0] * len(valid_atoms) # 0이 9개 들어있는 리스트를 만든다.
        idx = valid_atoms[sym] # sym에 해당하는 원소 기호가 몇 번째에 있는지?
        atm_one_hot[idx] = 1   # 해당되는 원소의 위치만 1로 바꾼다.

        # Check hybridization (7-dim)
        hybrid = atm.GetHybridization()
        hybrid_one_hot = [0] * 7 # [0, 0, 0, 0, 0, 0, 0]
        if hybrid == Chem.HybridizationType.SP3:
          hybrid_one_hot[0] = 1
        elif hybrid == Chem.HybridizationType.SP2:
          hybrid_one_hot[1] = 1
        elif hybrid == Chem.HybridizationType.SP:
          hybrid_one_hot[2] = 1
        elif hybrid == Chem.HybridizationType.S:
          hybrid_one_hot[3] = 1
        elif hybrid == Chem.HybridizationType.SP3D:
          hybrid_one_hot[4] = 1
        elif hybrid == Chem.HybridizationType.SP3D2:
          hybrid_one_hot[5] = 1
        else: # hybridization이 제대로 정의되지 않은 나머지의 모든 경우.
          hybrid_one_hot[6] = 1

        # aromatic 인지 아닌지?  (True/False)
        if atm.GetIsAromatic():
          arom = 1
        else:
          arom = 0

        # ring 안에 존재하는지 아닌지? (True/False)
        if atm.IsInRing():
          ring_flag = 1
        else:
          ring_flag = 0

        # Degree (공유 결합의 개수)  (6-dim, one-hot)
        # 0, 1, 2, 3, 4, >=5
        degree_one_hot = [0, 0, 0, 0, 0, 0]
        degree = atm.GetTotalDegree()
        if degree >= 5: # 5개 이상의 공유 결합을 가지는 원자.
          degree_one_hot[5]=1
        else:
          degree_one_hot[degree]=1

        # Number of hydrogens (5-dim, one-hot)
        # 결합되어 있는 수소의 개수.
        # 0, 1, 2, 3, >=4
        num_h = atm.GetTotalNumHs()
        hydrogen_one_hot = [0, 0, 0, 0, 0]
        if num_h >= 4:
          hydrogen_one_hot[4] = 1
        else:
          hydrogen_one_hot[num_h] = 1

        # Chirality (4-dim, one-hot)
        chiral = atm.GetChiralTag()
        if chiral == Chem.rdchem.ChiralType.CHI_OTHER:
          chiral_one_hot = [1, 0, 0, 0]
        # Counter-clock-wise (반시계)
        elif chiral == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
          chiral_one_hot = [0, 1, 0, 0]
        # Clockwise (시계방향)
        elif chiral == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
          chiral_one_hot = [0, 0, 1, 0]
        # Chirality 정의되지 않음.
        elif chiral == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
          chiral_one_hot = [0, 0, 0, 1]

        # 원자 특성 계산 [원자 symbol one-hot, 공유 결합 개수, 전체 valence의 개수 (explicit + implicit), is an atom aromatic (True/False)? ]
        # 더 추가 가능!
        # in total 25-dim.
        # 11-dim, 7-dim, 6-dim, 5-dim, 4-dim, 4-dim = 35-dim
        attr = atm_one_hot + \
                hybrid_one_hot + \
                degree_one_hot + \
                hydrogen_one_hot + \
                chiral_one_hot + \
                [arom, ring_flag, atm.GetFormalCharge(), atm.GetNumRadicalElectrons()]

        #print(atm_id, attr)
        node_attr.append(attr)

    #### node 속성 계산 완료 ####

    edge_index = []
    edge_attr = []
    edge_weight = []
    for edge_idx in range(n_bonds): # 전체 공유 결합에 대해서 loop을 돌린다.

        bond = mol2.GetBondWithIdx(edge_idx) # 각 공유 결합에 대해서 시작 atom과 끝 atom의 인덱스를 확인.
        edge_index.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
        edge_index.append([bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()]) # undirected graph를 만들기 위해서 순서를 바꿔서 edge를 2번 넣어준다.

        # BondType (4-dimensional one-hot)
        btype = bond.GetBondType() # 공유 결합의 종류.
        if btype == Chem.rdchem.BondType.SINGLE:
            bond_one_hot = [1, 0, 0, 0]
            edge_weight.extend([1.0, 1.0])
        elif btype == Chem.rdchem.BondType.AROMATIC:
            bond_one_hot = [0, 1, 0, 0]
            edge_weight.extend([1.5, 1.5])
        elif btype == Chem.rdchem.BondType.DOUBLE:
            bond_one_hot = [0, 0, 1, 0]
            edge_weight.extend([2.0, 2.0])
        elif btype == Chem.rdchem.BondType.TRIPLE:
            bond_one_hot = [0, 0, 0, 1]
            edge_weight.extend([3.0, 3.0])

        # BondStereo (6-dimensional one-hot)
        stype = bond.GetStereo()
        if stype == Chem.rdchem.BondStereo.STEREOANY:
          stereo_one_hot = [1, 0, 0, 0, 0, 0]
        elif stype == Chem.rdchem.BondStereo.STEREOCIS:
          stereo_one_hot = [0, 1, 0, 0, 0, 0]
        elif stype == Chem.rdchem.BondStereo.STEREOE:
          stereo_one_hot = [0, 0, 1, 0, 0, 0]
        elif stype == Chem.rdchem.BondStereo.STEREONONE:
          stereo_one_hot = [0, 0, 0, 1, 0, 0]
        elif stype == Chem.rdchem.BondStereo.STEREOTRANS:
          stereo_one_hot = [0, 0, 0, 0, 1, 0]
        elif stype == Chem.rdchem.BondStereo.STEREOZ:
          stereo_one_hot = [0, 0, 0, 0, 0, 1]

        # Is this bond included in a ring?
        if bond.IsInRing():
          ring_bond = 1
        else:
          ring_bond = 0

        # Is this bond a conjugated bond?
        if bond.GetIsConjugated():
          conjugate = 1
        else:
          conjugate = 0

        # In total 12-dimensional edge attribute
        # bond-type (4-dim), bondstereo (6-dim), (ring, conjugate)
        # Can you image more?
        attr = bond_one_hot + stereo_one_hot + [ring_bond, conjugate] # 12 차원의 공유 결합 속성.

        # 분자는 undirected graph이므로 edge가 두 번 정의된다.
        # 그러므로 동일한 attribute를 두 번 넣어주어야 한다.
        edge_attr.append(attr)
        edge_attr.append(attr)
    #### edge 속성 계산 완료 ####


    # PyTorch Tensor로 변환.
    edge_attr = torch.tensor(edge_attr, dtype = torch.float)
    node_attr = torch.tensor(node_attr, dtype = torch.float)
    edge_index = torch.tensor(edge_index, dtype = torch.long)
    edge_index = edge_index.t().contiguous()
    edge_weight = torch.tensor(edge_weight, dtype = torch.float)


    # 만일 3D 좌표 정보를 사용한다면
    if use_pos:
        val = AllChem.EmbedMolecule(mol2)
        if val !=0:
          print(f"Error while generating 3D: {Chem.MolToSmiles(mol)}")
          return None

        pos_list = [] # this is optional
        for atm_id in range(n_atoms):
          # Get Atomic Position.
          atm_pos = mol2.GetConformer(0).GetAtomPosition(atm_id)
          crd = [atm_pos.x, atm_pos.y, atm_pos.z]
          pos_list.append(crd)

        pos = torch.tensor(pos_list, dtype=torch.float)
    else:
      pos = None
    #print(edge_attr.shape)
    #print(node_attr.shape)
    #print(edge_index.shape)

    return edge_index, node_attr, edge_attr, pos, edge_weight

#분자들을 하나씩 PyG의 Data type으로 변환
fname='pocket2_converted_graph.joblib'
data_list = []

for mol, score in tqdm(ligand_data[:100000]):
  result = convert_mol_to_graph(mol)
  if result is None:
    continue

  edge_index, node_attr, edge_attr, pos, edge_weight = result
  y = torch.tensor([[score]], dtype = torch.float)
  dtmp = Data(x = node_attr,
              edge_index = edge_index,
              edge_attr = edge_attr,
              edge_weight = edge_weight,
              pos = pos, y = y
              )
  data_list.append(dtmp)
joblib.dump(data_list, fname)
print(f"{fname} is saved!")

import random
random.seed(12345)
random.shuffle(data_list)
n_data = len(data_list)

train_set, val_set, test_set = data_list[:int(n_data*0.7)], data_list[int(n_data*0.7):int(n_data*0.8)], data_list[int(n_data*0.8):]

print(f"Number of training set: {len(train_set)}")
print(f"Number of Validaation set: {len(val_set)}")
print(f"Number of test set: {len(test_set)}")

train_loader = DataLoader(train_set, batch_size=32, shuffle=True, drop_last = False)
test_loader = DataLoader(test_set, batch_size=32, shuffle=True, drop_last = False)
val_loader = DataLoader(val_set, batch_size=32, shuffle=False, drop_last = False)



#Graph Attentive Network (GAT)을 정의


class my_GAT_model(torch.nn.Module):
    def __init__(self, n_feat, n_hidden):
        super().__init__()
        self.conv1 = GATConv(n_feat, n_hidden) # 가장 기본적인 graph convolution model, https://pytorch-geometric.readthedocs.io/en/latest/modules/nn.html#torch_geometric.nn.conv.GCNConv
        self.conv2 = GATConv(n_hidden, n_hidden)
        self.conv3 = GATConv(n_hidden, n_hidden)

        self.lin1 = Linear(n_hidden, int(n_hidden/2))
        self.lin2 = Linear(int(n_hidden/2), 1)

    def forward(self, data):
        x, edge_index, batch, edge_attr, edge_weight = data.x, data.edge_index, data.batch, data.edge_attr, data.edge_weight

        x = self.conv1(x, edge_index, edge_attr) # node_attribute와 edge_index만 사용하고 있음. -> edge의 attribute는 사용하고 있지 않음.
        x = F.gelu(x)
        x = F.dropout(x, training=self.training)

        x = self.conv2(x, edge_index, edge_attr)
        x = F.gelu(x)

        x = self.conv3(x, edge_index, edge_attr)
        x = F.gelu(x)

        x = global_mean_pool(x, batch) # 전체의 node feature의 평균 값을 취한다. # [batch_size, hidden_channels]
        x = F.gelu(self.lin1(x))
        x = self.lin2(x)

        return x
		
model = my_GAT_model(37, 74)
model.to(device)

def test(loader):
    model.eval()
    error = 0.0
    out_all = []
    true = []
    for data in loader:  # Iterate in batches over the training/test dataset.
        out = model(data.to(device))
        tmp = (out - data.y)**2
        error += tmp.sum().item()  # Check against ground-truth labels.

        out_all.extend([x.item() for x in out])
        true.extend([x.item() for x in data.y])

    return error / len(loader.dataset), out_all, true  # Derive ratio of correct predictions.
	

def train():
  for idx, batch in enumerate(train_loader):
    out = model(batch.to(device))
    loss = criterion(out, batch.y)
    loss.backward()
    optimizer.step()  # Update parameters based on gradients.
    optimizer.zero_grad()  # Clear gradients.
    if idx%100 == 0:
      lr = scheduler.get_last_lr()
      print(f"IDX: {idx:5d}\tLoss: {loss:.4f}\tLR:{lr}")


optimizer = Adam(model.parameters(), lr=0.005)
criterion = torch.nn.MSELoss()


train_acc_list = []
val_acc_list = []

best_val_rmse = 9999.9

# Learning rate scheduler
scheduler = StepLR(optimizer, step_size=25, gamma=0.5)
model.train()

n_epoch = 10

for epoch in range(1, n_epoch):
    print("="*100)
    print("Epoch: ", epoch)

    train()
    scheduler.step()

    train_acc, out_tmp, true_tmp = test(train_loader)
    train_acc_list.append(train_acc)

    val_acc, val_pred, val_true = test(val_loader)
    val_acc_list.append(val_acc)

    # Save the best model
    if val_acc < best_val_rmse:
      best_epoch = epoch
      best_val_rmse = val_acc
      best_train_rmse = train_acc
      torch.save(model, "Best_GCN_model_v1.pt")

    print("-"*80)
    print(f'Epoch: {epoch:03d}, Train RMSE: {train_acc:.4f}, Val RMSE: {val_acc:.4f}')
    print("-"*80)
    print(f"Best Epoch: {best_epoch:03d}")
    print(f"Best val   RMSE: {best_val_rmse:.4f}")
    print(f"Best train RMSE: {best_train_rmse:.4f}")
    print('='*80)


test_acc, out_all, true_all = test(test_loader)
print(f"Test set RMSE: {test_acc:.4f}")

# plotting to investigate the correlation of docking and predicted scores
import matplotlib.pyplot as plt

experimental = [x for x in true_all]
prediction = [x for x in out_all]

plt.figure(figsize=(8, 8))
plt.scatter(experimental, prediction, marker = '.')
plt.plot(range(-12, 0), range(-12, 0), 'r--')
plt.xlabel("docking score", fontsize='xx-large')
plt.ylabel("predicted score", fontsize='xx-large')
plt.xlim(-12, 0)
plt.ylim(-12, 0)
