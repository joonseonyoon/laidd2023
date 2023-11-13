#Mordred descriptor 이용한 분자 표현 및 random forest regressor 를 이용한 학습

import rdkit.Chem as Chem
import pandas as pd
from tqdm import tqdm
import rdkit.Chem.AllChem as AllChem
from mordred import Calculator, descriptors

pocket2_db = pd.read_excel('C:/Users/user/Desktop/laidd/docking/Pocket_2_binding_affinity.xlsx')
mols = [Chem.MolFromSmiles(x) for x in pocket2_db["SMILES_CODE"]]

calc = Calculator(descriptors, ignore_3D = True)

mordred = calc.pandas(mols)
mordred = mordred.astype(float).dropna(axis=1)

X = mordred
y = pocket2_db["BINDING_AFFINITY"]

import sklearn.model_selection
X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.2, random_state=42)

import sklearn.ensemble
from sklearn.ensemble import RandomForestRegressor

my_model = RandomForestRegressor(random_state=42, n_jobs=30) #병렬 처리

my_model.fit(X_train, y_train)

y_pred2 = my_model2.predict(X_test)

plt.scatter(y_test, y_pred2)
plt.xlabel("Calculated binding affinity")
plt.ylabel("Predicted binding affinity")
plt.grid()
plt.plot(range(-10, 0), range(-10, 0), "r--", label = "y=x")
plt.legend()


from sklearn.metrics import mean_squared_error
print(f"Mean Squared Error of this model is {mean_squared_error(y_test, y_pred2):.3f}")

import numpy as np
print(f"Pearson's correlation coefficient of our model is {np.corrcoef(y_test,y_pred2)[0,1]:.3f}")

with open('rfregressor_pocket2.pkl', 'wb') as file:
    pickle.dump(my_model, file)
	
# Load the saved model
with open("rfregressor_pocket2.pkl", 'rb') as file:
    model = pickle.load(file)

