# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:14:50 2022

@author: Angelina
"""

# loading the required modules
import numpy as np
import pandas as pd
import keras
import matplotlib.pyplot as plt
import timeit
from keras.models import Sequential
from keras.layers import Dense
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import AdaBoostRegressor
from sklearn.metrics import mean_squared_error, make_scorer, r2_score
from sklearn.metrics import accuracy_score
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from keras.callbacks import EarlyStopping

random_state = 0

# Datensatz laden
# train_data = pd.read_csv('C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/train_data/2022-06-28-11-05-08-train_data_dump.csv', header = None)
train_data = pd.read_csv('C:/Users/Valentin/Documents/10. Semester/WR2/Projekt 3 - FNN CLF/WR-2-3-FFN-CLF/resources/train_data/training_dataset_ver1.csv', header = None)
coeff = train_data.iloc[:,:-1].values
label = train_data.iloc[:,-1:].values

#xtrain, xtest, ytrain, ytest=train_test_split(x, y, test_size=0.15)

# Splitte Daten in Test und Trainingsdatensatz auf
coeff_train, coeff_test, label_train, label_test = train_test_split(coeff,
                                                                    label,
                                                                    test_size=0.2)


# Input skalieren
sc_st = MinMaxScaler()
sc_st.fit(coeff)
coeff_train = sc_st.transform(coeff_train)
coeff_test = sc_st.transform(coeff_test)

# Definiere Neuronales Netz
def build_model(learning_rate,n_hidden_layers,layer_size,activation):
    inputs = keras.Input(shape = (4800,), name = 'input_layer')
    model = keras.Sequential(inputs)
    # model = Sequential()
    #for i in range(hp.Int("n_layers", 1, 3)):
    for i in range(n_hidden_layers):
       # model.add(Dense(units=hp.Int(f"units_{i}", min_value=32, max_value=128, step=32),
       #                        activation = 'sigmoid'
        #    ))
         model.add(Dense(units=layer_size,activation = activation))
         
    model.add(Dense(1, activation = activation))
        
    model.compile(loss = keras.losses.MeanSquaredError(),
                  optimizer = keras.optimizers.Adam(lr=learning_rate),
                  metrics = ["accuracy"])
    return model

model = KerasClassifier(build_model, verbose = 0)

# Zu testende Parameter des Modells
#parameters = {
# 'batch_size': [16, 32],
# 'epochs': [100, 200],
# 'learning_rate' : [0.001, 0.005],
# 'layer_size' : [2],
# 'output_shape' : [100]
# }
# Weitere parameters
# 'activation' : ['sigmoid','relu']
# 'data_percentage' : [1/3,1/2,1]
# 'dropout'
# 'patience'
# 'min_delta'

activation      = ['sigmoid','relu']
batch_size      = [8,16,32,64]
epochs          = [100,200]
learning_rate   = [0.001,0.01,0.1]
n_hidden_layers = [1,2,3,4]
layer_size      = [100,500,1000]



parameters = {
    'activation' : activation,
        
    'batch_size': batch_size,
    'epochs': epochs,
    'learning_rate' : learning_rate,
    'n_hidden_layers' : n_hidden_layers,
    'layer_size' : layer_size
    }

# Definiere score-function
score = make_scorer(mean_squared_error)
cv = StratifiedKFold(n_splits = 4, shuffle = True, random_state = random_state)

gridsearch = GridSearchCV(estimator = model,
                          param_grid = parameters,
                          scoring=score,
                          cv=cv,
                          #return_train_score=True,
                          verbose = 4)

callback = EarlyStopping(monitor='val_loss', min_delta=0.0001, patience=10)

# Wende Grid Search auf Daten an
start = timeit.timeit()
gridsearch.fit(coeff_train, np.ravel(label_train), callbacks=[callback], validation_split=0.1, verbose = 0) 
end = timeit.timeit()


best_estim=gridsearch.best_estimator_


# Wende bestes Modell auf Daten an
# best_estim ist bereits auf coeff_train gefittet
# best_estim.fit(coeff_train, np.ravel(label_train), callbacks=[callback], validation_split=0.1)

# Vorhersage
labeltr_pred=best_estim.predict(coeff_train)

# Berechne Fehler
mse = mean_squared_error(labeltr_pred,label_train)

print('Benoetigte Zeit:',end - start)
print('Best Parameters:',gridsearch.best_params_)
print('Best Estimator:',best_estim)

r2 = r2_score(labeltr_pred,label_train)
print("Traindata: MSE: %.2f" % mse)
print("Traindata: R2: %.2f" % r2)
print("Traindata: ACC: %.2f" % accuracy_score(labeltr_pred,label_train))

# Testen
labelpred=best_estim.predict(coeff_test)

# Berechne Fehler
mse = mean_squared_error(label_test, labelpred)
r2 = r2_score(label_test, labelpred)
print("Testdata: MSE: %.2f" % mse)
print("Testdata: R2: %.2f" % r2)
print("Testdata: ACC: %.2f" % accuracy_score(labelpred,label_test))

x_ax = range(len(label_test))
plt.scatter(x_ax, label_test, s=5, color="blue", label="original")
plt.plot(x_ax, labelpred, lw=0.8, color="red", label="predicted")
plt.legend()
plt.show()

