# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:14:50 2022

@author: Angelina
"""

# loading the required modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf
import joblib
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error, make_scorer, r2_score
from sklearn.metrics import accuracy_score, log_loss
from sklearn.metrics import confusion_matrix
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from tensorflow.keras.callbacks import EarlyStopping

tf.config.optimizer.set_jit(True)

random_state = 1

# Datensatz laden
# train_data1000 = pd.read_csv('C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/train_data/training_dataset_ver1.csv', header = None)
train_data1000 = pd.read_csv('C:/Users/Valentin/Documents/10. Semester/WR2/Projekt 3 - FNN CLF/WR-2-3-FFN-CLF/resources/train_data/training_dataset_ver1.csv', header = None)
train_data = train_data1000.loc[~np.any(train_data1000 == 1000,axis = 1)] 
X = train_data.iloc[:,:-1].values
y = train_data.iloc[:,-1:].values

# Splitte Daten in Test und Trainingsdatensatz auf
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Input skalieren
normalizer = MinMaxScaler()
normalizer.fit(X)
X_train = normalizer.transform(X_train)
X_test = normalizer.transform(X_test)

# metrics = [tf.keras.metrics.BinaryAccuracy(),
#            tf.keras.metrics.BinaryCrossentropy(),
#            tf.keras.metrics.AUC(),
#            tf.keras.metrics.Precision(),
#            tf.keras.metrics.Recall(),
#            tf.keras.metrics.TrueNegatives()
#            ]
metrics = None

# Definiere Neuronales Netz
def build_model(learning_rate,n_hidden_layers,layer_size,activation,dropout_rate):
    model = Sequential()
    
    for i in range(n_hidden_layers):
        model.add(tf.keras.layers.Dropout(dropout_rate, seed=random_state+i))
        model.add(Dense(units=layer_size,activation = activation))
         
    # Output Layer
    model.add(Dense(1, activation = 'sigmoid'))
        
    model.compile(loss = 'binary_crossentropy',
                  optimizer = tf.keras.optimizers.Adam(lr=learning_rate),
                  metrics = metrics)
    return model

model = KerasClassifier(build_model, verbose = 0)

# Schnelle Testaparameter
# parameters = {
#     'activation' : ['sigmoid'],  
#     'batch_size': [32],
#     'epochs': [10],
#     'learning_rate' : [0.01],
#     'n_hidden_layers' : [1],
#     'layer_size' : [2,6,10],
#     'dropout_rate' : [0]
#     }


# Weitere parameters
# 'data_percentage' : [1/3,1/2,1]
# init = ['glorot_uniform', 'normal', 'uniform']

#activation      = ['sigmoid','relu']
#batch_size      = [8,16,32,64]
#epochs          = [100,200,300]
#learning_rate   = [0.001,0.01,0.1]
#n_hidden_layers = [1,2,3,4]
#layer_size      = [100,200,500]
#dropout_rate    = [0, 0.1, 0.2, 0.5]

activation      = ['relu']
batch_size      = [8]
epochs          = [100]
learning_rate   = [0.001]
n_hidden_layers = [1]
layer_size      = [500]
dropout_rate    = [0.2]

parameters = {
    'activation' : activation,  
    'batch_size': batch_size,
    'epochs': epochs,
    'learning_rate' : learning_rate,
    'n_hidden_layers' : n_hidden_layers,
    'layer_size' : layer_size,
    'dropout_rate' : dropout_rate
    }


cv = StratifiedKFold(n_splits = 4, shuffle = True, random_state = random_state)

def false_negatives_loss(y_true,y_pred):
    # Berechnet den Anteil an FN/(FN+TP)
    fn = sum(np.all([np.array(y_true) == 1,np.array(y_pred) == 0],0)) 
    n_positives = sum(np.array(y_true) == 1)
    return fn[0]/n_positives[0]

# Definiere score-function
fn_loss = make_scorer(false_negatives_loss,greater_is_better = False)
float32neg_log_loss = make_scorer(log_loss, eps = 1e-7)
def confusion_matrix_scorer(clf, X, y):
     y_pred = clf.predict(X)
     cm = confusion_matrix(y, y_pred)
     return {'tn': cm[0, 0], 'fp': cm[0, 1],
             'fn': cm[1, 0], 'tp': cm[1, 1]}
randomsearch = RandomizedSearchCV(estimator = model,
                                  param_distributions = parameters,
                                  n_iter = 30,
                                  scoring={'accuracy':'accuracy',
                                           'ce':float32neg_log_loss,
                                           #'conf_matrix': confusion_matrix_scorer,
                                           'fn_rate': fn_loss},   
                                  refit = 'accuracy',
                                  cv=cv,
                                  random_state = random_state,
                                  return_train_score=True,
                                  verbose = 4)

# callback = EarlyStopping(monitor='val_loss', min_delta=0.0001, patience=10)
# callback = EarlyStopping(monitor='val_loss', min_delta=0.0001, patience=20)

# Wende Grid Search auf Daten an
randomsearch.fit(X_train, y_train, verbose = 0) 

result_df1 = pd.DataFrame.from_dict(randomsearch.cv_results_)
result_df2 = pd.DataFrame.from_dict(randomsearch2.cv_results_)
result_df3 = pd.DataFrame.from_dict(randomsearch3.cv_results_)
result_df = pd.concat([result_df1,result_df2,result_df3])
result_df = result_df.drop_duplicates(subset = ['param_' + parameter for parameter in parameters.keys()])
best_estim=randomsearch.best_estimator_

print('Train Accuracy: %0.2f' % randomsearch.score(X_train, y_train))
print('Test Accuracy: %0.2f' % randomsearch.score(X_test, y_test))
print('Best Parameters:',randomsearch.best_params_)
print('Best Estimator:',best_estim)
print('Best score:',randomsearch.best_score_)



# Vorhersage
y_test_pred=best_estim.predict(X_test)
# print("Scores des besten Classifier auf den Testdaten:")
# for m in metrics:
#     m.update_state(y_test,y_test_pred)
#     with tf.Session() as sess:
#         m.result().eval()
#     print(m.name + ": %.4f" % m.result().numpy())





 # Berechne Fehler
# mse = mean_squared_error(labeltr_pred,y_train)

# print('Best Parameters:',randomsearch.best_params_)
# print('Best Estimator:',best_estim)

# r2 = r2_score(labeltr_pred,y_train)
# print("Traindata: MSE: %.2f" % mse)
# print("Traindata: R2 : %.2f" % r2)
# print("Traindata: ACC: %.2f" % accuracy_score(labeltr_pred,y_train))

# # Testen
# labelpred=best_estim.predict(X_test)

# # Berechne Fehler
# mse = mean_squared_error(y_test, labelpred)
# r2 = r2_score(y_test, labelpred)
# print("Testdata: MSE: %.2f" % mse)
# print("Testdata: R2: %.2f" % r2)
# print("Testdata: ACC: %.2f" % accuracy_score(labelpred,y_test))

# x_ax = range(len(y_test))
# plt.scatter(x_ax, y_test, s=5, color="blue", label="original")
# plt.plot(x_ax, labelpred, lw=0.8, color="red", label="predicted")
# plt.legend()
# plt.show()

# Save trained model
joblib_file = "C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/trained_model/trained_model.pkl"
#joblib.dump(best_estim, joblib_file)
best_estim.save(joblib_file)


