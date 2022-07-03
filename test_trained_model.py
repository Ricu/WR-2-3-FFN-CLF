# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:42:09 2022

@author: Angelina
"""
#import numpy as np
#import pandas as pd
#from sklearn.metrics import mean_squared_error, make_scorer, r2_score
#from sklearn.metrics import accuracy_score
#from sklearn.externals import joblib

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
from sklearn.metrics import accuracy_score, log_loss, precision_score, recall_score
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
from tensorflow.keras.callbacks import EarlyStopping
from matplotlib import rcParams
from sklearn.metrics import confusion_matrix

# Datensatz laden
train_data1000 = pd.read_csv('C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/train_data/2022-06-29-05-12-03-train_data_dump.csv', header = None)
#train_data1000 = pd.read_csv('C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/train_data/training_dataset_ver1.csv', header = None)
# train_data1000 = pd.read_csv('C:/Users/Valentin/Documents/10. Semester/WR2/Projekt 3 - FNN CLF/WR-2-3-FFN-CLF/resources/train_data/training_dataset_ver1.csv', header = None)
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

metrics = [tf.keras.metrics.BinaryAccuracy(name='accuracy'),
           tf.keras.metrics.BinaryCrossentropy(),
           tf.keras.metrics.AUC(),
           tf.keras.metrics.Precision(name='precision'),
           tf.keras.metrics.Recall(name='recall'),
           tf.keras.metrics.TrueNegatives()
           ]
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

#model = KerasClassifier(build_model,learning_rate= 0.001,n_hidden_layers=1,layer_size=500,activation='relu',dropout_rate=0.2, verbose = 0)


activation      = 'relu'
batch_size      = 8
epochs          = 100
learning_rate   = 0.001
n_hidden_layers = 1
layer_size      = 500
dropout_rate    = 0.2

model = build_model(learning_rate,n_hidden_layers,layer_size,activation,dropout_rate)

history = model.fit(X_train, y_train, epochs = epochs, batch_size = batch_size, verbose = 0)  


# Plotts
rcParams['figure.figsize'] = (18, 8)
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
plt.plot(
    np.arange(1, 101), 
    history.history['loss'], label='Loss'
)
plt.plot(
    np.arange(1, 101), 
    history.history['accuracy'], label='Accuracy'
)
plt.plot(
    np.arange(1, 101), 
    history.history['precision'], label='Precision'
)
plt.plot(
    np.arange(1, 101), 
    history.history['recall'], label='Recall'
)
plt.title('Evaluation metrics', size=20)
plt.xlabel('Epoch', size=14)
plt.legend();

# Predict and classify train data
predictions = model.predict(X_train)
prediction_classes = [
    1 if prob > 0.45 else 0 for prob in np.ravel(predictions)
]

# Model evaluation
print("Train MSE: %.2f", mean_squared_error(y_train, prediction_classes))
print('Train Accuracy: %0.2f', accuracy_score(y_train, prediction_classes))
print('Train Confusion Matrix:',confusion_matrix(y_train, prediction_classes))
print('Train Precision: %0.2f', precision_score(y_train, prediction_classes))
print('Train Recall: %0.2f', recall_score(y_train, prediction_classes))

# Predict and classify test data
predictions = model.predict(X_test)
prediction_classes = [
    1 if prob > 0.45 else 0 for prob in np.ravel(predictions)
]

# Model evaluation
print("Test MSE: %.2f", mean_squared_error(y_test, prediction_classes))
print('Test Accuracy: %0.2f', accuracy_score(y_test, prediction_classes))
print('Test Confusion Matrix:',confusion_matrix(y_test, prediction_classes))
print('Test Precision: %0.2f', precision_score(y_test, prediction_classes))
print('Test Recall: %0.2f', recall_score(y_test, prediction_classes))


# Save trained model
file_name = "C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/trained_model/trained_model.pkl"
model.save(file_name)

# Load trained model
file_name = "C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/trained_model/trained_model.pkl"
trained_model = keras.models.load_model(file_name)


# Load test data
test_data = pd.read_csv("C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/test_data/2022-06-30-16-13-08-test_data_1_dump.csv", header = None)
X_test = test_data.iloc[:,:-1].values
y_test = test_data.iloc[:,-1:].values

# Predict and classify test data
predictions = model.predict(X_test)
prediction_classes = [
    1 if prob > 0.45 else 0 for prob in np.ravel(predictions)
]

# Model evaluation
print("Test MSE: %.2f", mean_squared_error(y_test, prediction_classes))
print('Test Accuracy: %0.2f', accuracy_score(y_test, prediction_classes))
print('Test Confusion Matrix:',confusion_matrix(y_test, prediction_classes))
print('Test Precision: %0.2f', precision_score(y_test, prediction_classes))
print('Test Recall: %0.2f', recall_score(y_test, prediction_classes))


# Save predicted labels
file_name = "C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/trained_model/predicted_labels_1.csv"
np.savetxt(file_name, prediction_classes, delimiter=",")

