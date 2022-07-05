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
train_data1000 = pd.read_csv('C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/train_data/training_dataset_ver1.csv', header = None)
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
predictions_train = model.predict(X_train)
prediction_classes_train_list = [
    1 if prob > 0.45 else 0 for prob in np.ravel(predictions_train)
]
prediction_classes_train = np.array(prediction_classes_train_list)

# Model evaluation
print("Train MSE: %.2f" % mean_squared_error(y_train, prediction_classes_train))
print('Train Accuracy: %.2f' % accuracy_score(y_train, prediction_classes_train))
print('Train Confusion Matrix:', confusion_matrix(y_train, prediction_classes_train))
print('Train Precision: %.2f' % precision_score(y_train, prediction_classes_train))
print('Train Recall: %.2f' % recall_score(y_train, prediction_classes_train))

# Predict and classify test data
predictions_test1 = model.predict(X_test)
prediction_classes_test1_list = [
    1 if prob > 0.45 else 0 for prob in np.ravel(predictions_test1)
]
prediction_classes_test1 = np.array(prediction_classes_test1_list)

# Model evaluation
print("Test MSE: %.2f" % mean_squared_error(y_test, prediction_classes_test1))
print('Test Accuracy: %.2f' % accuracy_score(y_test, prediction_classes_test1))
print('Test Confusion Matrix:', confusion_matrix(y_test, prediction_classes_test1))
print('Test Precision: %.2f' % precision_score(y_test, prediction_classes_test1))
print('Test Recall: %.2f' % recall_score(y_test, prediction_classes_test1))


# Save trained model
file_name = "C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/trained_model/trained_model.pkl"
model.save(file_name)

# Load trained model
file_name = "C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/trained_model/trained_model.pkl"
trained_model = tf.keras.models.load_model(file_name)


# Load test data
test_data = pd.read_csv("C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/test_data/2022-07-05-00-04-29-test_data_1_dump.csv", header = None)
X_test2 = test_data.iloc[:,:-1].values
y_test2 = test_data.iloc[:,-1:].values

# Predict and classify test data
predictions_test2 = trained_model.predict(X_test2)
prediction_classes_test2_list = [
    1 if prob > 0.45 else 0 for prob in np.ravel(predictions_test2)
]
prediction_classes_test2 = np.array(prediction_classes_test2_list)

# Model evaluation
print("Test MSE: %.2f" % mean_squared_error(y_test2, prediction_classes_test2))
print('Test Accuracy: %.2f' %  accuracy_score(y_test2, prediction_classes_test2))
print('Test Confusion Matrix:', confusion_matrix(y_test2, prediction_classes_test2))
print('Test Precision: %.2f' % precision_score(y_test2, prediction_classes_test2))
print('Test Recall: %.2f' % recall_score(y_test2, prediction_classes_test2))


# Save predicted labels
file_name = "C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/trained_model/predicted_labels_1.csv"
np.savetxt(file_name, prediction_classes_test2, delimiter=",")

# Load test data
test_data3 = pd.read_csv("C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/test_data/2022-07-05-00-04-55-test_data_2_dump.csv", header = None)
X_test3 = test_data3.iloc[:,:-1].values
y_test3 = test_data3.iloc[:,-1:].values

# Predict and classify test data
predictions_test3 = trained_model.predict(X_test3)
prediction_classes_test3_list = [
    1 if prob > 0.45 else 0 for prob in np.ravel(predictions_test3)
]
prediction_classes_test3 = np.array(prediction_classes_test3_list)

# Model evaluation
print("Test MSE: %.2f" % mean_squared_error(y_test3, prediction_classes_test3))
print('Test Accuracy: %.2f' %  accuracy_score(y_test3, prediction_classes_test3))
print('Test Confusion Matrix:', confusion_matrix(y_test3, prediction_classes_test3))
print('Test Precision: %.2f' % precision_score(y_test3, prediction_classes_test3))
print('Test Recall: %.2f' % recall_score(y_test3, prediction_classes_test3))


# Save predicted labels
file_name = "C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/trained_model/predicted_labels_2.csv"
np.savetxt(file_name, prediction_classes_test3, delimiter=",")

# Load test data
test_data4 = pd.read_csv("C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/test_data/2022-07-05-00-13-40-test_data_3_dump.csv", header = None)
X_test4 = test_data4.iloc[:,:-1].values
y_test4 = test_data4.iloc[:,-1:].values

# Predict and classify test data
predictions_test4 = trained_model.predict(X_test4)
prediction_classes_test4_list = [
    1 if prob > 0.45 else 0 for prob in np.ravel(predictions_test4)
]
prediction_classes_test4 = np.array(prediction_classes_test4_list)

# Model evaluation
print("Test MSE: %.2f" % mean_squared_error(y_test4, prediction_classes_test4))
print('Test Accuracy: %.2f' %  accuracy_score(y_test4, prediction_classes_test4))
print('Test Confusion Matrix:', confusion_matrix(y_test4, prediction_classes_test4))
print('Test Precision: %.2f' % precision_score(y_test4, prediction_classes_test4))
print('Test Recall: %.2f' % recall_score(y_test4, prediction_classes_test4))


# Save predicted labels
file_name = "C:/Users/Angelina/Documents/GitHub/WR-2-3-FFN-CLF/resources/trained_model/predicted_labels_3.csv"
np.savetxt(file_name, prediction_classes_test4, delimiter=",")

