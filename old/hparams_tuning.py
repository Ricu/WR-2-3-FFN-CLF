# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 23:58:30 2022

@author: Valentin

Hyperparameter Tuning

"""
import numpy as np
import tensorflow as tf
import keras_tuner
import pandas as pd
import datetime
import os
# from tensorboard.plugins.hparams import api as hp
from  tensorflow import keras
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers
from kerastuner import RandomSearch
from kerastuner.engine.hyperparameters import HyperParameters


# Datensatz laden
train_data = pd.read_csv('2022-06-28-11-05-08-train_data_dump.csv', header=None)
coefficients = train_data.iloc[:,:-1].values
label = train_data.iloc[:,-1:].values
print(label)

# Splitte Daten in Test und Trainingsdatensatz auf
X_train, X_test, y_train, y_test = train_test_split(X,
                                                    y,
                                                    test_size=0.25,
                                                    random_state=42)

log_dir = "logs/hparam_tuning/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

# Preprocessing (normalisieren)
normalize = layers.Normalization().adapt(X_train)


def build_model(hp):
# Simples Modell
    inputs = keras.Input(shape = (3362,), name = 'input_layer')
    model = keras.Sequential(inputs)
    for i in range(hp.Int("n_layers", 1, 3)):
        model.add(layers.dense(units=hp.Int(f"units_{i}", min_value=32, max_value=128, step=32),
                               activation = 'sigmoid'
            ))
    model.add(layers.dense(1, activation = 'sigmoid'))
        
    model.compile(loss = keras.losses.MeanSquaredError(),
                  optimizer = keras.optimizers.Adam(hp.Choice('learning_rate',values=[1e-2, 1e-3])),
                  metrics = ["accuracy"])
    return model


#hp = keras_tuner.HyperParameters()

# Random Search to find best hyperparameters
tuner_search = RandomSearch(build_model,
                            objective = 'val_accuracy',
                            max_trials = 10, directory = 'output', project_name = "sortingEdges")
tuner_search.search(train_edges,train_labels,epochs=10,validation_split=0.2)

# Best Hyperparameters Determined
model = tuner_search.get_best_models(num_models=1)[0]
model.summary()

# Retraining considering the best model hyperparameters
history = model.fit(train_edges,
                    train_labels,
                    epochs = 10,
                    callbacks = tf.keras.callbacks.TensorBoard(log_dir),
                    #verbose = 0,
                    validation_split = 0.2)

# Make predicitions
predictions = model.predict(test_edges)

accuracy = clf.evaluate(X_test,y_test)



























##### HParam Tuning
# Hyperparametr Grid
HP_NUM_UNITS = hp.HParam('num_units', hp.Discrete([16, 32]))
HP_DROPOUT = hp.HParam('dropout', hp.RealInterval(0.1, 0.2))
HP_OPTIMIZER = hp.HParam('optimizer', hp.Discrete(['adam', 'sgd']))

METRIC_ACCURACY = 'accuracy'

with tf.summary.create_file_writer(log_dir).as_default():
  hp.hparams_config(
    hparams=[HP_NUM_UNITS, HP_DROPOUT, HP_OPTIMIZER],
    metrics=[hp.Metric(METRIC_ACCURACY, display_name='Accuracy')],
  )

def train_test_model(hparams):
  model = tf.keras.models.Sequential([
    tf.keras.layers.Flatten(),
    tf.keras.layers.Dense(hparams[HP_NUM_UNITS], activation=tf.nn.relu),
    tf.keras.layers.Dropout(hparams[HP_DROPOUT]),
    tf.keras.layers.Dense(10, activation=tf.nn.softmax),
  ])
  model.compile(
      optimizer=hparams[HP_OPTIMIZER],
      loss='sparse_categorical_crossentropy',
      metrics=['accuracy'],
  )

  model.fit(X_train,
            y_train,
            epochs=1,
            callbacks=[tf.keras.callbacks.TensorBoard(log_dir),  # log metrics
                       hp.KerasCallback(log_dir, hparams),  # log hparams
                       ]
            )
  accuracy = model.evaluate(X_test, y_test)
  return accuracy

def run(run_dir, hparams):
  with tf.summary.create_file_writer(run_dir).as_default():
    hp.hparams(hparams)  # record the values used in this trial
    accuracy = train_test_model(hparams)
    tf.summary.scalar(METRIC_ACCURACY, accuracy, step=1)


session_num = 0

for num_units in HP_NUM_UNITS.domain.values:
  for dropout_rate in (HP_DROPOUT.domain.min_value, HP_DROPOUT.domain.max_value):
    for optimizer in HP_OPTIMIZER.domain.values:
      hparams = {
          HP_NUM_UNITS: num_units,
          HP_DROPOUT: dropout_rate,
          HP_OPTIMIZER: optimizer,
      }
      run_name = "run-%d" % session_num
      print('--- Starting trial: %s' % run_name)
      print({h.name: hparams[h] for h in hparams})
      run('logs/hparam_tuning/' + run_name, hparams)
      session_num += 1
