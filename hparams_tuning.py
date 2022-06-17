# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 23:58:30 2022

@author: Valentin

Hyperparameter Tuning

"""

import tensorflow as tf
import pandas as pd
import datetime
from tensorboard.plugins.hparams import api as hp
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers

# Datensatz laden
train_data = pd.read_csv('2022-06-17-00-17-57-train_data_dump.csv', header=None)
X = train_data.iloc[:,:-1].values
y = train_data.iloc[:,-1:].values

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
    inputs = layers.Input(shape = (3362,), name = 'input_layer')
    model = tf.keras.Sequential(inputs)
    for i in range(hp.Int("n_layers", 1, 3)):
        model.add(layers.dense(units=hp.Int(f"units_{i}", 32, 128, step=32),
                               activation = 'sigmoid'
            ))
    model.add(layers.dense(1, activation = 'sigmoid'))
        
    model.compile(loss = tf.keras.losses.MeanSquaredError(),
                  ptimizer = tf.keras.optimizers.Adam(),
                  metrics = ["accuracy"])

clf.fit(X_train,
        y_train,
        epochs = 10,
        callbacks = tf.keras.callbacks.TensorBoard(log_dir),
        verbose = 0,
        validation_split = 0.2)
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
