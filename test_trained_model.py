# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:42:09 2022

@author: Angelina
"""
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, make_scorer, r2_score
from sklearn.metrics import accuracy_score
from sklearn.externals import joblib

# Load trained model
joblib_file = "C:\Users\Angelina\Documents\GitHub\WR-2-3-FFN-CLF\resources\trained_model\trained_model.pkl"
trained_model = joblib.load(joblib_file)


# Load test data
test_data = pd.read_csv("C:\Users\Angelina\Documents\GitHub\WR-2-3-FFN-CLF\resources\test_data\2022-06-24-13-10-55-test_data_1_dump.csv")
coeff_test = test_data.iloc[:,:-1].values
label_test = test_data.iloc[:,-1:].values

# Calculate the predictions
label_predict = trained_model.predict(coeff_test)

# Berechne Fehler und Accuracy
mse = mean_squared_error(label_predict,label_test)
acc = accuracy_score(label_predict,label_test)

print("MSE: %.2f" % mse)
print("ACC: %.2f" % acc)

# Save predicted labels
fname = "C:\Users\Angelina\Documents\GitHub\WR-2-3-FFN-CLF\resources\trained_model\predicted_labels.csv"
np.savetxt(fname, label_predict, delimiter=",")

