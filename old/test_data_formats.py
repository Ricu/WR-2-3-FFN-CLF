# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 17:37:04 2022

@author: Valentin
"""
from timeit import default_timer as timer
import pandas as pd
csv_mat = pd.read_csv('test_matrix.csv', sep = ',', header = None)

# start = timer()
# parquet_mat = pd.read_parquet('test_matrix.parquet')
# end = timer()
# print(end-start)

# from scipy.io import loadmat
# import h5py
# start = timer()
# f = h5py.File('test_matrix.mat','r')
# mat_mat = f.get
# mat_mat = loadmat('test_matrix.mat')
# end = timer()
# print(end-start)

# import matlab.engine
# mat = matlab.engine.start_matlab()
# f = mat.load("dataset.mat", nargout=1)


