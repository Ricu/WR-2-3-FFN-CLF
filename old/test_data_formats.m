% TEST DATA FORMATS
clc; clear;
matrix = rand(100000,243);
fprintf('Matrix der Groesse %i,%i\n',size(matrix,1),size(matrix,2));
%%
f_mat = @() save('test_matrix.mat');
time_mat = timeit(f_mat);
fprintf('Dauer fuer .mat Format mit save(): %e\n', time_mat)
fmat = dir(which('test_matrix.mat'));
size_mat = fmat.bytes;
mat = [time_mat;size_mat];
%%
f_csv = @() writematrix(matrix,'test_matrix.csv');
time_csv = timeit(f_csv);
fprintf('Dauer fuer .csv Format mit writematrix(): %e\n', time_csv)
fcsv = dir(which('test_matrix.csv'));
size_csv = fcsv.bytes;
csv = [time_csv;size_csv];
%%
f_dat = @() writematrix(matrix,'test_matrix.dat');
time_dat = timeit(f_dat);
fprintf('Dauer fuer .dat Format mit writematrix(): %e \n', time_dat)
fdat = dir(which('test_matrix.dat'));
size_dat = fdat.bytes;
dat = [time_dat;size_dat];
%%
f_txt = @() writematrix(matrix,'test_matrix.txt');
time_txt = timeit(f_txt);
fprintf('Dauer fuer .txt Format mit writematrix(): %e\n', time_txt)
ftxt = dir(which('test_matrix.txt'));
size_txt = ftxt.bytes;
txt = [time_txt;size_txt];
%%
table_mat = array2table(matrix);
f_table = @() array2table(matrix);
time_table = timeit(f_table);
f_par = @() parquetwrite("test_matrix.parquet", table_mat);
time_par = timeit(f_par);
fprintf('Dauer fuer .parquet Format mit parquetwrite(): %e\n', time_par+time_table)
fpar = dir(which('test_matrix.parquet'));
size_par = fpar.bytes;
parquet = [time_par;size_par];


t = table({'runtime';'size'},mat,csv,dat,txt,parquet);
display(t)

