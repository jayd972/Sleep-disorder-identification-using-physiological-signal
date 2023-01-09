# Sleep-disorder-identification-using-physiological-signal
%%%%% Pre-Requisites %%%%%

This repository contains matlab codes of data preparation for binary classification of normal and disordered patients. Before using these codes download the CAP sleep data from physionet. Link to data: https://archive.physionet.org/cgi-bin/atm/ATM

Since the data is in the form of EDF format, first step will be to download the script file to read EDF documents. The same physionet provides the script for that. After downloading, change the file name to "edfread".

change the current directory path in ScoreReader_general file.


%%%%% General Explaination %%%%%

Physionet data contains two files : 1) Label / Annotation file and 2) Raw data/ EDF file

Annotation will be open using ScoreReader_general file which contains sleep stage information of a subjct. EDF file contains data aquired through various channels. Channel information for the CAP sleep database is given in Channels.mat.

We have done manual data saggerigation using sleep stage labels. After that we applied wavelet decompotion of 7 level to the data and calculated hjorth parameter based features for all sub bands. For wavelet decomposition we used bi orthogonal filterbank (file name: myExample - 5).

Run normal and all disorder files (file name starts with darji) to calculate the features of all the subjects. 

%%%%% Classification %%%%%

Last step is to feed all the features to classification algorithms. we obtained best results with ensembal bagged trees.


For more information on results, read our paper on : https://www.sciencedirect.com/science/article/abs/pii/S0010482522000166
