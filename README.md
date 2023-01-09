# Sleep-disorder-identification-using-physiological-signal
Note: This project was executed on MATLAB 2020.

%%%%% Pre-Requisites %%%%%

This repository contains MATLAB codes of data preparation for binary classification of normal and disordered patients. Before using these codes download the CAP sleep data from physionet. Link to data: https://archive.physionet.org/cgi-bin/atm/ATM

Since the data is in the form of EDF format, the first step will be to download the script file to read EDF documents. The same physionet provides the script for that. After downloading, change the file name to "edfread".

change the current directory path in the ScoreReader_general file.


%%%%% General Explanation %%%%%

Physionet data contains two files: 1) Label / Annotation file and 2) Raw data/ EDF file

The annotation will be opened using the ScoreReader_general file which contains the sleep stage information of a subject. EDF file contains data acquired through various channels. Channel information for the CAP sleep database is given in Channels.mat.

We have done manual data segregation using sleep stage labels. After that, we applied wavelet decomposition of 7 levels to the data and calculated Hjorth parameter-based features for all subbands. For wavelet decomposition, we used a bi-orthogonal filter bank (file name: myExample - 5).

Run normal and all disorder files (file name starts with darji) to calculate the features of all the subjects. 

%%%%% Classification %%%%%

The last step is to feed all the features to classification algorithms. we obtained the best results with ensemble bagged trees.


For more information on results, read our paper at: https://www.sciencedirect.com/science/article/abs/pii/S0010482522000166

