#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 13:19:27 2022

@author: kayla
"""
from sklearn import svm
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def run_SVM(data_file='Zs_hom.csv', label_file='hom_label.csv'):
    #load data
    #how to organize the data (I guess it should be all data, with the group label as 'target')
    #load Z_TNM_case.csv as np
    #create group label csv (and create cluster label csv) then add the label variable as 'target'
    data = pd.read_csv(data_file, header = None)
    label = pd.read_csv(label_file, header = None)
    #label = label.sample(frac=1).reset_index(drop=True)
    #label = pd.read_csv('cluster_label.csv')
    data['target'] = label
    
    #D = data.sample(frac=1)
    
    #split data into train/test splits
    train_dataset, test_dataset = train_test_split(data, test_size =.2)
    #test_dataset, valid_dataset = train_test_split(temp_test_dataset, test_size = 0.5)
    train_labels = train_dataset.pop('target')
    test_labels = test_dataset.pop('target')
    #valid_labels = valid_dataset.pop('target')
    
    #train model
    model = svm.SVC(C=1)
    #model = svm.LinearSVC (multi_c  lass='crammer_singer', max_iter=100000)
    model.fit(train_dataset, train_labels)
    #y_pred = model.predict(test_dataset)
    
    #evaluate
    from sklearn import metrics
    from sklearn.metrics import confusion_matrix
    y_predTrain = model.predict(train_dataset)
    print("Train Accuracy:", metrics.accuracy_score(train_labels, y_predTrain))
    cmTrain = confusion_matrix(train_labels, y_predTrain)
    #y_pred = model.predict(valid_dataset)
    #print("Valid Accuracy:", metrics.accuracy_score(valid_labels, y_pred))
    y_predTest = model.predict(test_dataset)
    print("Test Accuracy:", metrics.accuracy_score(test_labels, y_predTest))
    cmTest = confusion_matrix(test_labels, y_predTest)
    
    plot_conf_mat(cmTrain, 'Train Conference Matrix')
    plot_conf_mat(cmTest, 'Test Conference Matrix')

    
    #plot confusion matrix
def plot_conf_mat(cm, title):
    plt.figure()
    ax = plt.subplot()
    sns.heatmap(cm, annot=True, ax = ax)
    ax.set_xlabel('Predicted labels'); ax.set_ylabel('True labels');
    ax.set_title(title);
  
####okay something is wrong but what



def gen_sim_data(sim_name='sim_data.csv', label_name='sim_labels.csv'):
    dshape=(300,90)
    parlist=[[0,1], [15,2], [150, 30], [-5, 1], [-30,3], [-230, 43]]
    sim_data = np.concatenate([np.random.normal(i[0],i[1], size=dshape) for i in parlist],0)
    labels = []
    for i in range(6):
         labels += [i]*300
    labels = np.asarray(labels)
    
    np.savetxt(sim_name, sim_data, delimiter = ",")
    np.savetxt(label_name, labels, delimiter = "")
    return sim_name, label_name

#%% simulated data
sim_name, label_name = gen_sim_data()
run_SVM(data_file = sim_name, label_file = label_name)
    
#%%
run_SVM(data_file = 'Zs_sig.csv')