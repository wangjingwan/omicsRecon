# clear system python, or cant use tf
# module purge
# module load GCC/11.2.0 Anaconda3/2022.05
# /home/wangjingwan/.conda/envs/tf/bin/python
import pandas as pd
import numpy as np
import os
import random as rn
rn.seed(12345)
np.random.seed(42)
os.environ['PYTHONHASHSEED'] = '0'

import tensorflow.compat.v1 as tf
from tensorflow.keras.layers import Input, Dense, Dropout
from tensorflow.keras.models import Model
from tensorflow.keras import regularizers

from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import argparse


def read_csv_tsv(filename):
    if 'csv' in filename:
        tmp = pd.read_csv(filename, sep = ',',header = 0,index_col=0)
    else:
        tmp = pd.read_csv(filename, sep = '\t',header = 0,index_col=0)
    return tmp


def rSquared(true,predicted):
    cols = predicted.shape[1]
    rsq = np.zeros(shape=(cols), dtype = np.float32)
    for j in range(cols):
        rsq[j] = r2_score(true[:,j], predicted[:,j])
    return rsq


def pear(D,D_re):
    D = np.array(D)
    D_re = np.array(D_re)
    tmp = np.corrcoef(D.flatten(order='C'), D_re.flatten(order='C'))
    return tmp[0,1] 

######################################################
##################### 1.env prep #####################
######################################################
# The below tf.set_random_seed() will make random number generation
# in the TensorFlow backend have a well-defined initial state.
tf.disable_v2_behavior()
tf.set_random_seed(1234)
# Force TensorFlow to use single thread.
session_conf = tf.ConfigProto(intra_op_parallelism_threads=1,
                              inter_op_parallelism_threads=1,
                              allow_soft_placement=True)
sess = tf.compat.v1.Session(graph=tf.get_default_graph(), config=session_conf)
tf.compat.v1.keras.backend.set_session(sess)

######################################################
##################### 2.data prep ####################
######################################################
inputDir = '/data6/wangjingwan/0.omics/sprout_omics-main/tutorial/examples/methy/input/EPD/zero_std1/'
preprocessed_DNAMeth = read_csv_tsv(f'{inputDir}/train_DNAMeth.csv')
preprocessed_RNASeq = read_csv_tsv(f'{inputDir}/train_RNASeq.csv')
check_DNAMeth = read_csv_tsv(f'{inputDir}/test_DNAMeth.csv')
check_RNASeq = read_csv_tsv(f'{inputDir}/test_RNASeq.csv')
labels = read_csv_tsv('/data6/wangjingwan/0.omics/sprout_omics-main/tutorial/examples/methy/data/train/meta_NMT.tsv')
labels = labels['ensembl']
labels = labels.loc[preprocessed_DNAMeth.index]
outDir = inputDir




preprocessed_CNA  = pd.DataFrame(np.ones((preprocessed_RNASeq.shape))*0.51,columns=preprocessed_RNASeq.columns, index = preprocessed_RNASeq.index)
check_CNA = pd.DataFrame(np.ones((check_RNASeq.shape))*0.51,columns=check_RNASeq.columns, index = check_RNASeq.index)
z = pd.concat((check_DNAMeth, check_CNA),axis=1)
scalar = MinMaxScaler()
testset = scalar.fit_transform(z)
testy = check_RNASeq

######################################################
################### 3.start traning ##################
######################################################
x1 = preprocessed_DNAMeth
x2 = preprocessed_CNA
y = preprocessed_RNASeq

x1 = pd.DataFrame(x1)
x2 = pd.DataFrame(x2)
z = pd.concat((x1, x2),axis=1)

x_train, x_test, y_train, y_test, labels_train, labels_test = train_test_split(z, y, labels, test_size=0.1)

scalar = MinMaxScaler()
x_train = scalar.fit_transform(x_train)
x_test = scalar.transform(x_test)

noise_factor = 0.5
x_train_noisy = x_train + noise_factor * np.random.normal(0.0, 1.0, x_train.shape)
x_test_noisy = x_test + noise_factor * np.random.normal(0.0, 1.0, x_test.shape)

x_train_noisy = np.clip(x_train_noisy, 0., 1.)
x_test_noisy = np.clip(x_test_noisy, 0., 1.)

num_in_neurons = z.shape[1]
num_out_neurons = y.shape[1]

batch = 8
#batch = x_train_noisy.shape[0]
#with tf.device('/gpu:0'):
    # this is the size of our enconded representations
encoding_dim1 = 500
encoding_dim2 = 200
lambda_act = 0.0001
lambda_weight = 0.001
# this is our input placeholder
input_data = Input(shape=(num_in_neurons,))
# first encoded representation of the input
encoded = Dense(encoding_dim1, activation='relu', activity_regularizer=regularizers.l1(lambda_act), kernel_regularizer=regularizers.l2(lambda_weight), name='encoder1')(input_data)
# second encoded representation of the input
encoded = Dense(encoding_dim2, activation='relu', activity_regularizer=regularizers.l1(lambda_act), kernel_regularizer=regularizers.l2(lambda_weight), name='encoder2')(encoded)
# first lossy reconstruction of the input
decoded = Dense(encoding_dim1, activation='relu', name='decoder1')(encoded)
# the final lossy reconstruction of the input
decoded = Dense(num_in_neurons, activation='sigmoid', name='decoder2')(decoded)
# this model maps an input to its reconstruction
autoencoder = Model(inputs=input_data, outputs=decoded)
myencoder = Model(inputs=input_data, outputs=encoded)
autoencoder.compile(optimizer='sgd', loss='mse')
# training
print('training the autoencoder')
autoencoder.fit(x_train_noisy, x_train,
                epochs=25,
                batch_size=batch,
                shuffle=True,
                validation_data=(x_test_noisy, x_test))
autoencoder.trainable = False   #freeze autoencoder weights


num_hidden = encoding_dim2
# with tf.device('/gpu:0'):
x = autoencoder.get_layer('encoder2').output
x = Dropout(0.2)(x)             # adding 20% dropout
h = Dense(int(num_hidden * 3), activation='relu', name='hidden1')(x)
h = Dropout(0.5)(h)             # adding 50% dropout
h = Dense(int(num_hidden * 5), activation='relu', name='hidden2')(h)
h = Dropout(0.5)(h)             # adding 50% dropout
y = Dense(num_out_neurons, activation='linear', name='prediction')(h)
mlpRegressor = Model(inputs=autoencoder.inputs, outputs=y)
# Compile model
mlpRegressor.compile(loss='mse', optimizer='adam', metrics=['accuracy'])    # or loss='mae'
# Fit the model
print('training the MLP multi-output regressor')
mlpRegressor.fit(x_train, y_train, epochs=50, batch_size=batch)
y_pred = mlpRegressor.predict(x_test)
actual_mean = pd.DataFrame(y_test.mean(axis=0))
pred_mean = pd.DataFrame(y_pred.mean(axis=0))


print('MSE: (Actual Vs. Predicted)', mean_squared_error(y_test, y_pred))
print('r^2 value: (Mean of actual Vs. Mean of Predicted)', r2_score(actual_mean, pred_mean))


######################################################
################### 3.start testing ##################
######################################################
#with tf.device('/gpu:0'):
x = autoencoder.get_layer('encoder2').output
x = Dropout(0.2)(x)             # adding 20% dropout
h = Dense(int(num_hidden * 3), activation='relu', name='hidden1')(x)
h = Dropout(0.5)(h)             # adding 50% dropout
h = Dense(int(num_hidden * 5), activation='relu', name='hidden2')(h)
h = Dropout(0.5)(h)             # adding 50% dropout
y = Dense(num_out_neurons, activation='linear', name='prediction')(h)
mlpRegressor = Model(inputs=autoencoder.inputs, outputs=y)
# Compile model
mlpRegressor.compile(loss='mse', optimizer='adam', metrics=['accuracy'])    # or loss='mae'
# Fit the model
print('training the MLP multi-output regressor')
mlpRegressor.fit(x_train, y_train, epochs=50, batch_size=batch)
x_test = np.array(testset)
y_test = np.array(testy)
y_pred = mlpRegressor.predict(x_test)
actual_mean = pd.DataFrame(y_test.mean(axis=0))
pred_mean = pd.DataFrame(y_pred.mean(axis=0))


mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(actual_mean, pred_mean)
p = pear(np.array(check_RNASeq),y_pred)
res = pd.DataFrame()
for cell in result.index.tolist():
    cor = pear(result.loc[cell], check_RNASeq.loc[cell])
    tmp = pd.DataFrame(np.array([[cell, cor]]),index = ['zero_std1'])
    res = pd.concat((res, tmp))
    
evals = pd.DataFrame(np.array([[mse,r2,p]]),index = [outDir],columns=['mse','r2_mean','pear'])
evals.to_csv(f'{outDir}/methy_pred_evals.csv',sep = ',',header=True,index=True)

result = pd.DataFrame(y_pred,index = check_RNASeq.index, columns= check_RNASeq.columns)
result.to_csv(f'{outDir}/methy_pred_rna.csv',sep = ',',header=True,index=True)


#np.save(f'{args.outDir}/Embryo_pred_rna.npy',y_pred)
old = read_csv_tsv('/data6/wangjingwan/0.omics/sprout_omics-main/tutorial/examples/methy/input/cell_cor.csv')
draw = read_csv_tsv('compare2up1000.csv')
import seaborn as sns
import matplotlib.pyplot as plt
plt.figure(figsize=(3,4))
sns.boxplot(x="index", y='1', data=draw)
plt.xticks(rotation=90)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.xlabel('',fontsize=16)
plt.ylabel('',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# plt.axis('equal')
plt.tight_layout()
#plt.savefig('/data6/wangjingwan/0.omics/sprout_omics-main/tutorial/examples/methy/input/cell_cor.pdf')
plt.show()