import os
from subprocess import PIPE, run
import time
import pandas as pd
import numpy as np 
import pprint
import requests


Data = pd.read_csv('/Users/sergiomares/Desktop/Nunez/Jin file/TSS_CpG_crispriphenotype_table.txt', delimiter = '\t',header = 0)
Data = Data.sort_values(["gene",'average phenotype of strongest 3'], ascending = True).drop_duplicates(subset = 'gene', keep = 'first')
UCSC_TSS = pd.read_csv('UCSC_TSS.txt', delimiter = '\t',header = 0)
UCSC_TSS = UCSC_TSS.drop_duplicates() 

errors = []

for gene in (Data['gene']):
        try:
                if len(UCSC_TSS[UCSC_TSS['hg19.kgXref.geneSymbol'] == gene]["hg19.knownGene.txStart"]) > 1:
                        x = abs(UCSC_TSS[UCSC_TSS['hg19.kgXref.geneSymbol'] == gene]["hg19.knownGene.txStart"] - Data[Data['gene'] == gene]["Primary TSS, 3'"].iloc[-1]).min()
                        Data.loc[Data['gene'] == gene, "Primary TSS, 3'"] = (int(Data[Data['gene'] == gene]["Primary TSS, 3'"].iloc[-1] - int(x) + 1))
                else:
                        x = UCSC_TSS[UCSC_TSS['hg19.kgXref.geneSymbol'] == gene]["hg19.knownGene.txStart"].iloc[-1]
                        Data.loc[Data['gene'] == gene, "Primary TSS, 3'"] = int(x) + 1 
              

        except:
            errors.append(gene)

len(errors)

Data = Data.reset_index()

Promoter_sequences = pd.read_csv('3kb_Promoter.sequences-2.csv', sep=',', header = 0)
Promoter_sequences = pd.merge(how = 'outer', left = Promoter_sequences, right = Data, left_on = 'Gene', right_on = 'gene')
Promoter_sequences = Promoter_sequences.dropna()

cpgs = pd.read_csv('1-s2-S0092867421003536-mmc3.csv',sep = ',',  header = 0)
cpgs = pd.DataFrame(cpgs)

df = pd.merge(how = 'outer', left = Promoter_sequences, right = cpgs, left_on = 'Gene', right_on = 'gene')
df = df.loc[df.Gene.notna()]
df = df.reset_index()

irbs = pd.read_csv('wgEncodeHaibMethylRrbsK562HaibSitesRep1.bed', sep='\t', header = 0)
#Extract only the values with 100% Certainty of methylation 
irbs = irbs[irbs['Unnamed: 8'] == '255,0,0']

chromosome_list = np.unique(irbs.track)

tmp = np.zeros((len(df),3001))

for index, chromosomes in enumerate(chromosome_list):
    
    table = irbs[irbs['track'] == chromosomes]
    table2 = df[df.chromosome == chromosomes]

    for i, x in enumerate(table['name="SL725.1']):


        for o, z in enumerate(table2["Primary TSS, 3'"]):

            if (x > (z - 1500) and (x < (z + 1500))) == True:
                y = x - (z - 1500)
                #print((z - 1500), z,  (z + 1500))
                #print('CpG Islands found in:', x, 'Position on vector:', y,"in row", Promoter_sequences['gene'][o], o + 1)    

                tmp[df.level_0[df["Primary TSS, 3'"] == z ].iloc[-1] + 1][int(y)] = 1

## Function to show each position for a certain gene

import matplotlib.pyplot as plt

#### Nilah recommended creating a plot for methylation near TSS (position and its phenotype score)

def methylated_basepairs(int):
    if len(np.where(tmp[int] == 1)[0]) > 0:
        is_c(int)
        is_c(int)

        x = df.head(int).tail(1)['Gene'].iloc[-1]
        # print(x, "| Coordinates of TSS:", df[df['Gene'] == x]["Primary TSS, 3'"].iloc[-1])
        # print("Found", len([np.where(tmp[int] == 1)][0][0]) , "Methylated basepairs at positions:")
        # print([np.where(tmp[int] == 1)][0][0] + df[df['Gene'] == x]["Primary TSS, 3'"].iloc[-1] - 1500)
        # plt.title("Methylation profile for promoter region of gene: {}".format(x), size = 15)
        # plt.xlabel('Position of the bases | 1500 = Start of TSS ')
        # plt.ylabel('Methylation')
        # plt.plot(tmp[int])

        # #print(df[df['Gene'] == x]["chromosome"].iloc[-1],':',[np.where(tmp[int] == 1)][0][0][0] + df[df['Gene'] == x]["Primary TSS, 3'"].iloc[-1] - 1500, '-', [np.where(tmp[int] == 1)][0][0][0] + df[df['Gene'] == x]["Primary TSS, 3'"].iloc[-1] - 1500 + 1)
        # print(df[df['Gene'] == x]["chromosome"].iloc[-1],':', df[df['Gene'] == x]["Primary TSS, 3'"].iloc[-1] - 1500, '-', df[df['Gene'] == x]["Primary TSS, 3'"].iloc[-1] + 1500 )

        #return([np.where(tmp[int] == 1)][0][0])

        return(len([np.where(tmp[int] == 1)][0][0] + df[df['Gene'] == x]["Primary TSS, 3'"].iloc[-1]))

    else:
        return 0

# Are the predicted methylation sites C?  

def is_c(int):
    
    counter = 0
    total = 0

    try:
        x = [np.where(tmp[int] == 1)][0][0]
        y = df.head(int).tail(1)['Gene'].iloc[-1]

        if len(x) == 0:
            #print("No Methylated Cs")
            counter = 2
            
        for i in x:

            total += 1

            if list(df.Sequences[df['Gene'] == y].iloc[-1])[i] == 'C':
                #print(i, 'is a C')
                counter +=1 
            else:
                #print(i, 'is not a C, its a ' + list(df.Sequences[df['Gene'] == y].iloc[-1])[i])
                tmp[int][i] = 0
    except:
        print("Issue with int", int, "Gene:")

    #print('Correct Cs:',counter/total*100,'%')


methylated_cs = [methylated_basepairs(i) for i in range(20340)]

from sklearn.preprocessing import LabelEncoder, OneHotEncoder

sequences = df.Sequences

 # This removes empty sequences.

# The LabelEncoder encodes a sequence of bases as a sequence of integers.
integer_encoder = LabelEncoder()  
# The OneHotEncoder converts an array of integers to a sparse matrix where 
# each row corresponds to one possible value of each feature.
one_hot_encoder = OneHotEncoder(categories='auto')   
input_features = []

counter = 0 
for sequence in sequences:
    try:
        integer_encoded = integer_encoder.fit_transform(list(sequence))
        integer_encoded = np.array(integer_encoded).reshape(-1, 1)
        one_hot_encoded = one_hot_encoder.fit_transform(integer_encoded)
        input_features.append(one_hot_encoded.toarray())
    except:
        print('error with', sequence)

indexes = []

for i in range(len(input_features)):
    try:
        if input_features[i].shape != (3001, 4):
            indexes.append(i)
            input_features.pop(i)
            print(i)
    except:
        print('error with', i)    

np.set_printoptions(threshold=40)
input_features = np.stack(input_features)
print("Example sequence\n-----------------------")
print('DNA Sequence #1:\n',sequences[0][:10],'...',sequences[0][-10:])
print('One hot encoding of Sequence #1:\n',input_features[0].T)

tmp_1 = np.delete(tmp, indexes, axis = 0)
md_input_features = np.dstack((input_features, tmp_1))
df = df.drop(index = indexes)

from sklearn.model_selection import train_test_split

train_features_1, test_features, train_labels, test_labels = train_test_split(
    md_input_features, df['CRISPRoff_average'], test_size=0.25, random_state=42)

from tensorflow.keras.layers import Dense, MaxPooling1D, Conv1D, GlobalMaxPool1D
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input

model_in = Input(shape=(train_features_1[1].shape))
model = (Conv1D(filters = 15, kernel_size = 5))(model_in)

model1 = (Conv1D(filters = 3, kernel_size = 2))(model)
# model1 = (MaxPooling1D(pool_size=(7)))(model1)
# model1 =(Dense(14, activation='relu'))(model1)
# model1 = (MaxPooling1D(pool_size=(7)))(model1)
model1 = (Dense(1, activation='relu'))(model1)

model = Model(model_in, model1)
model.compile(loss='binary_crossentropy',optimizer='rmsprop', metrics=['accuracy'])
model.summary()

history = model.fit(train_features_1, train_labels, 
                    epochs=1, verbose=0, validation_split=0.25)