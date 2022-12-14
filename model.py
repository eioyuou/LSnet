#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
import pysam
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Input, Dense, Conv2D, DepthwiseConv2D, BatchNormalization, Dropout, GlobalAveragePooling2D, Reshape, multiply, add, Activation
#tf.config.set_visible_devices([], 'GPU')
tf.random.set_seed(123)
from os import listdir
from os.path import isfile, join
from tensorflow import keras
import tensorflow_addons as tfa
import math
from tensorflow.keras.utils import Sequence
from tensorflow.keras import backend as K
from collections import Counter
import os
import time


# In[2]:


def cbam_block(cbam_feature, ratio=8,kernel_size = (1,9)):

    cbam_feature = channel_attention(cbam_feature, ratio)
    cbam_feature = spatial_attention(cbam_feature,kernel_size)
    return cbam_feature

def channel_attention(input_feature, ratio=8):

    #channel_axis = 1 if K.image_data_format() == "channels_first" else -1
    channel = input_feature.shape[-1]
    filters = max(1, int(channel//ratio))
    shared_layer_one = tf.keras.layers.Dense(filters,
                             activation='relu',
                             kernel_initializer='he_normal',
                             use_bias=True,
                             bias_initializer='zeros')
    shared_layer_two = tf.keras.layers.Dense(channel,
                             kernel_initializer='he_normal',
                             use_bias=True,
                             bias_initializer='zeros')

    avg_pool = tf.keras.layers.GlobalAveragePooling2D()(input_feature)    
    avg_pool = tf.keras.layers.Reshape((1,1,channel))(avg_pool)
    avg_pool = shared_layer_one(avg_pool)
    avg_pool = shared_layer_two(avg_pool)

    max_pool = tf.keras.layers.GlobalMaxPooling2D()(input_feature)
    max_pool = tf.keras.layers.Reshape((1,1,channel))(max_pool)
    max_pool = shared_layer_one(max_pool)
    max_pool = shared_layer_two(max_pool)
   

    cbam_feature = tf.keras.layers.Add()([avg_pool,max_pool])
    cbam_feature = tf.keras.layers.Activation('sigmoid')(cbam_feature)


    return multiply([input_feature, cbam_feature])
def spatial_attention(input_feature,kernel_siz):
    kernel_size = kernel_siz

    channel = input_feature.shape[-1]
    cbam_feature = input_feature

    avg_pool = tf.keras.layers.Lambda(lambda x: K.mean(x, axis=3, keepdims=True))(cbam_feature)
    #assert avg_pool._keras_shape[-1] == 1
    max_pool = tf.keras.layers.Lambda(lambda x: K.max(x, axis=3, keepdims=True))(cbam_feature)
    #assert max_pool._keras_shape[-1] == 1
    concat = tf.keras.layers.Concatenate(axis=3)([avg_pool, max_pool])
    #assert concat._keras_shape[-1] == 2
    cbam_feature = tf.keras.layers.Conv2D(filters = 1,
                    kernel_size=kernel_size,
                    strides=1,
                    padding='same',
                    activation='sigmoid',
                    kernel_initializer='he_normal',
                    use_bias=False)(concat)	
    #assert cbam_feature._keras_shape[-1] == 1

    return multiply([input_feature, cbam_feature])




# In[3]:


def cnn_model():
    inputs = tf.keras.Input(shape=(200,9,1))
    x = tf.keras.layers.Conv2D(128, kernel_size =(2,9), padding='same', activation='relu')(inputs)
    x = tf.keras.layers.MaxPool2D((2,1))((x))
    x = tf.keras.layers.Conv2D(64, kernel_size =(2,1), padding='same', activation='relu')(x)
    x = tf.keras.layers.MaxPool2D((2,1))(x)
    x = cbam_block(x)
    x = tf.keras.layers.Conv2D(64, kernel_size =(2,1), padding='same', activation='relu')(x)
    x = tf.keras.layers.MaxPool2D((2,1))(x)

    x = tf.keras.layers.Conv2D(64, kernel_size =(2,1), padding='same', activation='relu')(x)
    x = tf.keras.layers.MaxPool2D((2,1))(x)
    x = cbam_block(x)
    x = tf.keras.layers.Conv2D(64, kernel_size =(2,1), padding='same', activation='relu')(x)
    x = tf.keras.layers.MaxPool2D((2,1))(x)

    x = tf.keras.layers.Conv2D(64, kernel_size =(2,1), padding='same', activation='relu')(x)
    x = tf.keras.layers.MaxPool2D((2,1))(x)
    x = cbam_block(x)
    x = tf.keras.layers.Conv2D(64, kernel_size =(2,1), padding='same', activation='relu')(x)
    x = tf.keras.layers.MaxPool2D((2,1))(x)
    x = tf.keras.layers.Flatten()(x)
    model = tf.keras.Model(inputs, x)
    return model


# In[15]:


def gru_model():
    inputs = tf.keras.Input((None, 200, 9,1))
    #inputs_x = tf.keras.layers.Reshape((None, 200, 6,1))(inputs)
    cnn_layer_object = cnn_model()
    encoded_frames = tf.keras.layers.TimeDistributed(cnn_layer_object)(inputs)
    encoded_sequence = tf.keras.layers.Bidirectional(tf.keras.layers.GRU(64, return_sequences = True
                                                               ))(encoded_frames)
    
    encoded_sequence = tf.keras.layers.Bidirectional(tf.keras.layers.GRU(64, return_sequences = True
                                                               ))(encoded_sequence)
    hidden_layer = tf.keras.layers.Dense(units=128, activation="relu",)(encoded_sequence)
    hidden_layer = tf.keras.layers.Dropout(0.4)(hidden_layer)
    hidden_layer = tf.keras.layers.Dense(units=128, activation="relu",)(hidden_layer)
    hidden_layer = tf.keras.layers.Dropout(0.4)(hidden_layer)
    outputs = tf.keras.layers.Dense(units=1, activation="sigmoid")(hidden_layer)
    model = tf.keras.Model(inputs, outputs)
    
  
    
    model.compile(loss=tf.keras.losses.BinaryCrossentropy(label_smoothing=0.01), 
            optimizer=tf.keras.optimizers.Adam(learning_rate = 3e-4),
            metrics=[keras.metrics.TruePositives(name='tp'),
              keras.metrics.FalsePositives(name='fp'),
              keras.metrics.TrueNegatives(name='tn'),
              keras.metrics.FalseNegatives(name='fn'), 
              keras.metrics.BinaryAccuracy(name='accuracy'),
              keras.metrics.Precision(name='precision'),
              keras.metrics.Recall(name='recall')])
    model.summary()
    return model


# In[16]:


def decode_flag(Flag):
    signal = {1 << 2: 0, 1 >> 1: 1, 1 << 4: 2, 1 << 11: 3, 1 << 4 | 1 << 11: 4}
    return signal[Flag] if(Flag in signal) else 0

def c_pos(cigar, refstart):
    number = ''
    numlist = [str(i) for i in range(10)]
    readstart = False
    readend = False
    refend = False
    readloc = 0
    refloc = refstart
    for c in cigar:
        if(c in numlist):
            number += c
        else:
            number = int(number)
            if(readstart == False and c in ['M', 'I', '=', 'X']):
                readstart = readloc
            if(readstart != False and c in ['H', 'S']):
                readend = readloc
                refend = refloc
                break

            if(c in ['M', 'I', 'S', '=', 'X']):
                readloc += number

            if(c in ['M', 'D', 'N', '=', 'X']):
                refloc += number
            number = ''
    if(readend == False):
        readend = readloc
        refend = refloc

    return refstart, refend, readstart, readend 






def del_signature(pre,bamfile):

    data = []
    for chr_name,start,end in pre:

        for read in bamfile.fetch(chr_name,start,end):
            aligned_length = read.reference_length
            if aligned_length == None:
                aligned_length= 0
            if (read.mapping_quality >= 0) and aligned_length >= 0:
                cigar = []
                sta = read.reference_start
                for ci  in read.cigartuples:
                    
                    #print(start)
                    if ci[0] in [0,3, 7, 8]:
                        sta += ci[1]
                    elif ci[0] == 2 :
                       #print(ci[1])
                        if ci[1] >=30 and (abs(sta-start) <200):
                            cigar.append([sta,sta+ci[1],ci[1]])
                        sta += ci[1]
                if len(cigar) !=0:
                    cigar = np.array(cigar)
                    cigar = cigar[np.argsort(cigar[:,0])]
                    a = mergecigar(cigar)

                    data.extend(a)
            #print(len(data))
            if(read.has_tag('SA') == True):
                code = decode_flag(read.flag)
                sapresent = True
                rawsalist = read.get_tag('SA').split(';')
                #print(rawsalist)
                for sa in rawsalist[:-1]:
                    sainfo = sa.split(',')
                    #print(sainfo)
                    tmpcontig, tmprefstart, strand, cigar = sainfo[0], int(sainfo[1]), sainfo[2], sainfo[3]
                    if(tmpcontig != chr_name):
                        continue
                    #print(code,strand)   
                    if((strand == '-' and (code %2) ==0) or (strand == '+' and (code %2) ==1)):
                        refstart_1, refend_1, readstart_1, readend_1 =  read.reference_start, read.reference_end,read.query_alignment_start,read.query_alignment_end
                        refstart_2, refend_2, readstart_2, readend_2 = c_pos(cigar, tmprefstart)
                        a = readend_1 - readstart_2
                        b = refend_1 - refstart_2
                        if(abs(a)<500):
                            if(abs(b-a)<30):
                                continue
                            if((b-a)<0):

                                if(min(refend_1, refstart_2)<(end) and min(refend_1, refstart_2)>(start)):
                                    #print(b-a)
                                    #print(start,end)
                                    data1 = [min(refend_1, refstart_2), min(refend_1, refstart_2)+abs(b-a),abs(b-a)]
                                    data1 = np.array(data1)
                                    #print(data1)
                                    #print(data1)
                                    data.append(data1)
    #print(data)
            #print(len(data))
      
    data = np.array(data)
    #print(len(data))  
    #print(data)
    if len(data) == 0:
        return data
    data = data[np.argsort(data[:,0])]
            
                
    return data


def mergecigar(infor):
    data = []
    i = 0
    while i>=0:
        count = 0
        if i >(len(infor)-1):
            break
        for j in range(i+1,len(infor)):
            #print(i,j)
            if abs(infor[j][0] - infor[i][1]) <= 40: 
                count = count + 1
                infor [i][1] = infor[j][1]
        #print(infor[i],count)
        lenth = abs(infor[i][0] - infor[i][1])
        #print(infor)
        data.append([infor[i][0],infor[i][1], lenth])


        if count == 0:
            i += 1
        else :
            i += (count+1)
    return data

def merge(infor):
    data = []
    i = 0
    while i>=0:
        dat = []
        
        count = 0
        if i >(len(infor)-1):
            break
        dat.append(infor[i])
        for j in range(i+1,len(infor)):
            #print(i,j)
            if( (abs(infor[i][0] -infor[j][0]) <= 1500) and (abs(infor[i][1] - infor[j][1])<= 1500)):
                count = count + 1
                dat.append(infor[j])
        #print(infor[i],count)
        dat = np.array(dat)
        data.append(dat)


        if count == 0:
            i += 1
        else :
            i += (count+1)
    return data


def mergedeleton_long(pre,index,chr_name,bamfile):
    data = []
    deletion = []
    for i in range(len(pre)):
        if pre[i] > 0.5:
            data.append([chr_name,index[i],index[i]+200])
    signature = del_signature(data,bamfile)
    #print(signature)
    dell = merge(signature)
    #print(dell)""
    for sig in dell:
            pp = np.array(sig)
            start = math.ceil(np.median(pp[:,0]))
            #en = math.ceil(pp[:1].mean)
            length = math.ceil(np.median(pp[:,2]))
            deletion.append([chr_name,start,length,len(pp),'DEL'])
    return deletion 


# In[17]:


def geno_model():
    inputs = tf.keras.Input(shape=(200,9,1))
    
    x = (Conv2D(64,kernel_size=(2,9),padding='same', activation = 'relu'))(inputs)
    x = (tf.keras.layers.MaxPool2D((2, 1)))(x)
    x = cbam_block(x)
    
    x = (Conv2D(64,kernel_size=(2,1),padding='same', activation = 'elu'))(x)
    x = (tf.keras.layers.MaxPool2D((2, 1)))(x)

    x = (Conv2D(64,kernel_size=(2,1),padding='same', activation = 'elu'))(x)
    x = (tf.keras.layers.MaxPool2D((2, 1)))(x)
    x = cbam_block(x)
        
    x = (Conv2D(64,kernel_size=(2,1),padding='same', activation = 'elu'))(x)
    x = (tf.keras.layers.MaxPool2D((2, 1)))(x)

    x = (Conv2D(64,kernel_size=(2,1),padding='same', activation = 'elu'))(x)
    x = (tf.keras.layers.MaxPool2D((2, 1)))(x)
    x = cbam_block(x)
    
    x = (Conv2D(64,kernel_size=(2,1),padding='same', activation = 'elu'))(x)
    x = (tf.keras.layers.MaxPool2D((2, 1)))(x)
    x = cbam_block(x)
    x = tf.keras.layers.Flatten()(x)
    x = tf.keras.layers.Dense(units=128, activation="relu",)(x)
    x = tf.keras.layers.Dropout(0.4)(x)
    x = tf.keras.layers.Dense(units=128, activation="relu",)(x)
    x = tf.keras.layers.Dropout(0.4)(x)
    outputs = tf.keras.layers.Dense(units=1,activation="sigmoid")(x)
    model = tf.keras.Model(inputs, outputs)
        
    model.compile(loss=tf.keras.losses.BinaryCrossentropy(label_smoothing=0.01), 
            optimizer=tf.keras.optimizers.Adam(learning_rate = 3e-4),
            metrics=[keras.metrics.TruePositives(name='tp'),
              keras.metrics.FalsePositives(name='fp'),
              keras.metrics.TrueNegatives(name='tn'),
              keras.metrics.FalseNegatives(name='fn'), 
              keras.metrics.BinaryAccuracy(name='accuracy'),
              keras.metrics.Precision(name='precision'),
              keras.metrics.Recall(name='recall')])


    return model


# In[18]:


def find_geno(dell,b):
    count = 0
    data = []
    geno_count = []
    ti = time.time()
    for i in dell[:]:
        contig,start,end,support,svlen = str(i[0]),int(i[1]),(int(i[2])),10,int(i[5])

        start_ge = (start//200)*200
        end_ge = (end//200+1)*200

        c = b[b[:,0]==contig]
        w = np.where((np.array(c[:,1],dtype=int)>= start_ge)&(np.array(c[:,2],dtype=int)<= end_ge) )[0]

        geno_count = np.array(c[w][:,-1],dtype = int)

        #geno_count = np.array(geno_count)

        if (len(geno_count[geno_count==1]) >= len(geno_count[geno_count==0])):
            data.append([contig,start,svlen,support,'DEL','1/1'])
        #elif (len(geno_count[geno_count==0]) >= len(geno_count[geno_count==1]) ):
        else:
            data.append([contig,start,svlen,support,'DEL','0/1'])
    return data


# In[19]:


def tovcf(data, contig2length, outputpath):
    top = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">\n"""
    body = ''
    for contig in contig2length:
        body += "##contig=<ID="+contig+",length="+str(int(contig2length[contig]))+">\n"
    tail = """##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, INS=Insertion">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of read support this record">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t.\n"""

    myvcf = top+body+tail
    #genomapper = {0:'0/1', 1:'1/1'}
    #svlist = [[rec[0], int(rec[1]), int(rec[2]), int(rec[3]), rec[-1]]]
        
    for rec in data:


        contig = rec[0]


        geno = rec[-1]
        
        recinfo = 'SVLEN=' + str(int(rec[2]))+';SVTYPE=' + 'DEL'+';END='+str(int(rec[1])+abs(int(rec[2])))+';RE='+str(int(rec[3])+1)+'\tGT\t'+str(geno)+'\n'
        myvcf += (contig +'\t'+ str(int(rec[1]))+'\t'+ '.'+'\t'+ '.'+'\t'+ '.'+'\t'+ str(int(rec[3])+1)+'\t'+ 'PASS'+'\t'+recinfo)


    with open(outputpath, "w") as f:
        f.write(myvcf)


# In[20]:


def batchdata(data,batch_size,step,window = 200):
    data = data[:,:]
    
    #data = data.reshape(-1,8)
    #data[:,-3:]=0
    #print(data)
    #print(data.shape)
    if step != 0:
        data = data.reshape(-1,9)[step:(step - window)]
    data = data.reshape(-1,200,9)

    size = data.shape[0]//batch_size
    size_ = data.shape[0]%batch_size
    #print(data.shape[0],size_)
    return data[:size*batch_size],data[size*batch_size:]
def predcit_step(base,predict):
    for i in range(len(predict)):
        if predict[i] >= 0.5:
            base[i] = 1
            base[i+1] = 1
    return base    


# In[10]:


def batchdata_geno(data,batch_size,step,window = 200):
    data = data[:,:]
    #data = data.reshape(-1,8)
    #data[:,-3:]=0
    #print(data)
    #print(data.shape)
    if step != 0:
        data = data.reshape(-1,9)[step:(step - window)]
    data = data.reshape(-1,200,9,1)

    size = data.shape[0]//batch_size
    size_ = data.shape[0]%batch_size
    #print(data.shape[0],size_)
    return data[:size*batch_size],data[size*batch_size:]


# In[11]:


import re
import pysam
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Input, Dense, Conv2D, DepthwiseConv2D, BatchNormalization, Dropout, GlobalAveragePooling2D, Reshape, multiply, add, Activation
#tf.config.set_visible_devices([], 'GPU')
tf.random.set_seed(123)
from os import listdir
from os.path import isfile, join
from tensorflow import keras
import tensorflow_addons as tfa
import math
from tensorflow.keras.utils import Sequence
from tensorflow.keras import backend as K
from collections import Counter
import os
import time


# In[12]:


def predict_funtion(deletion_predict_weight,genotype_predict_weight,datapath,bamfilepath,outvcfpath,contigg,support):
    start = 0
    infor_05 = []
    contig2length = {}
    bamfile = pysam.AlignmentFile(bamfilepath,'rb')
    
    if len(contigg) == 0:
        contig = []
        for count in range(len(bamfile.get_index_statistics())):
            contig.append(bamfile.get_index_statistics()[count].contig)
            contig2length[bamfile.get_index_statistics()[count].contig] = bamfile.lengths[count]
    else:
        contig = np.array(contigg).astype(str)
    for count in range(len(bamfile.get_index_statistics())):
            contig2length[bamfile.get_index_statistics()[count].contig] = bamfile.lengths[count]
    data = []
    resultlist = [['CONTIG', 'START', 'SVLEN', 'READ_SUPPORT', 'SVTYPE']]
    window_size = 200

    predict_del = gru_model()
    predict_del.summary()
    predict_del.load_weights(deletion_predict_weight)
    
    for ww in contig:
        chr_name = ww
        chr_length = contig2length[ww]
        if 'chr' in chr_name:
            chr_name1 = chr_name[-3:]
        else:
            chr_name1 = chr_name
        ider = math.ceil(chr_length / 10000000)
        start = 0
        print('+++++++',chr_name,'++++++++++')
        print('chr:',chr_name,ider)
        for i in range(ider):
            print('deletion_predict_chr:',chr_name,i,'/',ider)
            try:
                x_train_name = datapath +'/chr'+ chr_name1 + '_' + str(start)  +  '_' + str(start + 10000000) + '.npy'
                #print(x_train_name)
                x_t = np.load(x_train_name)

            except FileNotFoundError:
                start = start + 10000000
                continue
            else:
                index_name = datapath +'/chr'+ chr_name1 + '_' + str(start)  +  '_' + str(start + 10000000) + '_index.npy'
                index  = np.load(index_name)

                if len(x_t) == 0:
                    continue
                data1,data2 = batchdata(x_t,100,0)
                data1 = data1.reshape(-1,100,200,9)
                #print(data1.shape,data2.shape)
                datatmp1 = tf.data.Dataset.from_tensor_slices(data1).batch(80)
                predict1 = predict_del.predict(datatmp1)
                if len(data2)!=0:
                    data2 = data2.reshape(1,-1,200,9,1)
                    #datatmp2 = tf.data.Dataset.from_tensor_slices(data2)
                    predict2 = predict_del.predict(data2)
                else:
                    predict2 = np.array([])
                predict1 = predict1.flatten()
                predict2 = predict2.flatten()
                lis = []
                lis.extend(predict1)
                lis.extend(predict2)
                base = np.array(lis)
                data3,data4 = batchdata(x_t,100,100)
                data3 = data3.reshape(-1,100,200,9)
                datatmp3 = tf.data.Dataset.from_tensor_slices(data3).batch(80)
                predict3 = predict_del.predict(datatmp3)
                if len(data4)!=0:
                    data4 = data4.reshape(1,-1,200,9,1)

                    predict4 = predict_del.predict(data4)
                else:
                    predict4 = np.array([])
                predict3 = predict3.flatten()
                predict4 = predict4.flatten()         
                lis2 = []
                lis2.extend(predict3)
                lis2.extend(predict4)
                base1 = np.array(lis2)
                base = predcit_step(base,base1)
                contig, start = chr_name, start
                
                resultlist += mergedeleton_long(base,index,contig,bamfile)
                
                #print(start)
                start = start + 10000000
                #history = model.evaluate(x_t[:,:-1].reshape(-1, 200, 6), x_t[:,-1],batch_size = 2048)
                #print(chr_name,i)
    dell = []
    for read in resultlist[1:]:
        if read[3] >= int(support) and read[2] >= 50:
            
            dell.append([read[0],read[1],(read[1]+read[2]),'DEL','.',read[2]])
            
            
            
    start = 0
    if len(contigg) == 0:
        contig = []
        for count in range(len(bamfile.get_index_statistics())):
            contig.append(bamfile.get_index_statistics()[count].contig)
            contig2length[bamfile.get_index_statistics()[count].contig] = bamfile.lengths[count]
    else:
        contig = np.array(contigg).astype(str)
    for count in range(len(bamfile.get_index_statistics())):
            contig2length[bamfile.get_index_statistics()[count].contig] = bamfile.lengths[count]
    data = []
    window_size = 200
    predict_geno = geno_model()
    predict_geno.summary()
    predict_geno.load_weights(genotype_predict_weight)

    for qq in contig:
        chr_name = qq
        #print(chr_name)
        chr_length = contig2length[qq]
        if 'chr' in chr_name:
            chr_name1 = chr_name[-3:]
        else:
            chr_name1 = chr_name
        ider = math.ceil(chr_length / 10000000)
        start = 0
        print('+++++++',chr_name,'++++++++++')
       # print(chr_name,ider)
        for i in range(ider):
            print('genotype_predict_chr:',chr_name,i,'/',ider)
            try:
                x_train_name = datapath +'/chr'+ chr_name1 + '_' + str(start) +  '_' + str(start + 10000000) + '.npy'
                x_t = np.load(x_train_name)

            except FileNotFoundError:
                start = start + 10000000
                continue
            else:
                index_name = datapath +'/chr'+ chr_name1 + '_' + str(start) +  '_' + str(start + 10000000) + '_index.npy'
                index  = np.load(index_name)
                if len(x_t) == 0:
                    continue
        
                data1,data2 = batchdata_geno(x_t,100,0)   
       
                dataa = tf.data.Dataset.from_tensor_slices(data1).batch(1000)
                predict1 = predict_geno.predict(dataa)
                #predict1 = tf.argmax(predict1, axis = 1).numpy()
                if len(data2)!=0:
                    data2 = data2.reshape(-1,200,9,1)

                    predict2 = predict_geno.predict(data2)
                    #predict2 = tf.argmax(predict2, axis = 1).numpy()
                else:
                    predict2 = np.array([])
                predict1 = predict1.flatten()
                #print(predict1)
                predict2 = predict2.flatten()
                lis = []
                lis.extend(predict1)
                lis.extend(predict2)
                base = np.array(lis)
                contig, start = chr_name, start
                start = start + 10000000
                for kkk in range(len(base)):
                    if base[kkk] > 0.5:
                        data.append([chr_name,index[kkk],index[kkk]+200,1])
                    else:
                        data.append([chr_name,index[kkk],index[kkk]+200,0])       
    final_vcf = find_geno(dell,np.array(data))
    tovcf(final_vcf,contig2length,outvcfpath)
                


   


# In[ ]:




