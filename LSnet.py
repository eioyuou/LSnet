#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import tensorflow as tf
import sys
from model import predict_funtion
from model import gru_model
from feature import create_data_long_mul
from feature import create_data_longshort_mul
import tensorflow as tf
#tf.config.set_visible_devices([], 'GPU')
#model = gru_model()
mode = sys.argv[1]
debug = 0

if(mode == 'create_feature'):
    if(len(sys.argv) not in [4, 5, 6, 7]):
        debug = 1
    else:
        print('Produce data')
        if(len(sys.argv) == 7):
            bamfilepath_long,bamfilepath_short,outputpath, max_work,includecontig  = sys.argv[2], sys.argv[3],sys.argv[4],sys.argv[5], [str(contig) for contig in eval(sys.argv[6])]
        if(len(sys.argv) == 6):
            bamfilepath_long,bamfilepath_short,outputpath, max_work,includecontig = sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5]
        if(len(sys.argv) == 5):
            bamfilepath_long,bamfilepath_short,outputpath, max_work,includecontig = sys.argv[2], sys.argv[3], sys.argv[4],5,[]
        print('bamfile path ', bamfilepath_long,bamfilepath_short)
        print('output path ', outputpath)
        #print('max_worker set to ', str(max_worker))
        if(includecontig == []):
            print('All chromosomes within bamfile will be used')
        else:
            print('Following chromosomes will be used')
            print(includecontig)
        if bamfilepath_short == 'None':
            create_data_long_mul(bamfile_long_path = bamfilepath_long,outputpath=outputpath,contig=includecontig,max_work = max_work)
        else:   
            create_data_longshort_mul(bamfile_long_path = bamfilepath_long, bamfile_short_path = bamfilepath_short,outputpath=outputpath, contig=includecontig,max_work = max_work)
        print('\n\n')
        print('Completed')
        print('\n\n')
        
elif(mode == 'call_sv'):
    if(len(sys.argv) not in  [8, 9]):
        debug = 1
    else:
        print('testing')
        if(len(sys.argv) == 9):
            deletion_predict_weight,genotype_predict_weight,datapath,bamfilepath,outvcfpath,support,contigg = sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7], [str(contig) for contig in eval(sys.argv[8])]
        else:
            deletion_predict_weight,genotype_predict_weight,datapath,bamfilepath,outvcfpath,support,contigg = sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7], []
        
        print('bamfile path ', bamfilepath)
        print('weight path ', deletion_predict_weight,genotype_predict_weight)
        print('data file path ', datapath)
        print('vcf path',outvcfpath)
        if(contigg == []):
            print('All chromosomes within bamfile will be used')
        else:
            print('Following chromosomes will be used')
            print(contigg)
            

        predict_funtion(deletion_predict_weight,genotype_predict_weight,datapath,bamfilepath,outvcfpath,contigg,support)
        #predict_fn(datapath = datapath, weightpath = weightpath, bamfilepath = bamfilepath, includecontig=includecontig )
        print('\n\n')
        print('Completed, Result saved in current folder')
        print('\n\n')

else:
    debug = 1
if(debug ==1):
    print('\n\n')
    print('Useage')
    print('Produce data for call sv')
    print('python LSnet.py create_feature bamfile_path_long bamfile_path_short output_data_folder max_work includecontig(default:[](all chromosomes))')
    print('Call sv')
    print('python LSnet.py call_sv deletion_predict_weight,genotype_predict_weight,datapath,bamfilepath,outvcfpath,support, includecontig(default:[](all chromosomes)')
        

# In[ ]:




