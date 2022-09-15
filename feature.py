#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pysam 
import time
import numpy as np
from collections import Counter,defaultdict
import pandas as pd
import math
import multiprocessing
from multiprocessing import Process,Queue,Lock
from multiprocessing import Process
from multiprocessing.sharedctypes import Value, Array
import os


# In[2]:


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
def splitread(chr_name,bamfile):
    dada = []

    for read in bamfile.fetch(chr_name):

        #print(read.has_tag('SA'))
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
                        #print(tmpcontig)
                        refstart_1, refend_1, readstart_1, readend_1 =  read.reference_start, read.reference_end,read.query_alignment_start,read.query_alignment_end

                        refstart_2, refend_2, readstart_2, readend_2 = c_pos(cigar, tmprefstart)
                        #print(refstart_2, refend_2, readstart_2, readend_2)
                        a = readend_1 - readstart_2
                        b = refend_1 - refstart_2
                        #print(b-a)
                        if(abs(a)<2000):
                            if(abs(b-a)<30):
                                continue
                            if((b-a)<0 and abs(b-a)<200000):
                                #print(max(refend_1, refstart_2),end)
                                #if(max(refend_1, refstart_2)<end and max(refend_1, refstart_2)<start):

                                data1 = [min(refend_1, refstart_2), min(refend_1, refstart_2)+abs(b-a),abs(b-a)]

                                data2 = np.arange(data1[0],data1[1])
                                dada.extend(data2)
    data = pd.value_counts(dada)
    return data


# In[3]:


def feature_extraction_long(bamfile,chro,start,end,mapq):
    ref_pos = []
    del_count = []
    loci_clip_sm = []
    loci_clip_ms = []
    for read in bamfile.fetch(chro,start,end):
        aligned_length = read.reference_length
        if aligned_length == None:
            aligned_length= 0
        if (read.mapping_quality >= 0) and (aligned_length >= 0) :
            cigar = np.array(read.cigartuples)
            
            ref_pos +=(read.get_reference_positions())
            ref_pos_start = read.reference_start + 1 
            cigar_shape = cigar.shape
            for i in range(cigar_shape[0]):
                if cigar[i,0] == 0:  
                    ref_pos_start = cigar[i,1] + ref_pos_start  
                elif cigar[i,0] == 7:  
                    ref_pos_start = cigar[i,1] + ref_pos_start  
                elif cigar[i,0] == 8:
                    ref_pos_start = cigar[i,1] + ref_pos_start  
                elif cigar[i,0] == 2 :
                    for j  in range(cigar[i,1]):
                        del_count += [ref_pos_start]
                        ref_pos_start = ref_pos_start + 1

            if cigar[0,0] == 4 :
                loci_clip_sm.append(read.reference_start+1)
            if cigar[-1,0] == 4 :
                loci_clip_ms.append(read.reference_end)
 
    ref_pos1 = np.array(ref_pos) + 1
    data = pd.value_counts(ref_pos1)
    del_countt = pd.value_counts(del_count)
    loci_clip_sm = np.array(loci_clip_sm)
    loci_clip_ms = np.array(loci_clip_ms)
    
    loci_clip_sm = pd.value_counts(loci_clip_sm)
    loci_clip_ms = pd.value_counts(loci_clip_ms)
    return data,del_countt,loci_clip_sm,loci_clip_ms


# In[4]:


def feature_extraction_short(bamfile,chro,start,end,mapq):
    ref_pos = []
    del_count = []
    loci_clip_sm=[]
    loci_clip_ms=[]
    insert_size = []
    for read in bamfile.fetch(chro,start,end):
        aligned_length = read.reference_length
        if aligned_length == None:
            aligned_length= 0 
        if (read.mapping_quality >= mapq) and aligned_length > 0:
            cigar = np.array(read.cigartuples)

            ref_pos +=(read.get_reference_positions())

            cigar = np.array(read.cigartuples)
            cigar_shape = cigar.shape

            if cigar[0,0] == 4 :
                loci_clip_sm.append(read.reference_start+1)
            if cigar[-1,0] == 4 :
                loci_clip_ms.append(read.reference_end)
                    

 
    ref_pos1 = np.array(ref_pos) + 1
    data = pd.value_counts(ref_pos1)
    loci_clip_sm = np.array(loci_clip_sm)
    loci_clip_ms = np.array(loci_clip_ms)
    
    loci_clip_sm = pd.value_counts(loci_clip_sm)
    loci_clip_ms = pd.value_counts(loci_clip_ms)
    return data,loci_clip_sm,loci_clip_ms


# In[5]:


def compute(converage_long,loci_del_long,loci_clip_long_sm,loci_clip_long_ms,split_read,converage_short,loci_clip_short_sm,loci_clip_short_ms,start,end):
    s_e = np.arange(start,end)
    converage_long = converage_long.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_del_long = loci_del_long.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_clip_long_sm = loci_clip_long_sm.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_clip_long_ms= loci_clip_long_ms.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_del_split = split_read.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    data_mm = np.ones([len(s_e),1]) * np.mean(converage_long)
    converage_short = converage_short.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_clip_short_sm = loci_clip_short_sm.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_clip_short_ms = loci_clip_short_ms.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    #print(converage_long.shape,loci_del_long.shape,loci_clip_long_sm.shape,loci_clip_long_ms.shape,loci_del_split.shape,data_mm.shape,converage_short.shape,loci_clip_short_sm.shape,loci_clip_short_ms.shape)
    infor = np.concatenate([converage_long,loci_del_long,loci_clip_long_sm,loci_clip_long_ms,loci_del_split,data_mm,converage_short,loci_clip_short_sm,loci_clip_short_ms],axis = 1)
    return infor


# In[6]:


def fun(data):
    oshape = data.shape
    data = data.reshape(-1,9).astype('float32')
    data_long = data[:,:6]
    data_short = data[:,6:]
    #print(data_long.shape,data_short.shape)
    #a = a.astype('float32')
    data_long_mean = np.mean(data_long,axis = 1).reshape(-1,1)
    data_long = data_long - data_long_mean
    sqrt_long = (np.sqrt(data_long.var(axis =1))+1e-10).reshape(-1,1)
    #print(a)
    data_long = data_long/sqrt_long
    
    data_short_mean = np.mean(data_short,axis = 1).reshape(-1,1)
    data_short = data_short - data_short_mean
    sqrt_short = (np.sqrt(data_short.var(axis =1))+1e-10).reshape(-1,1)
    #print(a)
    data_short = data_short/sqrt_short
    datao = np.concatenate([data_long,data_short],axis = 1)
    #print(datao.shape)
    return datao.reshape(oshape)


# In[7]:


def labeldata(vcfpath, contig, start, end, window_size, index):
  goldl = []
  if('chr' in contig):
    contig = contig[3:]
  for rec in pysam.VariantFile(vcfpath).fetch():

    if(rec.contig != contig):
      continue            
    if((rec.info['SVTYPE'] == 'DEL')):
      goldl.append([rec.start, rec.stop, rec.stop - rec.start, 1])
        
  
    
  goldl = (pd.DataFrame(goldl).sort_values([0, 1]).values).astype('float64')


  y = []
  for rec in index:
        
    if(((goldl[:,1:2] > rec) & (goldl[:,:1] < (rec+window_size))).sum() != 0):
      y.append((((goldl[:,1:2] > rec) & (goldl[:,:1] < (rec+window_size))) * goldl[:,3:]).sum())


    else:
      y.append(0)
  return (np.array(y)>0).astype('float32')


# In[8]:


def create_data_longshort(bamfile_long_path,bamfile_short_path,outputpath,contig):
    time_st = time.time()
    bamfile_short = pysam.AlignmentFile(bamfile_short_path,'rb', threads = 20)
    ref_name_short = bamfile_short.get_reference_name
    chr_length_short = bamfile_short.lengths
    bamfile_long = pysam.AlignmentFile(bamfile_long_path,'rb', threads = 20)
    ref_name_long = bamfile_long.get_reference_name
    chr_length_long = bamfile_long.lengths
    contig2length = {}
    window = 200
    if len(contig) == 0:
        contig = []
        for count in range(len(bamfile_long.get_index_statistics())):
            contig.append(bamfile_long.get_index_statistics()[count].contig)
            contig2length[bamfile_long.get_index_statistics()[count].contig] = bamfile_long.lengths[count]
    else:
        contig = np.array(contig).astype(str)
    for count in range(len(bamfile_long.get_index_statistics())):
            contig2length[bamfile_long.get_index_statistics()[count].contig] = bamfile_long.lengths[count]
    for ww in contig:
        chr_name_long = ww
        chr_name_short = ww
        chr_length = contig2length[ww]
        ider = math.ceil(chr_length/10000000)
        #max_length = ider * 10000000
        start = 0
        end = 10000000
        s = 0
        print(chr_name_long,ider)
        time_q = time.time()
        split_read = splitread(chr_name_long,bamfile_long)
        for n in range (ider):
            time_s = time.time()
            read_mes = []
            x_data = []
            index = []
            print('chr',chr_name_long,'start:',start,'end:',end,n,'/',ider)


            loci_cover_long ,loci_del_long,loci_clip_long_sm,loci_clip_long_ms = feature_extraction_long(bamfile_long,chr_name_long,start,end,20)
            loci_cover_short,loci_clip_short_sm,loci_clip_short_ms = feature_extraction_short(bamfile_short,chr_name_short,start,end,20)

            #print(loci_cover_long.shape ,loci_del_long.shape,loci_clip_long.shape,loci_cover_short.shape,loci_clip_short.shape,insert_size_short.shape)
            xx= compute(loci_cover_long ,loci_del_long,loci_clip_long_sm,loci_clip_long_ms,split_read,loci_cover_short,loci_clip_short_sm,loci_clip_short_ms,start,end)
            xx_test = xx[:,:5].reshape(-1,1000)
            xx = xx.reshape(-1,1800)
            #print(time.time()-time_s)
            for k in range(len(xx_test)):
            
                if xx_test[k].any() != 0:
                    x_data.append(xx[k])
                    index.append(s)
                s += window
            x_data = np.array(x_data)
            index = np.array(index)
            #print(x_data.shape,index.shape,chr_name_long)
            #y_label = labeldata('/tf/home/gaoruntian/code/structual_variant_practice/result/HG002_SVs_Tier1_v0.6.vcf.gz',chr_name_long,start,end,200,index)
            start = start + 10000000
            end = end + 10000000 
            if len(x_data) == 0:
                continue
            x_data = fun(x_data) 
            if 'chr' in chr_name_long:
                filename_data = outputpath + '/'   + chr_name_long + '_' + str(start-10000000)+ '_' + str(end-10000000) +'.npy'
                filename_index = outputpath + '/'  + chr_name_long + '_' + str(start-10000000)+ '_' + str(end-10000000) +'_index.npy'
            else:
                filename_data = outputpath + '/chr'   + chr_name_long + '_' + str(start-10000000)+ '_' + str(end-10000000) +'.npy'
                filename_index = outputpath + '/chr'  + chr_name_long + '_' + str(start-10000000)+ '_' + str(end-10000000) +'_index.npy'

            np.save(filename_data,x_data)
            np.save(filename_index,index)
 
            time_e = time.time()
            print(time_e - time_s)
        print(time.time()-time_q) 
    print(time.time()-time_st)


# In[9]:


def compute_long(converage_long,loci_del_long,loci_clip_long_sm,loci_clip_long_ms,split_read,start,end):
    s_e = np.arange(start,end)
    converage_long = converage_long.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_del_long = loci_del_long.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_clip_long_sm = loci_clip_long_sm.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_clip_long_ms= loci_clip_long_ms.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    loci_del_split = split_read.reindex(index = s_e ).fillna(value=0).values.reshape(-1,1)
    data_mm = np.ones([len(s_e),1]) * np.mean(converage_long)
    short = np.zeros([len(s_e),1])
    
    infor = np.concatenate([converage_long,loci_del_long,loci_clip_long_sm,loci_clip_long_ms,loci_del_split,data_mm,short,short,short],axis = 1)
    return infor


# In[10]:


def create_data_long(bamfile_long_path,outputpath,contig):
    time_st = time.time()
    bamfile_long = pysam.AlignmentFile(bamfile_long_path,'rb', threads = 20)
    ref_name_long = bamfile_long.get_reference_name
    chr_length_long = bamfile_long.lengths
    contig2length = {}
    window = 200
    if len(contig) == 0:
        contig = []
        for count in range(len(bamfile_long.get_index_statistics())):
            contig.append(bamfile_long.get_index_statistics()[count].contig)
            contig2length[bamfile_long.get_index_statistics()[count].contig] = bamfile_long.lengths[count]
    else:
        contig = np.array(contig).astype(str)
    for count in range(len(bamfile_long.get_index_statistics())):
        contig2length[bamfile_long.get_index_statistics()[count].contig] = bamfile_long.lengths[count]
    for ww in contig:
        chr_name_long = ww
        chr_name_short = ww
        chr_length = contig2length[ww]
        ider = math.ceil(chr_length/10000000)
        #max_length = ider * 10000000
        start = 0
        end = 10000000
        s = 0
        print(chr_name_long,ider)
        time_q = time.time()
        split_read = splitread(chr_name_long,bamfile_long)
        for n in range (ider):
            time_s = time.time()
            read_mes = []
            x_data = []
            index = []
            print('chr',chr_name_long,'start:',start,'end:',end,n,'/',ider)


            loci_cover_long ,loci_del_long,loci_clip_long_sm,loci_clip_long_ms = feature_extraction_long(bamfile_long,chr_name_long,start,end,20)
            #loci_cover_short,loci_clip_short_sm,loci_clip_short_ms = feature_extraction_short(bamfile_short,chr_name_short,start,end,20)

            #print(loci_cover_long.shape ,loci_del_long.shape,loci_clip_long.shape,loci_cover_short.shape,loci_clip_short.shape,insert_size_short.shape)
            xx= compute_long(loci_cover_long ,loci_del_long,loci_clip_long_sm,loci_clip_long_ms,split_read,start,end)
            xx_test = xx[:,:5].reshape(-1,1000)
            xx = xx.reshape(-1,1800)
            #print(xx.shape)

            for k in range(len(xx_test)):
                if xx_test[k].any() != 0:
                    #print(xx_test[k])
                    x_data.append(xx[k])
                    index.append(s)
                s += window
            #print(x_data)
            x_data = np.array(x_data)
            index = np.array(index)
            
            #print(x_data.shape,index.shape)
            start = start + 10000000
            end = end + 10000000 
            if len(x_data) == 0:
                continue
            x_data = fun(x_data)
            #y_label = y_label.reshape(-1,1)
            #print(x_data.shape,y_label.shape)

            #x_data = np.concatenate([x_data,y_label],axis = 1)
            #print(x_data.shape)
            if 'chr' in chr_name_long:
                filename_data = outputpath + '/'   + chr_name_long + '_' + str(start-10000000)+ '_' + str(end-10000000) +'.npy'
                filename_index = outputpath + '/'  + chr_name_long + '_' + str(start-10000000)+ '_' + str(end-10000000) +'_index.npy'
            else:
                filename_data = outputpath + '/chr'   + chr_name_long + '_' + str(start-10000000)+ '_' + str(end-10000000) +'.npy'
                filename_index = outputpath + '/chr'  + chr_name_long + '_' + str(start-10000000)+ '_' + str(end-10000000) +'_index.npy'
            np.save(filename_data,x_data)
            np.save(filename_index,index)

            time_e = time.time()
            print(time_e - time_s)
        print(time.time()-time_q) 
    print(time.time()-time_st)


# In[15]:


def create_data_long_mul(bamfile_long_path,outputpath,contig,max_work = 5):

    bamfile_long = pysam.AlignmentFile(bamfile_long_path,'rb', threads = 20)
    contig2length = {}
    if len(contig) == 0:
        contig = []
        for count in range(len(bamfile_long.get_index_statistics())):
            contig.append(bamfile_long.get_index_statistics()[count].contig)
            contig2length[bamfile_long.get_index_statistics()[count].contig] = bamfile_long.lengths[count]
    else:
        contig = np.array(contig).astype(str)
    count = 0
    while((count) < len(contig)):
            if(len(multiprocessing.active_children()) < int(max_work)): 
                    j = contig[count]

                    p = Process(target=create_data_long, args=(bamfile_long_path,outputpath,[j])) #实例化进程对象
                    p.start()
                    count += 1
            else:
                  time.sleep(2)    


# In[16]:


def create_data_longshort_mul(bamfile_long_path,bamfile_short_path,outputpath,contig,max_work = 5):
    bamfile_long = pysam.AlignmentFile(bamfile_long_path,'rb', threads = 20)
    contig2length = {}
    if len(contig) == 0:
        contig = []
        for count in range(len(bamfile_long.get_index_statistics())):
            contig.append(bamfile_long.get_index_statistics()[count].contig)
            contig2length[bamfile_long.get_index_statistics()[count].contig] = bamfile_long.lengths[count]
    else:
        contig = np.array(contig).astype(str)
    count = 0
    while((count) < len(contig)):
            if(len(multiprocessing.active_children()) < int(max_work)): 
                    j = contig[count]

                    p = Process(target=create_data_longshort, args=(bamfile_long_path,bamfile_short_path,outputpath,[j])) #实例化进程对象
                    p.start()
                    count += 1
            else:
                  time.sleep(2) 







# In[ ]:




