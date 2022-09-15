# LSnet
An approach for detecting and genotyping deletions with a deeplearning network.
## Installation
### Requirements
* python 3.9, numpy, pandas, Matplotlib, TensorFlow 2.7, pysam,tensorflow_addons 0.17.1
### 1. Create a virtual environment  
```
#create
conda create -n LSnet python=3.9
#activate
conda activate LSnet
#deactivate
conda deactivate
```   
### 2. clone LSnet
* After creating and activating the LSnet virtual environment, download LSnet from github:
```　 
git clone https://github.com/eioyuou/LSnet.git
cd LSnet
```
### 3. Install 
```　
conda activate LSnet
conda install numpy, pandas, Matplotlib, TensorFlow 2.7, pysam
pip install tensorflow-addons
```
## Usage
### 1.Produce data for call SV
```　 
python LSnet.py create_feature bamfile_path_long bamfile_path_short output_data_folder max_work(default:5) includecontig   
    
bamfile_path_long is the path of the alignment file about the reference and the long read set;  
bamfile_path_short is the path of the alignment file about the reference and the short read set(If there is no short read data, use None);  
output_data_folder is a folder which is used to store evaluation data;  
max_work is the number of threads to use;  
includecontig is the list of contig to preform detection.(default: [], all contig are used)  
   
eg: python LSnet.py create_feature ./long_read.bam /short_read.bam ./outpath 5 [12,13,14,15,16,17,18,19,20,21,22] 
    python LSnet.py create_feature ./CCS.bam None ./outpath 5 [12,13,14,15,16,17,18,19,20,21,22]
``` 
### 2.Call deletion 
```　 
python LSnet.py call_sv deletion_predict_weight,genotype_predict_weight,datapath,bamfilepath,outvcfpath,support,includecontig   
   
deletion_predict_weight and genotype_predict_weight are the paths of the model weights;  
datapath is a folder which is used to store evaluation data;  
bamfilepath is the path of the alignment file about the reference and the long read set;  
outvcfpath is the path of output vcf filel;  
support is min support reads;  
includecontig is the list of contig to preform detection.(default: [], all contig are used)  
   
eg: python LSnet.py call_sv ./deletion.h5 ./genotype.h5 ./datapath ./long_read.bam ./outvcfpath 10 [12,13,14,15,16,17,18,19,20,21,22]  
```  
## Tested data 
The data can be downloaded from:  

### HG002 long read data
```
https://ftp.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/HG002_PB_70x_RG_HP10XtrioRTG.bam
```
### HG002 short read data
```
https://ftp.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/HG002_PB_70x_RG_HP10XtrioRTG.bam)
```   
### HG002 ccs data
``` 
https://ftp.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/alignment/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
 ```  
### NA19240 long read data
```
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160905_smithm_pacbio_aligns/NA19240_bwamem_GRCh38DH_YRI_20160905_pacbio.bam
```   
### NA19240 short read data
```
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/YRI/NA19240/high_cov_alignment/NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.cram
```   

