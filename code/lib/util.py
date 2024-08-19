

import argparse
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
import sys
import re
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import warnings
import os
from lib.indelProcess import *
from lib.snvProcess import *

warnings.filterwarnings("ignore", message="A worker stopped while some jobs were given to the executor.", category=UserWarning)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
# 忽略特定警告
from Bio import BiopythonWarning
warnings.filterwarnings("ignore", category=BiopythonWarning, message="Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.")

class NeoantigenSNVfinding:
    def __init__(self,input_file,df_pkl,output_file,output_fa_len) -> None:
        ### 初始化输出参数
        self.input_file = input_file
        self.df_pkl=df_pkl
        self.output_file=output_file
        self.output_fa_len=output_fa_len
        # 调用初始化方法
        self.mutation_result = pd.read_csv(self.input_file)
        self.NM_pkl_df=pd.read_pickle(self.df_pkl)
        self.file_dir=os.path.dirname(output_file)
        # 检查数据是否为空
        if  self.mutation_result.empty:
            print("Error: Mutation data is empty. Exiting script.")
            sys.exit(1)
    def call_peptide(self,row):
        if row['type']=='snv':
            return self.snv_peptide(row)
        if row['type']=='indel':
            return self.indel_peptide(row)
    def snv_peptide(self,row):
        if '-' in row['base']:
            return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
        return snv_process(self,row)
    def indel_peptide(self,row):
        if '-' in row['base']:
            return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
        if 'delins' in row['base']:
            return delins_peptide(self,row)
        if 'inv' in row['base']:
            return inv_peptide(self,row)
        if 'del' in row['base']:
            return del_peptide(self,row)
        if 'ins' in row['base']:
            return ins_peptide(self,row)
        if 'dup' in row['base']:
            return dup_peptide(self,row)
##### 2024年08月02日 这个函数更新是为了 找到两个序列不同的那个划窗位置的起始和终止位置
# 定义一个处理函数
def process_rows(row,len_fa):
    WT = row['WT_AA']
    MT = row['MT_AA']
    start, end = find_difference_range(WT, MT)
    seq = get_seq(row,start,end,len_fa)
    return pd.Series({'seq':seq})
    # return pd.Series({'start_diff': start, 'end_diff': end , 'seq':seq})
def get_seq(row,start,end,len_fa):
    if start==None:
        return None
    start=max(start-len_fa+1,0)  ### 因为加入是15个长度，突变位置前面+14个 所以要-1
    end=end+len_fa-1  ### 因为加入是15个长度，突变位置后面+14个 所以要-1
    if 'fs' in row['AA']:
        end =None
    seq=row['MT_AA'][start:end]
    return seq
def find_difference_range(WT, MT):
    if WT==None or MT==None:
        return None,None
    len1, len2 = len(WT), len(MT)
    min_length = min(len1, len2)
    max_length = max(len1, len2)
    if len1 == len2:
        # 找到第一个不同的字符位置
        start = 0
        while start < min_length and WT[start] == MT[start]:
            start += 1
        start_str=WT[start:]
        end_str=MT[start:]
        end=max_length-start
        while end != 0  and start_str[end-1] == end_str[end-1]:
            end -= 1
        end=end+start
            
    if len1 <= len2:
        start = 0
        while start < min_length and WT[start] == MT[start]:
            start += 1
        start_str='@'*(len2-len1)+WT[start:]
        end_str=MT[start:]
        end=max_length-start
        while end != 0  and start_str[end-1] == end_str[end-1]:
            end -= 1
        end=end+start
        
    if len1 >= len2:
        start = 0
        while start < min_length and WT[start] == MT[start]:
            start += 1
        start_str=WT[start:]
        end_str='@'*(len1-len2)+MT[start:]
        end=max_length-start
        while end != 0  and start_str[end-1] == end_str[end-1]:
            end -= 1
        end=end+start
        end=end-(len1-len2)
    return start,end
##### 从seq序列切割#### 下面的代码
import argparse
import pandas as pd
import sys
import re
import warnings
import os
import math
# 如果提供了 --warning_file 参数，则捕获警告信息但不输出到 stderr

#----------------------------------------------------------------------------------------------------------------
def slice_seq(row,len_peptide,file_dir):
#### 2024年07月11日为了stop gain 单独增加
    if row['type']=='snv':
        if row['protein_variant_type_annovar']=='stopgain':
            seq = row['seq']
            num = 1
            if not seq:
                return
            seq = seq.split('*')[0]+'X'
        else:
            seq = row['seq']
            num = 1
            if not seq:
                return
            seq = seq.split('*')[0]
        start = 0
        with open(os.path.join(file_dir,'MT.SNV.{len}.fa'.format(len=len_peptide)), 'a') as file:
            if len(seq) < len_peptide:
                return
            else:
                while num <=len(seq):
                    end = start +len_peptide
                    window = seq[start:end]
                    # pdb.set_trace()
                    if len(window) < len_peptide:
                        break
                    file.write('>'+row['transcript'] + '-' + row['AA'] + '-' + str(len_peptide) + '-' + str(num) + '\n')
                    file.write(window + '\n')
                    start += 1
                    num += 1
                return
    if row['type']=='indel':
            #### 2024年07月31日日为了stop gain 单独增加
        if row['protein_variant_type_annovar']=='stopgain':
            seq = row['seq']
            num = 1
            if not seq:
                return
            seq = seq.split('*')[0]+'X'
        else:
            seq = row['seq']
            num = 1
            if not seq:
                return
            seq = seq.split('*')[0]
        start = 0
        with open(os.path.join(file_dir,'MT.INDEL.{len}.fa'.format(len=len_peptide)), 'a') as file:
            if len(seq) < len_peptide:
                return
            else:
                while num <=len(seq):
                    end = start +len_peptide
                    window = seq[start:end]
                    # pdb.set_trace()
                    if len(window) < len_peptide:
                        break
                    file.write('>'+row['transcript'] + '-' + row['base'] + '-' + str(len_peptide) + '-' + str(num) + '\n')
                    file.write(window + '\n')
                    start += 1
                    num += 1
                return
            
# 删除指定文件（如果存在）
def remove_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)
    # 创建一个新的空文件
    with open(file_path, 'w') as file:
        pass  # 使用 pass 来创建一个空文件
