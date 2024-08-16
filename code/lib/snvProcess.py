#----------------------------------------------------------------------------------------------------------------
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
import sys
from lib.tools import *


def snv_process(self, row):
    sequence_name=row['transcript']
    tubian_leixing=row['AA']
    tubian_leixing_base=row['base']
    if '-' in tubian_leixing_base:
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None})
    if tubian_leixing == 'p.M1?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None})    
    if tubian_leixing == 'p.?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None})
    # pdb.set_trace()
    ####################
    weizhi,yuanshi,zhihou=extract_mutation_info(tubian_leixing)
    protein_sequence=get_gene_info_by_NM(self.NM_pkl_df,sequence_name).translate()
    mutable_seq = MutableSeq(protein_sequence)
    # 进行突变
    mutable_seq[weizhi - 1] = zhihou
    # 转换为不可变的氨基酸序列
    mutated_protein_sequence = Seq(str(mutable_seq))
    # if weizhi<=self.output_fa_len:
    #     front_AA=str(mutated_protein_sequence[0:weizhi])
    # else:
    #     front_AA=str(mutated_protein_sequence[weizhi-self.output_fa_len:weizhi])
    
    # if zhihou=='*':
    #     back_AA='*'
    # else:
    #     back_AA=str(mutated_protein_sequence[weizhi-1:weizhi+self.output_fa_len-1])
        
    WT_AA=str(protein_sequence).split('*')[0]+'*'
    if '*' in str(mutated_protein_sequence):
        protein_seq=str(mutated_protein_sequence)
        protein_seq=protein_seq.split('*')[0]+'*'
        MT_AA=protein_seq
    else : 
        MT_AA=str(protein_sequence)
    #### 2024年08月06日11:15 更新 front_base 和 back_base 的start和end 使用find_difference_range 函数获取
    start,end=find_difference_range(WT_AA,MT_AA)
    if start<=self.output_fa_len:
        front_AA=str(mutated_protein_sequence[0:start+1])
    else: 
        front_AA=str(mutated_protein_sequence[start+1-self.output_fa_len:start+1])
    if zhihou=='*':
        back_AA='*'
    else:
        back_AA=str(mutated_protein_sequence[end-1:end+self.output_fa_len-1])
    #####################################2024年08月06日11:15 更新结束
    weizhi,yuanshi,zhihou=extract_mutation_info_for_acid(tubian_leixing_base)

    base_sequence=get_gene_info_by_NM(self.NM_pkl_df,sequence_name)
    mutable_seq_base = MutableSeq(base_sequence)
    # pdb.set_trace()
    # 进行突变
    mutable_seq_base[weizhi - 1] = zhihou
    # 转换为不可变的氨基酸序列
    mutated_protein_sequence = Seq(str(mutable_seq_base))
    if weizhi<=self.output_fa_len*3:
        front_base_str=str(mutated_protein_sequence[0:weizhi-1])
    else:
        front_base_str=str(mutated_protein_sequence[weizhi-1-self.output_fa_len*3:weizhi-1])
    back_base_str=str(mutated_protein_sequence[weizhi:weizhi+self.output_fa_len*3])
    str_base_WT=str(base_sequence)
    str_base_MT=str(mutated_protein_sequence)
    return pd.Series({"front_base": front_base_str, "back_base": back_base_str,'front_AA':front_AA,'back_AA':back_AA,'WT_base':str_base_WT,'MT_base':str_base_MT,'WT_AA':WT_AA,'MT_AA':MT_AA})

