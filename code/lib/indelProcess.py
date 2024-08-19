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

def delins_peptide(self, row):
    found_sequence = None
    sequence_name=row['transcript']
    nucleotide_sequence=get_gene_info_by_NM(self.NM_pkl_df,sequence_name)
    yuanshi_nucleotide_sequence=nucleotide_sequence
    nucleotide_sequence=MutableSeq(nucleotide_sequence)
    tubian_leixing=row['AA']
    if tubian_leixing == 'p.M1?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    if tubian_leixing == 'p.?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    weizhi=extract_mutation_info_for_indel(tubian_leixing)-1
    if 'delins' in row['base']:
        matches = re.findall(r'\d+', row['base'])
        if len(matches) == 2:
            position1 = int(matches[0])-1
            position2 = int(matches[1])
            del nucleotide_sequence[position1:position2]
            inserted_bases=row['base'].split('ins')[-1]
            nucleotide_sequence[position1:position1] = inserted_bases
            nucleotide_sequence=Seq(str(nucleotide_sequence))
            protein_sequence = nucleotide_sequence.translate()
            ###########
            if position1<=self.output_fa_len*3:
                front_base_str=str(nucleotide_sequence[0:position1])
            else:
                front_base_str=str(str(nucleotide_sequence[position1-self.output_fa_len*3:position1]))
            back_base_str=str(nucleotide_sequence[position1+len(inserted_bases):position1+len(inserted_bases)+self.output_fa_len*3])
            str_base_WT=str(yuanshi_nucleotide_sequence)
            str_base_MT=str(nucleotide_sequence)
            ###########
            # if weizhi<=self.output_fa_len:
            #     if '*' in protein_sequence[0:weizhi+1]:
            #         front_AA=str(protein_sequence[0:weizhi+1])+'*'
            #     else:
            #         front_AA=str(protein_sequence[0:weizhi+1])
            # elif '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
            #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]).split('*')[0]+'*')
            # else:
            #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]))
            # if '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
            #     back_AA=''
            # elif '*' in protein_sequence[weizhi:weizhi+self.output_fa_len]:
            #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len].split('*')[0]+'*')
            # else:
            #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len])
            WT_AA=str(yuanshi_nucleotide_sequence.translate()).split('*')[0]+'*'
            if '*' in str(protein_sequence):
                protein_seq=str(protein_sequence)
                protein_seq=protein_seq.split('*')[0]+'*'
                MT_AA=protein_seq
            else : 
                MT_AA=str(protein_sequence)
            #### 2024年08月06日11:15 更新 front_base 和 back_base 的start和end 使用find_difference_range 函数获取
            start,end=find_difference_range(WT_AA,MT_AA)
            if start<=self.output_fa_len:
                front_AA=str(MT_AA[0:start+1])
            else: 
                front_AA=str(MT_AA[start+1-self.output_fa_len:start+1])
            ## back_AA=str(MT_AA[end-1:end+self.output_fa_len-1]) 2024年08月13日 对end的输出进行更改
            if 'fs' in row['AA']:
                end =len(MT_AA)
            else:
                end = start 
            back_AA=str(MT_AA[end:end+self.output_fa_len])
            #####################################2024年08月06日11:15 更新结束
        if len(matches) == 1:
            position1 = int(matches[0])-1
            del nucleotide_sequence[position1]
            inserted_bases=row['base'].split('ins')[-1]
            nucleotide_sequence[position1:position1] = inserted_bases
            nucleotide_sequence=Seq(str(nucleotide_sequence))
            protein_sequence = nucleotide_sequence.translate()
            #############
            if position1<=self.output_fa_len*3:
                front_base_str=str(nucleotide_sequence[0:position1])
            else:
                front_base_str=str(str(nucleotide_sequence[position1-self.output_fa_len*3:position1]))
            back_base_str=str(nucleotide_sequence[position1+len(inserted_bases):position1+len(inserted_bases)+self.output_fa_len*3])
            str_base_WT=str(yuanshi_nucleotide_sequence)
            str_base_MT=str(nucleotide_sequence)
            #############
            # if weizhi<=self.output_fa_len:
            #     if '*' in protein_sequence[0:weizhi+1]:
            #         front_AA=str(protein_sequence[0:weizhi+1])+'*'
            #     else:
            #         front_AA=str(protein_sequence[0:weizhi+1])
            # elif '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
            #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]).split('*')[0]+'*')
            # else:
            #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]))
            # if '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
            #     back_AA=''
            # elif '*' in protein_sequence[weizhi:weizhi+self.output_fa_len]:
            #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len].split('*')[0]+'*')
            # else:
            #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len])
            WT_AA=str(yuanshi_nucleotide_sequence.translate()).split('*')[0]+'*'
            if '*' in str(protein_sequence):
                protein_seq=str(protein_sequence)
                protein_seq=protein_seq.split('*')[0]+'*'
                MT_AA=protein_seq
            else : 
                MT_AA=str(protein_sequence)
            #### 2024年08月06日11:15 更新 front_base 和 back_base 的start和end 使用find_difference_range 函数获取
            start,end=find_difference_range(WT_AA,MT_AA)
            if start<=self.output_fa_len:
                front_AA=str(MT_AA[0:start+1])
            else: 
                front_AA=str(MT_AA[start+1-self.output_fa_len:start+1])
            ## back_AA=str(MT_AA[end-1:end+self.output_fa_len-1]) 2024年08月13日 对end的输出进行更改
            if 'fs' in row['AA']:
                end =len(MT_AA)
                
            else:
                end = start 
            back_AA=str(MT_AA[end:end+self.output_fa_len])
            #####################################2024年08月06日11:15 更新结束
            
    return pd.Series({"front_base": front_base_str, "back_base": back_base_str,'front_AA':front_AA,'back_AA':back_AA,'WT_base':str_base_WT,'MT_base':str_base_MT,'WT_AA':WT_AA,'MT_AA':MT_AA})

def inv_peptide(self, row):
    found_sequence = None
    sequence_name=row['transcript']
    nucleotide_sequence=get_gene_info_by_NM(self.NM_pkl_df,sequence_name)
    yuanshi_nucleotide_sequence=nucleotide_sequence
    nucleotide_sequence=MutableSeq(nucleotide_sequence)
    matches = re.findall(r'\d+', row['base'])
    tubian_leixing=row['AA']
    if tubian_leixing == 'p.M1?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    if tubian_leixing == 'p.?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    weizhi=extract_mutation_info_for_indel(tubian_leixing)-1

    if len(matches) == 2:
        position1 = int(matches[0])-1
        position2 = int(matches[1])
        # 提取要颠倒的区域
        region_to_reverse = nucleotide_sequence[position1:position2]
        # 颠倒区域
        # 定义互补碱基的映射字典
        complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
        # 颠倒区域并取互补配对的核酸
        reversed_region = "".join([complement_dict[base] for base in region_to_reverse[::-1]])
        # 构建新的 DNA 序列，将颠倒后的区域替换回原始序列
        nucleotide_sequence = nucleotide_sequence[:position1] + reversed_region + nucleotide_sequence[position2:]
        nucleotide_sequence=Seq(str(nucleotide_sequence))
        protein_sequence = nucleotide_sequence.translate()
        #############
        if position1<=self.output_fa_len*3:
            front_base_str=str(nucleotide_sequence[0:position1])
        else:
            front_base_str=str(str(nucleotide_sequence[position1-self.output_fa_len*3:position1]))
        back_base_str=str(nucleotide_sequence[position2:position2+self.output_fa_len*3])
        str_base_WT=str(yuanshi_nucleotide_sequence)
        str_base_MT=str(nucleotide_sequence)
        #############
        # if weizhi<=self.output_fa_len:
        #     if '*' in protein_sequence[0:weizhi+1]:
        #         front_AA=str(protein_sequence[0:weizhi+1])+'*'
        #     else:
        #         front_AA=str(protein_sequence[0:weizhi+1])
        # elif '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]).split('*')[0]+'*')
        # else:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]))
        # if '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     back_AA=''
        # elif '*' in protein_sequence[weizhi:weizhi+self.output_fa_len]:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len].split('*')[0]+'*')
        # else:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len])
        WT_AA=str(yuanshi_nucleotide_sequence.translate()).split('*')[0]+'*'
        if '*' in str(protein_sequence):
            protein_seq=str(protein_sequence)
            protein_seq=protein_seq.split('*')[0]+'*'
            MT_AA=protein_seq
        else : 
            MT_AA=str(protein_sequence)
            
        #### 2024年08月06日11:15 更新 front_base 和 back_base 的start和end 使用find_difference_range 函数获取
        start,end=find_difference_range(WT_AA,MT_AA)
        if start<=self.output_fa_len:
            front_AA=str(MT_AA[0:start+1])
        else: 
            front_AA=str(MT_AA[start+1-self.output_fa_len:start+1])
        ## back_AA=str(MT_AA[end-1:end+self.output_fa_len-1]) 2024年08月13日 对end的输出进行更改
        if 'fs' in row['AA']:
            ## end =len(MT_AA)
            end = start 
        else:
            end = start 
        back_AA=str(MT_AA[end:end+self.output_fa_len])
        #####################################2024年08月06日11:15 更新结束
    return pd.Series({"front_base": front_base_str, "back_base": back_base_str,'front_AA':front_AA,'back_AA':back_AA,'WT_base':str_base_WT,'MT_base':str_base_MT,'WT_AA':WT_AA,'MT_AA':MT_AA})

def del_peptide(self, row):
    found_sequence = None
    sequence_name=row['transcript']
    nucleotide_sequence=get_gene_info_by_NM(self.NM_pkl_df,sequence_name)
    yuanshi_nucleotide_sequence=nucleotide_sequence
    nucleotide_sequence=MutableSeq(nucleotide_sequence)
    matches = re.findall(r'\d+', row['base'])
    tubian_leixing=row['AA']
    if tubian_leixing == 'p.M1?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    if tubian_leixing == 'p.?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    weizhi=extract_mutation_info_for_indel(tubian_leixing)-1

    if len(matches) == 2:
        position1 = int(matches[0])-1
        position2 = int(matches[1])
        del nucleotide_sequence[position1:position2]
        nucleotide_sequence=Seq(str(nucleotide_sequence))
        protein_sequence = nucleotide_sequence.translate()
        ##############
        if position1<=self.output_fa_len*3:
            front_base_str=str(nucleotide_sequence[0:position1])
        else:
            front_base_str=str(str(nucleotide_sequence[position1-self.output_fa_len*3:position1]))
        back_base_str=str(nucleotide_sequence[position1:position1+self.output_fa_len*3])
        str_base_WT=str(yuanshi_nucleotide_sequence)
        str_base_MT=str(nucleotide_sequence)
        ##############
        # if weizhi<=self.output_fa_len:
        #     if '*' in protein_sequence[0:weizhi+1]:
        #         front_AA=str(protein_sequence[0:weizhi+1])+'*'
        #     else:
        #         front_AA=str(protein_sequence[0:weizhi+1])
        # elif '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]).split('*')[0]+'*')
        # else:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]))
        # if '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     back_AA=''
        # elif '*' in protein_sequence[weizhi:weizhi+self.output_fa_len]:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len].split('*')[0]+'*')
        # else:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len])
        WT_AA=str(yuanshi_nucleotide_sequence.translate()).split('*')[0]+'*'
        if '*' in str(protein_sequence):
            protein_seq=str(protein_sequence)
            protein_seq=protein_seq.split('*')[0]+'*'
            MT_AA=protein_seq
        else : 
            MT_AA=str(protein_sequence)
        #### 2024年08月06日11:15 更新 front_base 和 back_base 的start和end 使用find_difference_range 函数获取
        start,end=find_difference_range(WT_AA,MT_AA)
        if start<=self.output_fa_len:
            front_AA=str(MT_AA[0:start+1])
        else: 
            front_AA=str(MT_AA[start+1-self.output_fa_len:start+1])
        ## back_AA=str(MT_AA[end-1:end+self.output_fa_len-1]) 2024年08月13日 对end的输出进行更改
        if 'fs' in row['AA']:
            ## end =len(MT_AA)
            end = start 
        else:
            end = start 
        back_AA=str(MT_AA[end:end+self.output_fa_len])
        #####################################2024年08月06日11:15 更新结束
    if len(matches) == 1:
        position = int(matches[0])-1
        prim_seq=found_sequence
        del nucleotide_sequence[position]
        nucleotide_sequence=Seq(str(nucleotide_sequence))
        protein_sequence = nucleotide_sequence.translate()
        if position<=self.output_fa_len*3:
            front_base_str=str(nucleotide_sequence[0:position])
        else:
            front_base_str=str(str(nucleotide_sequence[position-self.output_fa_len*3:position]))
        back_base_str=str(nucleotide_sequence[position:position+self.output_fa_len*3])
        str_base_WT=str(yuanshi_nucleotide_sequence)
        str_base_MT=str(nucleotide_sequence)
        ##################
        # if weizhi<=self.output_fa_len:
        #     if '*' in protein_sequence[0:weizhi+1]:
        #         front_AA=str(protein_sequence[0:weizhi+1])+'*'
        #     else:
        #         front_AA=str(protein_sequence[0:weizhi+1])
        # elif '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]).split('*')[0]+'*')
        # else:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]))
        # if '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     back_AA=''
        #     # back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len].split('*')[0]+'*')
        # elif '*' in protein_sequence[weizhi:weizhi+self.output_fa_len]:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len].split('*')[0]+'*')
        # else:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len])
        WT_AA=str(yuanshi_nucleotide_sequence.translate()).split('*')[0]+'*'
        if '*' in str(protein_sequence):
            protein_seq=str(protein_sequence)
            protein_seq=protein_seq.split('*')[0]+'*'
            MT_AA=protein_seq
        else : 
            MT_AA=str(protein_sequence)
        #### 2024年08月06日11:15 更新 front_base 和 back_base 的start和end 使用find_difference_range 函数获取
        start,end=find_difference_range(WT_AA,MT_AA)
        if start<=self.output_fa_len:
            front_AA=str(MT_AA[0:start+1])
        else: 
            front_AA=str(MT_AA[start+1-self.output_fa_len:start+1])
        ## back_AA=str(MT_AA[end-1:end+self.output_fa_len-1]) 2024年08月13日 对end的输出进行更改
        if 'fs' in row['AA']:
            ## end =len(MT_AA)
            end = start 
        else:
            end = start 
        back_AA=str(MT_AA[end:end+self.output_fa_len])
        #####################################2024年08月06日11:15 更新结束
        ################
    return pd.Series({"front_base": front_base_str, "back_base": back_base_str,'front_AA':front_AA,'back_AA':back_AA,'WT_base':str_base_WT,'MT_base':str_base_MT,'WT_AA':WT_AA,'MT_AA':MT_AA})

def ins_peptide(self, row):
    found_sequence = None
    sequence_name=row['transcript']
    nucleotide_sequence=get_gene_info_by_NM(self.NM_pkl_df,sequence_name)
    yuanshi_nucleotide_sequence=nucleotide_sequence
    nucleotide_sequence=MutableSeq(nucleotide_sequence)
    matches = re.findall(r'\d+', row['base'])
    tubian_leixing=row['AA']
    if tubian_leixing == 'p.M1?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    if tubian_leixing == 'p.?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    weizhi=extract_mutation_info_for_indel(tubian_leixing)-1

    if len(matches) == 2:
        position_1 = int(matches[0])
        position_2 = int(matches[1])
        inserted_bases=row['base'].split('ins')[-1]
        nucleotide_sequence[position_1:position_1] = inserted_bases
        nucleotide_sequence=Seq(str(nucleotide_sequence))
        protein_sequence = nucleotide_sequence.translate()
        ############
        if position_1<=self.output_fa_len*3:
            front_base_str=str(nucleotide_sequence[0:position_1])
        else:
            front_base_str=str(str(nucleotide_sequence[position_1-self.output_fa_len*3:position_1]))
        back_base_str=str(nucleotide_sequence[position_1+len(inserted_bases):position_1+len(inserted_bases)+self.output_fa_len*3])
        str_base_WT=str(yuanshi_nucleotide_sequence)
        str_base_MT=str(nucleotide_sequence)
        ##############
        # if weizhi<=self.output_fa_len:
        #     if '*' in protein_sequence[0:weizhi+1]:
        #         front_AA=str(protein_sequence[0:weizhi+1])+'*'
        #     else:
        #         front_AA=str(protein_sequence[0:weizhi+1])
        # elif '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]).split('*')[0]+'*')
        # else:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]))
        # if '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     back_AA=''
        # elif '*' in protein_sequence[weizhi:weizhi+self.output_fa_len]:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len].split('*')[0]+'*')
        # else:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len])
        WT_AA=str(yuanshi_nucleotide_sequence.translate()).split('*')[0]+'*'
        if '*' in str(protein_sequence):
            protein_seq=str(protein_sequence)
            protein_seq=protein_seq.split('*')[0]+'*'
            MT_AA=protein_seq
        else : 
            MT_AA=str(protein_sequence)
        #### 2024年08月06日11:15 更新 front_base 和 back_base 的start和end 使用find_difference_range 函数获取
        start,end=find_difference_range(WT_AA,MT_AA)
        if start<=self.output_fa_len:
            front_AA=str(MT_AA[0:start+1])
        else: 
            front_AA=str(MT_AA[start+1-self.output_fa_len:start+1])
        ## back_AA=str(MT_AA[end-1:end+self.output_fa_len-1]) 2024年08月13日 对end的输出进行更改
        if 'fs' in row['AA']:
            ## end =len(MT_AA)
            end = start 
        else:
            end = start 
        back_AA=str(MT_AA[end:end+self.output_fa_len])
        #####################################2024年08月06日11:15 更新结束
    return pd.Series({"front_base": front_base_str, "back_base": back_base_str,'front_AA':front_AA,'back_AA':back_AA,'WT_base':str_base_WT,'MT_base':str_base_MT,'WT_AA':WT_AA,'MT_AA':MT_AA})

def dup_peptide(self, row):
    found_sequence = None
    sequence_name=row['transcript']
    nucleotide_sequence=get_gene_info_by_NM(self.NM_pkl_df,sequence_name)
    yuanshi_nucleotide_sequence=nucleotide_sequence
    nucleotide_sequence=MutableSeq(nucleotide_sequence)
    matches = re.findall(r'\d+', row['base'])
    tubian_leixing=row['AA']
    if tubian_leixing == 'p.M1?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    if tubian_leixing == 'p.?':
        return pd.Series({"front_base": None, "back_base": None,'front_AA':None,'back_AA':None,'WT_base':None,'MT_base':None,'WT_AA':None,'MT_AA':None,'beforeStart':None,'Start2End':None,'afterEnd':None})
    weizhi=extract_mutation_info_for_indel(tubian_leixing)-1

    if len(matches) == 2:
        position_1 = int(matches[0])-1
        position_2 = int(matches[1])
        dup_bases=nucleotide_sequence[position_1:position_2]
        nucleotide_sequence[position_1:position_1] = dup_bases
        nucleotide_sequence=Seq(str(nucleotide_sequence))
        protein_sequence = nucleotide_sequence.translate()
        ################
        if position_1<=self.output_fa_len*3:
            front_base_str=str(nucleotide_sequence[0:position_1+1])
        else:
            front_base_str=str(str(nucleotide_sequence[position_1+1-self.output_fa_len*3:position_1+1]))
        back_base_str=str(nucleotide_sequence[position_1+len(dup_bases)+1:position_1+len(dup_bases)+1+self.output_fa_len*3])
        str_base_WT=str(yuanshi_nucleotide_sequence)
        str_base_MT=str(nucleotide_sequence)
        ################
        # if weizhi<=self.output_fa_len:
        #     if '*' in protein_sequence[0:weizhi+1]:
        #         front_AA=str(protein_sequence[0:weizhi+1])+'*'
        #     else:
        #         front_AA=str(protein_sequence[0:weizhi+1])
        # elif '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]).split('*')[0]+'*')
        # else:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]))
        # if '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     back_AA=''
        # elif '*' in protein_sequence[weizhi:weizhi+self.output_fa_len]:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len].split('*')[0]+'*')
        # else:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len])
        WT_AA=str(yuanshi_nucleotide_sequence.translate()).split('*')[0]+'*'
        if '*' in str(protein_sequence):
            protein_seq=str(protein_sequence)
            protein_seq=protein_seq.split('*')[0]+'*'
            MT_AA=protein_seq
        else : 
            MT_AA=str(protein_sequence)
        #### 2024年08月06日11:15 更新 front_base 和 back_base 的start和end 使用find_difference_range 函数获取
        start,end=find_difference_range(WT_AA,MT_AA)
        if start<=self.output_fa_len:
            front_AA=str(MT_AA[0:start+1])
        else: 
            front_AA=str(MT_AA[start+1-self.output_fa_len:start+1])
        ## back_AA=str(MT_AA[end-1:end+self.output_fa_len-1]) 2024年08月13日 对end的输出进行更改
        if 'fs' in row['AA']:
            ## end =len(MT_AA)
            end = start 
        else:
            end = start 
        back_AA=str(MT_AA[end:end+self.output_fa_len])
        #####################################2024年08月06日11:15 更新结束
    if len(matches) == 1:
        position_1 = int(matches[0])-1
        dup_bases=nucleotide_sequence[position_1]
        nucleotide_sequence[position_1:position_1] = dup_bases
        nucleotide_sequence=Seq(str(nucleotide_sequence))
        protein_sequence = nucleotide_sequence.translate()
        if position_1<=self.output_fa_len*3:
            front_base_str=str(nucleotide_sequence[0:position_1+1])
        else:
            front_base_str=str(str(nucleotide_sequence[position_1+1-self.output_fa_len*3:position_1+1]))
        back_base_str=str(nucleotide_sequence[position_1+len(dup_bases)+1:position_1+len(dup_bases)+1+self.output_fa_len*3])
        str_base_WT=str(yuanshi_nucleotide_sequence)
        str_base_MT=str(nucleotide_sequence)
        ################
        # if weizhi<=self.output_fa_len:
        #     if '*' in protein_sequence[0:weizhi+1]:
        #         front_AA=str(protein_sequence[0:weizhi+1])+'*'
        #     else:
        #         front_AA=str(protein_sequence[0:weizhi+1])
        # elif '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len:weizhi+2]).split('*')[0]+'*')
        # else:
        #     front_AA=str(str(protein_sequence[weizhi-self.output_fa_len:weizhi+2]))
            
        # if '*' in protein_sequence[weizhi-self.output_fa_len-1:weizhi+1]:
        #     back_AA=''
        # elif '*' in protein_sequence[weizhi:weizhi+self.output_fa_len]:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len].split('*')[0]+'*')
        # else:
        #     back_AA=str(protein_sequence[weizhi:weizhi+self.output_fa_len])
        WT_AA=str(yuanshi_nucleotide_sequence.translate()).split('*')[0]+'*'
        if '*' in str(protein_sequence):
            protein_seq=str(protein_sequence)
            protein_seq=protein_seq.split('*')[0]+'*'
            MT_AA=protein_seq
        else : 
            MT_AA=str(protein_sequence)
        #### 2024年08月06日11:15 更新 front_base 和 back_base 的start和end 使用find_difference_range 函数获取
        start,end=find_difference_range(WT_AA,MT_AA)
        if start<=self.output_fa_len:
            front_AA=str(MT_AA[0:start+1])
        else: 
            front_AA=str(MT_AA[start+1-self.output_fa_len:start+1])
        ## back_AA=str(MT_AA[end-1:end+self.output_fa_len-1]) 2024年08月13日 对end的输出进行更改
        if 'fs' in row['AA']:
            ## end =len(MT_AA)
            end = start 
        else:
            end = start 
        back_AA=str(MT_AA[end:end+self.output_fa_len])
        #####################################2024年08月06日11:15 更新结束
    return pd.Series({"front_base": front_base_str, "back_base": back_base_str,'front_AA':front_AA,'back_AA':back_AA,'WT_base':str_base_WT,'MT_base':str_base_MT,'WT_AA':WT_AA,'MT_AA':MT_AA})