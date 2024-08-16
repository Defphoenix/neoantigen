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
#----------------------------------------------------------------------------------------------------------------

def extract_mutation_info_for_acid(mutation_description):
    """
    从突变描述中提取位置、改变前的碱基和改变后的碱基信息。

    Args:
        mutation_description (str): 突变描述，如 'c.7088G>T'。

    Returns:
        tuple: 包含位置、改变前的碱基和改变后的碱基的信息。
    """
    # 使用正则表达式匹配和提取位置、改变前的碱基和改变后的碱基信息
    mutation_description=mutation_description.replace('*', 'X')
    match = re.match(r'c\.(\d+)([ACGT])>([ACGT])', mutation_description)
    if match:
        position = int(match.group(1))
        original_base = match.group(2)
        mutated_base = match.group(3)
        if mutated_base=='X':
            mutated_base='*'
        return position, original_base, mutated_base
    else:
        return None
    
#----------------------------------------------------------------------------------------------------------------    
def extract_mutation_info_for_indel(mutation_description):
    """
    从突变描述中提取位置、原始氨基酸和突变后的氨基酸信息。

    Args:
        mutation_description (str): 突变描述，如 'p.R809W'。

    Returns:
        tuple: 包含位置、原始氨基酸和突变后的氨基酸的信息。
    """
    match = re.match(r'p\.([A-Z*])(\d+)', mutation_description)
    if match:
        original_aa = match.group(1)
        position = int(match.group(2))
        
        return position 
    else:
        return None
#---------------------------------------------------------------------------------------------------------------- 
def get_gene_info_by_NM(NM_pkl_df,NM_number):
    selected_sequence = NM_pkl_df.loc[NM_pkl_df['NM_ID_no_version'] == NM_number, 'CDS_RNA_Sequence'].values[0]+ NM_pkl_df.loc[NM_pkl_df['NM_ID_no_version'] == NM_number, 'CDS_Rest_of_Sequence'].values[0]
    rna_seq=Seq(selected_sequence)
    rna_seq=rna_seq.back_transcribe()
    return rna_seq

#----------------------------------------------------------------------------------------------------------------

def extract_mutation_info(mutation_description):
    """
    从突变描述中提取位置、原始氨基酸和突变后的氨基酸信息。

    Args:
        mutation_description (str): 突变描述，如 'p.R809W'。

    Returns:
        tuple: 包含位置、原始氨基酸和突变后的氨基酸的信息。
    """
    mutation_description=mutation_description.replace('*', 'X')
    match = re.match(r'p\.([A-Z])(\d+)([A-Z])', mutation_description)
    if match:
        original_aa = match.group(1)
        position = int(match.group(2))
        mutated_aa = match.group(3)
        if mutated_aa == 'X':
            mutated_aa = '*'
        return position, original_aa, mutated_aa
    else:
        return None
#----------------------------------------------------------------------------------------------------------------

def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory '{directory}' created.")
    else:
        print(f"Directory '{directory}' already exists.")

#----------------------------------------------------------------------------------------------------------------

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
    ### 2024年08月06日 更新 因为del的关系，例如SSSSNB 其中最后两个SS丢了 那么start和end会相同，这里强制让end为start+1
    if start==end:
        end=start+1
    return start,end