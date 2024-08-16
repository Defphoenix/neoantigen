##2024年06月17日 更新
import argparse
import pandas as pd
import sys
import re
import warnings
import os
import pdb
import math
# 如果提供了 --warning_file 参数，则捕获警告信息但不输出到 stderr

#----------------------------------------------------------------------------------------------------------------
## waning信息存储
def warning_handler(message, category, filename, lineno, file=None, line=None):
    # 如果提供了 --warning_file 参数，则写入警告信息
    if args.warning_file:
        with open(args.warning_file, 'a') as f:
            log_entry = f"{category.__name__}: {message}\n"
            f.write(log_entry)
            
# 删除指定文件（如果存在）
def remove_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)
    # 创建一个新的空文件
    with open(file_path, 'w') as file:
        pass  # 使用 pass 来创建一个空文件
#----------------------------------------------------------------------------------------------------------------
def main(args):
    mutation_result=pd.read_pickle(args.input_file)
    if len(mutation_result['type']=='snv') >0:
        snv_df_mutilist=mutation_result[mutation_result['type']=='snv']
        snv_df_mutilist=snv_df_mutilist[['chr','ref','alt','gene','type','protein_variant_type_annovar','transcript','base','AA']]
        snv_df_mutilist['type']=='SNV'
    if len(mutation_result['type']=='indel') >0:
        indel_df_mutilist=mutation_result[mutation_result['type']=='indel']
        indel_df_mutilist=indel_df_mutilist[['chr','ref','alt','gene','type','protein_variant_type_annovar','transcript','base','AA']]
        indel_df_mutilist['type']=='Indel'
    # pdb.set_trace()
    remove_file(os.path.join(args.output_dir,args.sample_name+'.snv.mutlist'))
    remove_file(os.path.join(args.output_dir,args.sample_name+'.indel.mutlist'))
    snv_df_mutilist.to_csv(os.path.join(args.output_dir,args.sample_name+'.snv.mutlist'),sep='\t',header=False,index=False)
    indel_df_mutilist.to_csv(os.path.join(args.output_dir,args.sample_name+'.indel.mutlist'),sep='\t',header=False,index=False)
    # 使用示例
#----------------------------------------------------------------------------------------------------------------
warnings.showwarning = warning_handler
#----------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process mutation data.')
    parser.add_argument('--input_file', help='Input CSV file mutation data with peptide')
    parser.add_argument('--sample_name', help='Output fa file sample_name')
    parser.add_argument('--output_dir', help='Output fa file peptide_len')
    parser.add_argument('--warning_file', help='Optional file to write warnings. If provided, warnings will be written to this file. If not provided, warnings will not be displayed.') 
    args = parser.parse_args()

    # 运行你的 main 函数
    mutation_result=main(args)


