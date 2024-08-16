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
from lib.util import *
from lib.checkPre import *

def main():
    parser = argparse.ArgumentParser(description='Process mutation data.')
    subparsers = parser.add_subparsers(dest="command", help="Sub-commands for different input modes")
    
    # Project subcommand
    parser_project = subparsers.add_parser('snv', help="Input mode for project")
    parser_project.add_argument('--input_file',required=True, help='Input CSV file with mutation data')
    parser_project.add_argument('--df_pkl',required=True,  help='Input PKL file with NM sequence')
    parser_project.add_argument('--output_file',required=True,  help='Output CSV file for processed data')
    parser_project.add_argument('--output_fa_len',required=True, type=int, help='Output peptide_len in csv')
    parser_project.add_argument('--output_peptide_fa_len',required=True, type=int, help='Output peptide_len in fa file')
    parser_project.add_argument('--check_pre', type=bool, default=False,help='get check before run ')
    parser_project.add_argument('--warning_file', help='Optional file to write warnings. If provided, warnings will be written to this file. If not provided, warnings will not be displayed.') 
    
    # File subcommand
    parser_file = subparsers.add_parser('fusion', help="输入fusion的文件结果，输出fusion")
    parser_file.add_argument('--csvfile', required=True, help="Path of the CSV file")
    parser_file.add_argument('--description', help="Description of the project")
    parser_file.add_argument('--output_result', type=str, help="Output directory for the results")
    
    args = parser.parse_args()
    
    if args.command == "snv":
        if args.check_pre==True:
            if checkPree(args):
                handle_snv(args)
            else:
                print('not a csv file')
    elif args.command == "fusion":
        handle_fusion(args)
#-----------------------------------------------------------------------------------------------
def handle_snv(args):
    #### 这里后面需要补个check pre
    ns=NeoantigenSNVfinding(args.input_file,args.df_pkl,args.output_file,args.output_fa_len)
    # ns.mutation_result.apply(ns.call_peptide,axis=1)
    ns.mutation_result=pd.concat([ns.mutation_result,ns.mutation_result.apply(ns.call_peptide,axis=1)],axis=1)
    ns.mutation_result=pd.concat([ns.mutation_result,ns.mutation_result.apply(process_rows, args=(args.output_peptide_fa_len,), axis=1)],axis=1)
    remove_file(os.path.join(ns.file_dir, 'MT.SNV.{len}.fa'.format(len=args.output_peptide_fa_len)))
    remove_file(os.path.join(ns.file_dir, 'MT.INDEL.{len}.fa'.format(len=args.output_peptide_fa_len)))
    ns.mutation_result.apply(slice_seq,args=(args.output_peptide_fa_len,ns.file_dir),axis=1)
    ns.mutation_result.to_pickle(args.output_file+'.pkl')
    ns.mutation_result.to_csv(args.output_file)
    return ns
#-----------------------------------------------------------------------------------------------
def handle_fusion(args):
    #### 这里后面需要补个check pre
    return 

# -------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
    
