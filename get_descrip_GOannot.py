#!/usr/bin/env python
#coding: utf-8
#Program:
#    This program show "to add protein descriptions to *gobar.xls"!
#History:
#2017/4/5     hui.wan     v1
#Example: python get_descrip_GOannot.py -i1 e_vs_v.diff.exp.xls.up.listGO.listlevel2-e_vs_v.diff.exp.xls.down.listGO.listlevel2-gobars.xls -i2 protein.xls 

import re,argparse
import pandas as pd
import numpy as np

#参数设置
parser = argparse.ArgumentParser(description='This is to add protein descriptions to *gobar.xls')

parser.add_argument("-i1", "--inputfile1", type=str, help='input *gobar.xls file', required=True)
parser.add_argument("-i2", "--inputfile2", type=str, help='input protein.xls file which contain Accession and Description column names', required=True)
parser.add_argument("-o", "--outputfile", type=str, help = "output gobar.xls", required=False, default = 'gobar.detail.xls')

args = vars(parser.parse_args())
inf1 = args["inputfile1"]
inf2 = args["inputfile2"]
ouf = args['outputfile']

#替换GOid成protein names
def get_description(l):
    if l is not np.nan:
        L = re.findall(r'(.*?)\(.*?\);?',l)
        Lnew =[i+'('+ prod[i]+')' for i in L]
        new = ';'.join(Lnew)
    else:
        new = np.nan
    return new

def get_int(l):
    try: 
	new = int(l)
    except:
	new = np.nan
    return new

#读取protein表格，将accession和description中的蛋白名字存成字典prod
pro = pd.read_csv(inf2,sep='\t',index_col = 0, usecols = ["Accession","Description"])
pros = pro.Description.str.extract(r"(.*?) OS=.*?",expand=True)
prod = pros[0].to_dict()

#替换gobar.xls的GO信息
go = pd.read_csv(inf1,sep = '\t')
old = go[['up proteins','down proteins']]
go["Up proteins"] = old['up proteins'].map(lambda s:get_description(s))
go["Down proteins"] = old['down proteins'].map(lambda s:get_description(s))
go.to_csv(ouf, sep='\t', na_rep= 0, index = False, columns = ["GO term","num of up protein","Up proteins","num of down protein","Down proteins"])

