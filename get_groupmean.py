#!/usr/bin/env python
#coding: utf-8
#Program:
#    This program show "to averaging the exp.txt according sample.config"!
#History:
#2017/03/28     hui.wan     v1

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='averaging the exp.txt according sample.config')

parser.add_argument("-e", "--expfile", type=str, help = "exp.txt", required=True)
parser.add_argument("-s", "--sample", type= str, help = "sample.config", required =True)
parser.add_argument("-o", "--output", type= str, help = "outfile", required =False, default = "group.exp.txt")



args = vars(parser.parse_args())

exp = args["expfile"]
spc = args["sample"]
out = args["output"]

Exp = pd.read_csv(exp, sep = '\t', index_col = 0)
Spc = pd.read_csv(spc, sep = '\t', header=None, index_col =1,names=["group"])

exp_group = Exp.groupby(Spc.group.to_dict(),axis = 1).mean()

exp_group.to_csv(out,sep = '\t', header = True, index = True)
