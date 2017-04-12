#!/usr/bin/env python
#coding: utf-8
#Program:
#    This program show "to get KOG xls from kog.annot.xls and KOG.class.catalog.xls"!
#History:
#2017/03/27     hui.wan     v1

import re,os, argparse
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='get KOG/COG total information from two excel')

parser.add_argument("-type", "--type", type = str, help='COG/KOG',required=True)
parser.add_argument("-a", "--annotfile", type = str, help = "KOG/COG.annot.xls",required=True)
parser.add_argument("-c", "--classfile", type = str, help = "KOG/COG.class.catalog.xls", required =True)
parser.add_argument("-o", "--outputfile", type = str, help = "KOG/COG.xls", required =True)

args = vars(parser.parse_args())

tp=args["type"]
ann = args["annotfile"]
cla = args["classfile"]
out = args["outputfile"]


annot = pd.read_csv(ann, sep ='\t')
catal = pd.read_csv(cla, sep ='\t')

catal['class']= [i[0] for i in catal.Functional_categories.str.split(' ',n = 1)]
catal['Fuctional_name']= [i[1] for i in catal.Functional_categories.str.split(' ',n = 1)]

annot['class'] = ['[' +i+']' for i in annot.Function.str[1]]

new = pd.merge(annot,catal,on='class')
new.columns = [u'Accession', tp, u'Function', u'Name', u'Class', u'#Type', u'Functional_categories', '_'.join([tp,"stat"]), u'Functional_name']
new.to_csv(out, sep ='\t', columns = ["#Type","Function","Functional_name","Accession",tp, "Name"],index = False,header = True)
