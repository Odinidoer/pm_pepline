#!/usr/bin/env python
#coding:utf-8
#author:jun.yan@majorbio.com
#last_modified:20170111

import argparse,sys,os

parser = argparse.ArgumentParser(description = "used to calc and plot correlation between groups")
parser.add_argument("-f","--file",dest = "file",type = str,required = True,help = "please input file,normally 'exp.txt'")
parser.add_argument("-s","--start",dest = "start",type = int,required =	True,help = "choose groups col begin in input file ;113?3:2 ")
parser.add_argument("-t","--type",dest = "type",type = str,default = "pearson",help = "correlation methods: pearson, kendall or spearman; 'pearson' as default")
parser.add_argument("-o","--out",dest = "out",type = str,default = "corrplot",help = "output correlation file: '*.pdf' and '*.txt'; 'corrplot' as default")
args = parser.parse_args()

with open('correlation.r','w') as cor:
	cor.write('''
library("corrplot")
data <-read.table("%s",check.names = F,header = T,sep = "\\t")
corr <- cor(data[%s:dim(data)[2]],method = "%s")#choose method:"pearson", "kendall", "spearman"
write.table(corr, file = "%s.txt", row.names = T, col.names = NA, quote = F, sep = "\\t") 
pdf("%s.pdf")
corrplot(corr,method = "circle", tl.col="black")
dev.off()
''' %(args.file,args.start,args.type,args.out,args.out))

os.system('R --restore --no-save < correlation.r')
#os.remove('correlation.r')
