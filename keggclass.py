#!/usr/bin/env python
#coding: utf-8
#Program:
#    This program show "to class KEGG annotation and enrichment results according organisms"!
#History:
#2017/1/16     hui.wan     v2
import argparse,re,os,shutil,linecache,commands
import rpy2.robjects as robjects    #robjects.r['Rcode']
from rpy2.robjects import r 

parser = argparse.ArgumentParser(description='To extract kegg and pathways of specific organism')

parser.add_argument("-i", "--inputfile", help = "input file or folder",required=True)
parser.add_argument("-orgC", "--orgnismClass", help = "organism classification: animals/plants/bacteria/all",default="all",required=False)
parser.add_argument("-annot", "--annotation", help = "do annotation:T/F",required=True)
parser.add_argument("-enrich", "--enrichment", help = "do enrichment: T/F",required=True)

args = vars(parser.parse_args())
inf = args["inputfile"]
org=args["orgnismClass"]
annot=args["annotation"]
enrich=args["enrichment"]

#check organism argparse 
organs = ['all','animals','plants','bacteria']
try:
    if org in organs:
        org='.'.join([org,'ko_20161230.list'])
        orgPath='/'.join(['/mnt/ilustre/users/hui.wan/script/wanhui/keggOrgClass',org])
        with open(orgPath,'r') as koList:
            pathko=[line.split('\t')[0] for line in koList.readlines()]
            ko=[s.split(':')[1] for s in pathko]
            
except:
    raise IOError("%s is a incorrect input organism,just put animals or plants or bacteria." %orgnismClass)
    exit


    

class File2file(object):
    def __init__(self,oldfilepath,newfilepath):
        self.old = oldfilepath
        self.new = newfilepath
        
#match can use or not use ,multiMatch function can replace it 
    def match(self,x):
        with open(self.old,'r') as old:
            for line in old.readlines():
                if x in line:
                    with open(self.new,'a+') as new:
                        if not line in new.readlines():
                            new.write(line)
                            
    def multiMatch(self):
        with open(self.old,'r') as old:
            for line in old.readlines():
                ko_old_all = re.findall(r'ko\d+', line)
                ko_old_new= filter(lambda x:x in ko,ko_old_all)
                with open(self.new,'a+') as new:
                    newline = line.replace(";".join(["path:"+i for i in ko_old_all]),";".join(["path:"+i for i in ko_old_new]))
                    if not newline in new.readlines():
                        new.write(newline)
                            
    def addheader(self,nrow):
        header=linecache.getlines(self.old)[0:nrow]
        with open(self.new,'a') as f:
            f.write(''.join(header))
            

if annot == "T":
    oldfr = inf
    #creat new folder and filename
    oldpath="/".join([oldfr,'pathway_table.xls'])
    oldkegg="/".join([oldfr,'kegg_table.xls'])
    newfr=".".join([oldfr,'filter'])
    newpath="/".join([newfr,'pathway_table.xls'])
    newkegg="/".join([newfr,'kegg_table.xls'])
        
    #output new folder is exit?
    if os.path.exists(newfr):
        shutil.rmtree(newfr)
        os.mkdir(newfr)
    else:
        os.mkdir(newfr)
    
    #deal excel
    path=File2file(oldpath,newpath)
    kegg=File2file(oldkegg,newkegg)
    path.addheader(1)
    kegg.addheader(1)
    map(path.match,ko)
    #path.multiMatch()
    kegg.multiMatch()

    #deal png and html file
    def getpic(x):
        listdir=os.listdir(oldfr)
        for file in listdir:
            if x in file:
                oldpic="/".join([oldfr,file])
                newpic="/".join([newfr,file])
                shutil.copyfile(oldpic,newpic)
    map(getpic,ko)


if enrich == "T":
    oldfile = inf
    newfile=".".join([oldfile,"filter"])
    tempfile = ".".join([oldfile,"temp"])
    newfr = []
    enrichtemp = File2file(oldfile,tempfile)
    map(enrichtemp.match,ko)
    
    #recaculate fdr
    with open(tempfile,'r') as f:
        fr = [line.strip() for line in f.readlines()]
        pval = [i.split('\t')[5] for i in fr]
        plen = len(pval)
        FDR = list(robjects.r['p.adjust'](pval,method="fdr",n=plen))
        FDRold = [i.split('\t')[6] for i in fr]
        for i in range(len(FDRold)):
            newfr.append(fr[i].replace(FDRold[i],str(FDR[i])))
    os.remove(tempfile)
    
    #write newfile
    if os.path.isfile(newfile):
        os.remove(newfile)
    enrich = File2file(oldfile,newfile)
    enrich.addheader(5)
    with open(newfile,'a') as f:
        f.write('\n'.join(newfr))
        f.write('''
--------------------

#Term	Database	ID	Input number	Background number	P-Value	Corrected P-Value	Input	Hyperlink

--------------------

#Term	Database	ID	Input number	Background number	P-Value	Corrected P-Value	Input	Hyperlink

--------------------
''')

