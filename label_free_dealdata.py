#!/usr/bin/env python
#encoding:utf-8
#jun.yan@majorbio.com
#20170310

from itertools import combinations
import argparse
import commands
import os
import re

parser = argparse.ArgumentParser(description="deal label free datas")
parser.add_argument(
    "-i", 
    dest="exp_txt",
    type=str,
    required=True,
    help="please input 'exp.txt'")
parser.add_argument(
    "-g",
    dest="config",
    type=str,
    required=True,
    help="please input 'sample.config'")
parser.add_argument(
    "-o",
    dest="out_dir",
    type=str,
    required=True,
    help="out dir for all results")
args = parser.parse_args()

class label(object):
    def __init__(self, line, sample_config, loc_config, header, con_list):
        self.line = line
        self.sample_config = sample_config
        self.loc_config = loc_config
        self.items = line.split('\t')
        self.combs = list(combinations(con_list, 2))
        self.header = header

    def get_result(self, data):
        data = str(data)
        if re.search(r'\d', data):
            if float(data) == 0:
                return 0
            else:
                return 1
        else:
            return 0

    def get_venn(self):
        for key in self.sample_config.keys():
            venn_file = 'venn/%s.venn.list' % key
            if self.get_result(self.items[int(self.loc_config[key])]):
                with open(venn_file, 'a') as venn:
                    venn.write(self.items[0] + '\n')				

    def get_diff(self):
        for comb in self.combs:
            diff_vs_w = open('diff/%s_vs_%s.exp.txt' %
                             (str(comb[0]), str(comb[1])), 'a')
            diff_list_w = open('diff/%s_vs_%s.diff.exp.xls.list' %
                               (str(comb[0]), str(comb[1])), 'a')
            diff01 = open('diff/%s_vs_%s.01.xls' %
						(str(comb[0]), str(comb[1])), 'a')
            diff10 = open('diff/%s_vs_%s.10.xls' %
						(str(comb[0]), str(comb[1])), 'a')					
            comb0 = [
                s for s in self.sample_config.keys()
                if self.sample_config[s] == comb[0]
            ]
            comb0.sort()
            comb1 = [
                s for s in self.sample_config.keys()
                if self.sample_config[s] == comb[1]
            ]
            comb1.sort()
            get_header_cmd = '''grep "Accession" diff/%s_vs_%s.exp.txt''' %(
            str(comb[0]), str(comb[1]))
            get_header = commands.getoutput(get_header_cmd)
            if not 'Accession' in get_header:
                diff_vs_w.write('Accession\t' + '\t'.join(comb0[
                    x] for x in range(len(comb0))) + '\t' + '\t'.join(comb1[
                        x] for x in range(len(comb1))) + '\n')
						
            get_header_cmd = '''grep "Accession" diff/%s_vs_%s.01.xls''' %(
            str(comb[0]), str(comb[1]))
            get_header = commands.getoutput(get_header_cmd)
            if not 'Accession' in get_header:
                diff01.write('Accession\t' + '\t'.join(comb0[
                    x] for x in range(len(comb0))) + '\t' + '\t'.join(comb1[
                        x] for x in range(len(comb1))) + '\n')	
						
            get_header_cmd = '''grep "Accession" diff/%s_vs_%s.10.xls''' %(
            str(comb[0]), str(comb[1]))
            get_header = commands.getoutput(get_header_cmd)
            if not 'Accession' in get_header:
                diff10.write('Accession\t' + '\t'.join(comb0[
                    x] for x in range(len(comb0))) + '\t' + '\t'.join(comb1[
                        x] for x in range(len(comb1))) + '\n')
			
            for i in range(0,2):
                com = [comb0,comb1][i]
                sum = 0
                locals()['markin%s' %i] = 0
                for key in com:
                    if self.get_result(self.items[int(self.loc_config[key])]):
                        locals()['markin%s' %i] += 1
                        sum += float(self.items[int(self.loc_config[key])])
                if locals()['markin%s' %i] >= len(com) * 2 / 3:
                    locals()['markout%s' %i] = 1
                    for key in com:
                        if not self.get_result(
                                self.items[int(self.loc_config[key])]):
                            self.items[int(self.loc_config[
                                key])] = sum / locals()['markin%s' %i]
                else:
                    locals()['markout%s' %i] = 0						
            if locals()['markout%s' %('0')] and locals()['markout%s' %('1')]:
                diff_vs_w.write(self.items[0] + '\t' + '\t'.join([
                    str(self.items[int(self.loc_config[x])]) for x in comb0
                ]) + '\t' + '\t'.join([
                    str(self.items[int(self.loc_config[x])]) for x in comb1
                ]) + '\n')
            elif locals()['markout%s' %('0')] or locals()['markout%s' %('1')]:				
                diff_list_w.write(self.items[0] + '\n')
                if locals()['markout%s' %('0')]:
					diff10.write(self.items[0] + '\t' + '\t'.join([
                    str(self.items[int(self.loc_config[x])]) for x in comb0
                ]) + '\t' + '\t'.join([
                    str(self.items[int(self.loc_config[x])]) for x in comb1
                ]) + '\n')
                if locals()['markout%s' %('1')]:
					diff01.write(self.items[0] + '\t' + '\t'.join([
                    str(self.items[int(self.loc_config[x])]) for x in comb0
                ]) + '\t' + '\t'.join([
                    str(self.items[int(self.loc_config[x])]) for x in comb1
                ]) + '\n')
				
            diff_vs_w.close()
            diff_list_w.close()
            diff01.close()
            diff10.close()

            com_config = 'diff/' + comb[0] + '_vs_' + comb[1] + '.config'
            with open(com_config, 'w') as ccc:
                for i in range(len(comb0)):
                    ccc.write(comb[0] + '\t' + comb0[i] + '\n')
                for i in range(len(comb1)):
                    ccc.write(comb[1] + '\t' + comb1[i] + '\n')

    def get_cluster(self):
        markout = 0
        for values in set(self.sample_config.values()):
            markin = 0
            keys = [
                s for s in self.sample_config.keys()
                if self.sample_config[s] == values
            ]
            for key in keys:
                if self.get_result(self.items[int(self.loc_config[key])]):
                    markin += 1
            if markin >= len(keys) * 2 / 3:
                markout += 1
        if markout == len(set(self.sample_config.values())):
            cluster_w = open('cluster/All.cluster.txt', 'a')
            get_header_cmd = '''grep "Accession" cluster/All.cluster.txt'''
            get_header = commands.getoutput(get_header_cmd)
            if not 'Accession' in get_header:
                cluster_w.write(self.header)
            cluster_w.write('\t'.join([
                str(self.items[int(i)])
                for i in range(len(self.sample_config.keys()) + 1)
            ]) + '\n')
            cluster_w.close()

config = {}
con_list = []
with open(args.config, 'r') as con:
    for line in con.readlines():
        items = line.strip('\n').split('\t')
        config[items[1]] = items[0]
        if not items[0] in con_list:
            con_list.append(items[0])

if os.path.isdir(args.out_dir):
    cmd = "/bin/rm -rf %s" % args.out_dir
    os.system(cmd)
os.mkdir(args.out_dir)
os.chdir(args.out_dir)
os.mkdir('venn')
os.mkdir('cluster')
os.mkdir('diff')

with open('../' + args.exp_txt, 'r') as exp:
    loc = {}
    headers = exp.readline()
    items = headers.strip('\n').split('\t')
    for i in range(len(items)):
        loc[items[i]] = i
    for line in exp.readlines():
        if re.search(r'.+?\n', line):#blank line
            lafe = label(line.strip('\n'), config, loc, headers, con_list)            
            lafe.get_venn()#must first
            lafe.get_diff()
            lafe.get_cluster()      
