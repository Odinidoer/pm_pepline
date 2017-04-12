#!/usr/bin/env python
#coding:utf-8
#author:jun.yan@majorbio.com
#last_modified:20170206

import re
import commands

def get_line_num(file):
	line_num = 0
	with open(file,'r') as file_r:
		for line in file_r.readlines():
			if re.search(r'.+?\n',line) and not re.search(r'ccession',line):
				line_num = line_num+1
	return line_num
	
file_inf = commands.getoutput('''for i in *.DE.list;do echo "${i%%.*}"* ;done''')
file_list = file_inf.split('\n')

with open('all_diff_up_down.xls','a')as file_each_a:
	file_each_a.write("name\tall_num\tdiff_num\tup_num\tdown_num\n")
for file in file_list:
	file_each_name = re.search(r'(.*?)\.',file).group(1)
	all = diff = up = down = 0
	file = file.strip('\n')
	file_detail = file.split(' ')
	for file_each in file_detail:
		if re.search('diff\.exp\.xls$',file_each):
		#	print(file_each)
			all = get_line_num(file_each)
		if re.search('DE',file_each):
			diff = diff+get_line_num(file_each)
		if re.search('up',file_each):
			up = up+get_line_num(file_each)
		if re.search('down',file_each):
			down = down+get_line_num(file_each)
	with open('all_diff_up_down.xls','a')as file_each_a:
		file_each_a.write(file_each_name+'\t'+str(all)+'\t'+str(diff)+'\t'+str(up)+'\t'+str(down)+'\n')
