#!/usr/bin/env python
#coding:utf-8
#author:jun.yan@majorbio.com
#last_modified:20170213

import argparse
import commands
import os
import re
import sys

parser = argparse.ArgumentParser(description="change string_data to humaned")
parser.add_argument("-i1", "--DE_list", dest="DE_list", type=str, required=True, help="please input DE.list,normally 'A_vs_B.DE.list'")
parser.add_argument("-i2", "--annotation", dest="annotation", type=str, required=True, help="please input DE.list,normally 'A_vs_B.protein_annotations.tsv'")
parser.add_argument("-i3", "--interaction", dest="interaction", type=str, required=True, help="please input interactions.tsv,normally 'A_vs_B.interactions.tsv'")
parser.add_argument("-o1", "--annotation_out", dest="annotation_out", type=str, help="out annotations.tsv name, default:'A_vs_B.protein_annotations.tsv_new'")
parser.add_argument("-o2", "--interaction_out", dest="interaction_out", type=str, help="out interactions.tsv name, default:'A_vs_B.interactions.tsv_new'")
args = parser.parse_args()
if not args.annotation_out:
	args.annotation_out = '%s_new' %args.annotation
if not args.interaction_out:
	args.interaction_out = '%s_new' %args.interaction	

def get_idmapping(de_list,protein_annotations):
	os.environ['delist'] = de_list
	os.environ['proteinannotations'] = protein_annotations
	cmd = '''awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$1}NR!=FNR{for (i=1;i<=NF;i++){if (a[$i]){mark=i}}{if (mark){print $1,$mark}}}' $delist $proteinannotations'''
	idmapping = commands.getoutput(cmd)
	return idmapping

def set_hash(idmapping):
	items = idmapping.split('\n')
	idm_hash = {}
	for item in items:
		ids = item.split('\t')
		idm_hash[ids[0]] = ids[1]
	return idm_hash
		
idm = get_idmapping(args.DE_list, args.annotation)
print(idm)
idm_hash = set_hash(idm)

new_annotation = open('%s' %args.annotation_out,'a')
new_annotation.write("#node\taccession\tannotation\tdomain_summary_url\n")
with open(args.annotation, 'r')as annotation:
	annotation.readline()
	for line in annotation.readlines():
		node = line.split('\t')[0]
		domain_summary_url = line.split('\t')[2]
		annotation = line.split('\t')[3]
		new_annotation.write('%s\t%s\t%s\t%s\n' %(node, idm_hash[node], annotation, domain_summary_url))
new_annotation.close()

new_interaction = open('%s' %args.interaction_out,'a')
new_interaction.write("#node1\tnode2\tnode1_accession_id\tnode2_accession_id\tneighborhood_on_chromosome\tgene_fusion\tphylogenetic_cooccurrence\thomology\tcoexpression\texperimentally_determined_interaction\tdatabase_annotated\tautomated_textmining\tcombined_score\n")
with open (args.interaction, 'r')as interaction:
	interaction.readline()
	for line in interaction.readlines():
		items = line.split('\t')
		node1 = items[0]
		node2 = items[1]
		new_line = re.sub(r'%s\t%s\t.*?\t.*?\t.*?\t.*?\t' %(node1, node2), r'%s\t%s\t%s\t%s\t' %(node1, node2, idm_hash[node1], idm_hash[node2]), line)
		new_interaction.write(new_line)
new_interaction.close()
