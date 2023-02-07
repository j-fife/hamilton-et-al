import pandas as pd 
import numpy as np 
import scipy.stats
import gseapy as gp
import glob
import copy
import sys


def combine_results(prefix):
	main_df = None 
	subdirs = glob.glob(prefix + "*")
	for f in subdirs:
		csvs = glob.glob(f + "/*.csv")
		if len(csvs) == 1:
			c = pd.read_csv(csvs[0])
			if main_df is None:
				main_df = c
			else:
				main_df = pd.concat([main_df, c])
	return main_df


def run_gsea(rnk, modules, prefix, use_filter = False):
	gene_members = set(rnk[0].values)
	counter = 0
	gene_sets = {}
	for index, row in modules.iterrows():
		
		members = row["Members"].split(" ")
		members = list(filter(lambda x : x in gene_members, members))
		gene_sets[f"module_{index}"] = members

	gene_sets = {}
	for index, row in modules.iterrows():
		members = row["Members"].split(" ")
		members = list(filter(lambda x : x in gene_members, members))
		if use_filter:
			found = False 
			for m in members:
				if m in screen_genes:
					found = True 
			if found:
				gene_sets[f"module_{index}"] = members
		else:
			gene_sets[f"module_{index}"] = members

	pre_res = gp.prerank(
		rnk=rnk,
		gene_sets = gene_sets, 
		processes=4,
		permutation_num=iterations, # reduce number to speed up testing
		outdir= prefix + "/", 
		format='png', 
		min_size = 4,
		seed=6
    )
