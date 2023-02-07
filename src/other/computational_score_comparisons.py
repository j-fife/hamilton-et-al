
import glob
import pandas as pd 
import numpy as np 
import pickle
import scipy.stats
import sys

def get_distribution(score_name, genes, use_preload = True):
	if use_preload:
		all_values = pickle.load(open("./score_preloads/" + score_name + ".pickle", "rb"))
		return all_values

	all_values = []
	for i, gene in enumerate(genes):
		print(i, gene)
		if gene in gene_set_2:
			df1 = pd.read_csv("~/ukbb_220k_src/profiles/" + gene + ".csv")
		else:
			df1 = pd.read_csv(path + gene + ".csv")
		
		df1 = df1.set_index("most_severe_variant")
		df2 = pd.read_csv("~/lipid_variant_predictions/whole_genome_all_score_files_2/" + gene + ".csv", index_col = 0)
		df3 = df1[regression_cols].join(df2[regression_cols_9])
		to_predict = df3.loc[df3["Missense"] == 1]
		to_predict = to_predict[[score_name]].dropna() 
		vals = to_predict[score_name].values
		all_values.extend(vals)
	pickle.dump(all_values, open("./score_preloads/" + score_name + ".pickle", "wb" ))


def separate_variants(dist, col_name, low_cutoff = 25, high_cutoff = 75):
	vals = []
	for i, gene in enumerate(genes_of_interest):
		ldl_vals = []
		perc_vals = []

		print(i, gene)
		if gene in gene_set_2:
			df1 = pd.read_csv("~/ukbb_220k_src/profiles/" + gene + ".csv")
			gene_profiles = pd.read_csv("~/ukbb_220k_src/profiles/" + gene + ".csv")
		else:
			df1 = pd.read_csv(path + gene + ".csv")
			gene_profiles = pd.read_csv(path + gene + ".csv")

		df1 = df1.set_index("most_severe_variant")
		df2 = pd.read_csv("~/lipid_variant_predictions/whole_genome_all_score_files_2/" + gene + ".csv", index_col = 0)
		df3 = df1[regression_cols].join(df2[regression_cols_9])
		df3 = df3.loc[df3["Missense"] == 1]
		for index, row in df3.iterrows():
			val = row[col_name]
			percentile = scipy.stats.percentileofscore(dist, val)
			
			carriers = gene_profiles.loc[
				(gene_profiles["most_severe_variant"] == index) & 
				(gene_profiles["ldl"] > 0)
			]

			ldl_vals.extend(carriers["ldl"].values)
			perc_list = [percentile] * len(carriers)
			perc_vals.extend(perc_list)
			

		vals.append(
			[
				gene, 
				col_name, 
				"|".join(list(map(lambda x : str(x), perc_vals))),  
				"|".join(list(map(lambda x : str(x), ldl_vals)))
			]
		)

	df = pd.DataFrame(vals, columns  = ["gene", "score", "percentile_vals", "ldl_vals"])
	df.to_csv("./score_distribution_data/" + col_name + "_ldl_values.csv")
			




if __name__ == "__main__":
	global regression_cols 
	global regression_cols_9
	global upper_cutoff 
	global lower_cutoff 
	global all_genes 
	global gene_set_2
	global genes_of_interest

	genes_of_interest=['PCSK9', 'LDLR', 'APOB', 'NPC1L1', 'ASGR1', 'ABCA1', 'ABCG5', 'TIMD4', 'SLC10A7', 'APOE', 'CFTR', 'HMGCR']
	regression_cols = ["CADD", "GERP", "allele_frequency", "Missense",  "phyloP", "Deleterious"]
	regression_cols_9 = [
		'VEST4_score', 
		'M-CAP_score', 
		'MPC_score', 
		'PrimateAI_score', 
		'GM12878_fitCons_score', 
		'MutationAssessor_score', 
		'FATHMM_score', 
		'MetaLR_score', 
		'MetaSVM_score'
	]
	path = "~/ukbb_220k_src/profiles_whole_genome/"
	all_genes = glob.glob(path + "*.csv")
	
	all_genes = list(map(lambda x : x.replace(path, "").replace(".csv", ""),  all_genes))
	
	gene_set_2 = glob.glob("~/ukbb_220k_src/profiles/*.csv")
	gene_set_2 = list(
		map(
			lambda x : x.replace("~/ukbb_220k_src/profiles/", "").replace(".csv", ""), gene_set_2
		)
	)
	gene_set_2 = set(gene_set_2)

	score_name = sys.argv[1] 

	dist = get_distribution(score_name, all_genes, use_preload = False)
	separate_variants(dist, score_name)

