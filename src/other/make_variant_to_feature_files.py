import pandas as pd
import glob
import sys
from vep_parser import parse_annotated_vep_file
from do_regressions import clean_and_filter_variants



def make_files(gene_list):
	for gene in gene_list:
		variant_details_a = parse_annotated_vep_file(vcf_path + gene + ".vcf")
		cleaned_df_a = clean_and_filter_variants(variant_details_a, normalize = False, clean_additional_cols = True)
		cleaned_df_a.to_csv(out_path + gene + ".csv")

if __name__ == "__main__":
	global vcf_path
	global out_path
	global variant_info_path

	# for the broader set...
	#vcf_path = "/net/ukbb/datasets/45878/hail_extracts/gene_annotated_vcfs/"
	vcf_path = "~/ldl_burden_ukbb/annotated_vcfs/"
	variant_info_path = "/net/ukbb/datasets/45878/hail_extracts/gene_files/"
	out_path = "~/lipid_variant_predictions/whole_genome_all_score_files/"

	contained_vcf_genes = glob.glob(vcf_path + "*.vcf")
	contained_vcf_genes = set(
        map(
            lambda x : x.replace(vcf_path, "").replace(".vcf", ""),
            contained_vcf_genes
        )
	)

	print(len(contained_vcf_genes))

	make_files(contained_vcf_genes)
