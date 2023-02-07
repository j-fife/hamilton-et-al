import pandas as pd 
import numpy as np 
import scipy.stats
import glob 
import sys
import pickle 
import json
sys.path.insert(1, "/net/ukbb/datasets/45878/tools/")
from vep_parser import parse_annotated_vep_file
from sklearn.metrics import roc_auc_score

def load_clinvar_df(fn = "~/lipid_variant_predictions/clinvar/variant_summary_2022-05.txt"):
    clinvar_df = pd.read_csv(fn, delimiter="\t")
    return clinvar_df
# deprecated, we no longer use a pretrained model 
def get_model():
    regression_model = pickle.load(open("~/lipid_variant_predictions/regression_model_9_3.pickle", "rb"))
    return regression_model

    
def get_combined_df(use_preload = False):
    if use_preload:
        to_return  = pd.read_csv("~/lipid_variant_predictions/clinvar/combined_clinvar_df_preload.csv", index_col = "Name")
        return to_return
    in_files = glob.glob("~/lipid_variant_predictions/clinvar/splits/annotated*.vcf")
    df_list = []
    for f in in_files:
        print(f)
        df = parse_annotated_vep_file(f)
        df_list.append(df)
    to_return = pd.concat(df_list)
    to_return.to_csv("~/lipid_variant_predictions/clinvar/combined_clinvar_df_preload.csv")
    return to_return 

def get_distribution(name):
    p = "~/lipid_variant_predictions/score_preloads/" + name + ".pickle"
    dist = pickle.load(open(p, "rb"))
    return dist 

def get_whole_genome_dict():
    variant_to_prediction = pickle.load(open("~/ukbb_220k_src/src/whole_genome_dict_checkpoints/checkpoint_completed_032521_9.pickle", "rb" ))
    return variant_to_prediction


def auc_delet_model(clinvar_df):
    model = get_model()
    variant_to_prediction = get_whole_genome_dict()
    dist = list(variant_to_prediction.values())
    annotated_subset = clinvar_annotated_df
    covered_variants = set(annotated_subset.index)
    needed_variants = set(clinvar_df.index)
    path_variants = set(clinvar_df.loc[
        clinvar_df["is_path"] 
    ].index)
    valid_variants = covered_variants.intersection(needed_variants)
    annotated_subset = annotated_subset.loc[list(valid_variants)]
    y_pred = []
    y_true = []
    annotated_subset["allele_frequency"] = annotated_subset["allele_frequency"].map(lambda x : -6 if pd.isna(x) else np.log10(x))
    annotated_subset["Missense"] = 1
    annotated_subset["Deleterious"] = 0
    annotated_subset = annotated_subset[regression_cols]
    annotated_subset = annotated_subset.dropna()
    annotated_subset = annotated_subset[~annotated_subset.isin([np.nan, np.inf, -np.inf]).any(1)]

    for counter, (index, row) in enumerate(annotated_subset.iterrows()):
        score = model.predict([row.values])
        pct = scipy.stats.percentileofscore(dist, score)
        y_pred.append(pct)
        annot = 1 if index in path_variants else 0
        y_true.append(annot)

    auc = roc_auc_score(y_true, y_pred)
    return auc, y_pred, y_true   

def get_labeled_variants():
    clinvar_df = load_clinvar_df()
    pathogenic_annotations = ["Pathogenic"]
    benign_annotations = ["Benign"]
    acceptable_annotations = pathogenic_annotations + benign_annotations
    clinvar_df = clinvar_df.loc[
        (clinvar_df["Assembly"] == "GRCh38") &
        (clinvar_df["Type"] == "single nucleotide variant") &
        (clinvar_df["ClinicalSignificance"].isin(acceptable_annotations))
    ]
    clinvar_df["is_path"] = clinvar_df["ClinicalSignificance"].apply(lambda x : x in pathogenic_annotations)
    return clinvar_df 


def calculate_auroc(score_name, clinvar_df):
    annotated_subset = clinvar_annotated_df.loc[
        ~clinvar_annotated_df[score_name].isna()
    ]
    covered_variants = set(annotated_subset.index)
    needed_variants = set(clinvar_df.index)

    path_variants = set(clinvar_df.loc[
        clinvar_df["is_path"] 
    ].index)
    valid_variants = covered_variants.intersection(needed_variants)
    annotated_subset = annotated_subset.loc[list(valid_variants)]
    y_pred = []
    y_true = []
    distribution = get_distribution(score_name)
    for counter, (index, item) in enumerate(annotated_subset[score_name].iteritems()):
        pct = scipy.stats.percentileofscore(distribution, item)
        y_pred.append(pct)
        annot = 1 if index in path_variants else 0
        y_true.append(annot)
    auc = roc_auc_score(y_true, y_pred)
    return auc, y_pred, y_true


if __name__ == "__main__":
    global clinvar_annotated_df 
    global regression_cols 
    additional_cols = [
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
    regression_cols = ["CADD", "GERP", "allele_frequency", "Missense", "phyloP", "Deleterious"] + additional_cols
    clinvar_annotated_df = get_combined_df(use_preload = True)
    clinvar_annotated_df = clinvar_annotated_df.loc[
        clinvar_annotated_df["vep_consequence"] == "missense_variant"
    ]
    clinvar_df = get_labeled_variants()
    clinvar_df["Name"] = clinvar_df["Chromosome"].astype(str) + "-" + clinvar_df["PositionVCF"].astype(str) + "-" + clinvar_df["ReferenceAlleleVCF"].astype(str) + "-" + clinvar_df["AlternateAlleleVCF"].astype(str)
    clinvar_df = clinvar_df.set_index("Name")
    score_name = sys.argv[1]

    if score_name == "Delet":
        auc, y_pred, y_true = auc_delet_model(clinvar_df)
    else:
        auc, y_pred, y_true = calculate_auroc(score_name, clinvar_df)
    
    to_return = {
        "auc" : auc,
        "y_pred" : y_pred,
        "y_true" : y_true
    }

    with open("~/lipid_variant_predictions/clinvar/auc_results/" + score_name + ".json", "w") as o:
        json.dump(to_return, o)
    o.close()



    
