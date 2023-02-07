
import pandas as pd
import numpy as np
import sys
import time
import json
import glob
import copy
import random
random.seed(0)
import datetime
import pickle
import matplotlib.pyplot as plt
import warnings 
import os
import scipy.stats
from sklearn.linear_model import LinearRegression

def format_p_value_string(p): 
    if p < 0.0001:
        sci_notation_pieces = "{:.2e}".format(p).split("e")
        exponent = int(sci_notation_pieces[1])
        decimal = sci_notation_pieces[0]
        p_string =  "$" + decimal + " x 10^{" + str(exponent) + "}" +  "$"
    else:
        p_string = str(round(p, 4))
    
    return p_string


def get_correlations(gene_list, score_name = "VEST4_score", use_deleterious_variants = True):
    patient_to_adj_ldl = pickle.load(open("~/lipid_variant_predictions/LDL_a_values.pickle", "rb"))
    mean_adj_ldl = np.mean(list(patient_to_adj_ldl.values()))
    valid_patients = set(patient_to_adj_ldl.keys())
    report_vals = []
    score_distribution = pickle.load(open("~/lipid_variant_predictions/score_preloads/" + score_name + ".pickle", "rb" ))
    score_distribution = list(filter(lambda x : x != 0, score_distribution))
    max_score = np.max(score_distribution)
    gene_set_2 = glob.glob("~/ukbb_220k_src/profiles/*.csv")
    gene_set_2 = list(
        map(
            lambda x : x.replace("~/ukbb_220k_src/profiles/", "").replace(".csv", ""), gene_set_2
        )
    )
    gene_set_2 = set(gene_set_2)
    for g in gene_list:
        # new
        if g in gene_set_2:
            p = pd.read_csv("~/ukbb_220k_src/profiles/" + g + ".csv", index_col = 0)
        else:
            p = pd.read_csv("~/ukbb_220k_src/profiles_whole_genome/" + g + ".csv", index_col = 0)
 
        v = pd.read_csv("~/lipid_variant_predictions/whole_genome_all_score_files_2/" + g + ".csv")
        if use_deleterious_variants:
            valid_variants = copy.deepcopy(v)
            delet_variants = set(valid_variants.loc[valid_variants["Deleterious"] == 1]["Name"].values)
        else:
            valid_variants = copy.deepcopy(v.loc[v["Missense"] == 1])
            delet_variants = set()
        covered_variants = set(valid_variants["Name"].values)
        p = p.loc[
            (~p["most_severe_variant"].isna()) &
            (p["ldl"] > 0)
        ]
        x = []
        y = []
        for index, row in p.iterrows():
            var = row["most_severe_variant"]
            if index in valid_patients:
                ldl = patient_to_adj_ldl[index]
            else:
                continue 
            if var in delet_variants:
                x.append(max_score)
                y.append(ldl)
            elif var in covered_variants:
                info = valid_variants.loc[valid_variants["Name"] == var].iloc[0]
                score = info[score_name]
                x.append(score)
                y.append(ldl)
        corr, p = scipy.stats.spearmanr(x, y)
        if use_deleterious_variants:
            col_names = [
                "gene",
                "Spearman correlation Missense + Deleterious",
                "Spearman p-value Missense + Deleterious"
            ]
            out_name = g + "_missense_deleterious.png"
        else:
            col_names = [
                "gene",
                "Spearman correlation Missense",
                "Spearman p-value Missense"
            ]

            out_name = g + "_missense.png"
        report_vals.append([g, corr, p])

        y = list(map(lambda x : x - mean_adj_ldl, y))

        if len(x) == 0:
            print("encountered empty: ", g)
            continue 
        
        fig, (ax1, ax2) = plt.subplots(1,2,figsize = (24, 8))
        y_ticklabels = ["-100","-50","0","+50","+100","+150","+200","+250","+300"]
        y_ticks = list(np.arange(-100, 301, 50))
        x_ticks = list(np.arange(0, 1.01, 0.1))
        x_ticklabels = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        ax = ax1
        ax.scatter(x, y)
        ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), c = "red", lw = 5)
        ax.set_title("$\it{" + g + "}$", fontsize = 25)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_ticklabels, fontsize = 20)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize = 20)
        ax.text(0.05, 275, r"Spearman $ \rho $ = " + str(round(corr, 2)), fontsize = 20)
        ax.text(0.05, 250, r"p-value = " + format_p_value_string(p), fontsize = 20)
        ax.set_xlabel("VEST4 Score", fontsize = 20)
        ax.set_ylabel("$\Delta$ LDL", fontsize = 20)

        ax = ax2 
        bin_size = 0.05
        y_ticklabels = ["-100","-50", "0", "+50", "+100",]
        y_ticks = list(np.arange(-100, 101, 50))
        x_ticks = list(np.arange(0, 1.01, 0.1))
        x_ticklabels = [0, 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        temp_df = pd.DataFrame(np.transpose([x, y]), columns = ["score", "ldl"])
        valid_x = []
        valid_y = []
        sizes = []


        ## final bin getting cut off 
        
        for i in np.arange(0, 1, bin_size):
            valid_x.append(i)
            subset = temp_df.loc[
                (temp_df["score"] >= i) &
                (temp_df["score"] < i + bin_size)
            ]
            valid_y.append(np.mean(subset["ldl"].values))
            sizes.append(len(subset)*5)
        ax.scatter(valid_x,valid_y, s = sizes)
        ax.set_title("$\it{" + g + "}$", fontsize = 25)
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_ticklabels, fontsize = 20)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticklabels, fontsize = 20)
        ax.set_xlabel("VEST4 Score", fontsize = 20)
        ax.set_ylabel("$\Delta$ LDL", fontsize = 20)
        fig.savefig("~/lipid_variant_predictions/VEST4_correlation_plots/" + out_name, bbox_incehes = "tight")
        fig.savefig("~/lipid_variant_predictions/VEST4_correlation_plots/" + out_name.replace(".png", ".eps"), fmt = "eps", bbox_incehes = "tight")
        plt.close(fig)

    df = pd.DataFrame(report_vals,columns = col_names)
    df = df.set_index("gene")
    return df

if __name__ == "__main__":
    global whole_genome_dict
    global max_score 
    # global dist 
    global model
    # model = get_model()
    path = "~/ukbb_220k_src/profiles_whole_genome/"
    with open("~/lipid_variant_predictions/490_genes_from_cell_screen.txt", "r") as f1:
        genes_490 = f1.readlines()
        genes_490 = set(list(map(lambda x : x.replace("\n", ""), genes_490)))
    f1.close()
    all_genes  = glob.glob(path + "*.csv")
    all_genes = list(
        map(
            lambda x : x.replace(path, "").replace(".csv", ""), 
            all_genes
        )
    )
    all_genes_2  = glob.glob("~/ukbb_220k_src/profiles/*.csv")
    all_genes_2 = list(
        map(
            lambda x : x.replace("~/ukbb_220k_src/profiles/", "").replace(".csv", ""), 
            all_genes_2
        )
    )
    all_genes = set(all_genes).union(set(all_genes_2))
    specialized_490 = all_genes.intersection(genes_490)
    t1 = datetime.datetime.now()
    delet_df = get_correlations(specialized_490)
    missense_df = get_correlations(specialized_490, use_deleterious_variants = False)
    large_df = delet_df.join(missense_df)
    large_df.to_csv("VEST4_correlations_screen_genes.csv")

