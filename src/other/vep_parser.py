import pandas as pd
import copy
import re
import sys
import numpy as np

def annotate_consequences(consequence_list):
    consequences_used = []
    is_coding = []
    ## getting varaint type and functional consequence
    for c in consequence_list:
        if "splice_acceptor" in c or "frameshift_variant" in c or "splice_donor_variant" in c or "start_lost" in c or "stop_gained" in c:
            consequences_used.append("Deleterious")
            is_coding.append(True)
        elif "missense" in c:
            consequences_used.append("Missense")
            is_coding.append(True)
        elif "inframe_insertion" in c or "inframe_deletion" in c:
            consequences_used.append("Inframe Indel")
            is_coding.append(True)
        elif "synon" in c:
            consequences_used.append("Synonymous")
            is_coding.append(False)
        else:
            consequences_used.append("Intronic")
            is_coding.append(False)
    return consequences_used, is_coding

def determine_most_severe_variant(variant_df):
    deleterious_variants = variant_df.loc[variant_df["Deleterious"] == 1]
    missense_variants = variant_df.loc[variant_df["Missense"] == 1]
    synonymous_variants = variant_df.loc[variant_df["Synonymous"] == 1]

    if len(deleterious_variants) > 0:
        variants_with_cadd = deleterious_variants.loc[~deleterious_variants["CADD"].isna()]
        variants_with_af = deleterious_variants.loc[~deleterious_variants["allele_frequency"].isna()]
        if len(variants_with_cadd) > 0:
            return variants_with_cadd.loc[variants_with_cadd["CADD"] == np.max(variants_with_cadd["CADD"])].iloc[0]

        elif len(variants_with_af) > 0:
            return variants_with_af.loc[variants_with_af["allele_frequency"] == np.min(variants_with_af["allele_frequency"])].iloc[0]

        else:
            return deleterious_variants.iloc[0]

    elif len(missense_variants) > 0:
        variants_with_cadd = missense_variants.loc[~missense_variants["CADD"].isna()]
        variants_with_af = missense_variants.loc[~missense_variants["allele_frequency"].isna()]
        if len(variants_with_cadd) > 0:
            return variants_with_cadd.loc[variants_with_cadd["CADD"] == np.max(variants_with_cadd["CADD"])].iloc[0]

        elif len(variants_with_af) > 0:
            return variants_with_af.loc[variants_with_af["allele_frequency"] == np.min(variants_with_af["allele_frequency"])].iloc[0]

        else:
            return missense_variants.iloc[0]

    elif len(synonymous_variants) > 0:
        variants_with_cadd = synonymous_variants.loc[~synonymous_variants["CADD"].isna()]
        variants_with_af = synonymous_variants.loc[~synonymous_variants["allele_frequency"].isna()]
        if len(variants_with_cadd) > 0:
            return variants_with_cadd.loc[variants_with_cadd["CADD"] == np.max(variants_with_cadd["CADD"])].iloc[0]

        elif len(variants_with_af) > 0:
            return variants_with_af.loc[variants_with_af["allele_frequency"] == np.min(variants_with_af["allele_frequency"])].iloc[0]

        else:
            return synonymous_variants.iloc[0]

    else:
        return variant_df.iloc[0]

def make_vcf(variant_list):
    header = ["#CHROM", "POS","ID","REF", "ALT", "QUAL", "FILTER", "INFO"]
    rows = []
    for v in variant_list:
        pieces = v.split("-")
        chrom = int(pieces[0])
        pos = int(pieces[1])
        ref = pieces[2]
        alt = pieces[3]
        rows.append([chrom, pos, ".", ref, alt, ".", "PASS", "."])
    vcf_df = pd.DataFrame(rows, columns = header)
    vcf_df.to_csv(f"./data/vcf_files/{gene}_SNPs.vcf", sep = "\t", index = False)

def parse_ampersand_value(value_string, position = None):
    if value_string == ".":
        return None 
    if "&" in value_string:
        found = False 
        values = value_string.split("&")
        if position is not None:
            try:
                possible = float(values[position])
                return possible
            except:
                for v in values:
                    try:
                        possible = float(v)
                        return possible
                    except:
                        continue
    else:
        try:
            v = float(value_string)
            return v
        except:
            return None

def parse_annotated_vep_file(
    input_file
):
    score_columns  = [
        'DANN_score', 
        'DEOGEN2_score', 
        'FATHMM_score', 
        'GM12878_fitCons_score', 
        'GenoCanyon_score', 
        'H1-hESC_fitCons_score', 
        'HUVEC_fitCons_score', 
        'LRT_score', 
        'M-CAP_score', 
        'MPC_score', 
        'MVP_score', 
        'MetaLR_score', 
        'MetaSVM_score', 
        'MutPred_score', 
        'MutationAssessor_score', 
        'MutationTaster_score', 
        'PROVEAN_score', 
        'Polyphen2_HDIV_score', 
        'Polyphen2_HVAR_score', 
        'PrimateAI_score', 
        'REVEL_score', 
        'SIFT4G_score', 
        'SIFT_score', 
        'VEST4_score', 
        'fathmm-MKL_coding_score', 
        'fathmm-XF_coding_score', 
        'integrated_fitCons_score',
        'SiPhy_29way_logOdds'
    ]

    rows_to_skip = 0
    canonical_position = 0
    headers = []
    seen_VEP_fields = set()
    info_line_index = None 
    with open(input_file, "r") as f:
        found = False
        line = True
        while not found and line:
            line = f.readline()
            if "Ensembl VEP. Format: " in line:
                info_line_index = rows_to_skip
            if "#CHROM" in line and "REF" in line:
                found = True
                break
            rows_to_skip += 1
            headers.append(line)
    f.close()
    df = pd.read_csv(input_file, delimiter = "\t", skiprows = rows_to_skip)
    info_sections = headers[info_line_index].split("Ensembl VEP. Format: ")[1].replace(">\n", "").split("|")[1:]
    columns = ["Name", "Chrom", "Pos", "Ref", "Alt"] + info_sections
    rows = []
    n = 399
    for index, row in df.iterrows():
        chrom, pos, ref, alt = row["#CHROM"], row["POS"], row["REF"], row["ALT"]
        try: 
            pos = int(row["POS"])
            if chrom not in ["X", "Y"]:
                chrom = int(chrom.replace("chr", ""))
        except:
            continue 
	
        name = str(chrom) + "-" + str(pos) + "-" + ref + "-" + alt
        pieces =  list(map(lambda x: None if x == "" else x, row["INFO"].split("|")[1:]))
        transcripts = [pieces[i:i + n] for i in range(0, len(pieces), n)]
        if len(transcripts) == 0:
            continue
        canonical_fields = None
        canonical_found = 0
        for c_index, c in enumerate(transcripts):
            d = {}
            canon_str = False
            for field, value in zip(info_sections, c):
                if field == "CANONICAL":
                    canon_str = True
                d[field] = value

            if canon_str:
                if (d["CANONICAL"] is not None) and  ("YES" in d["CANONICAL"]):
                    canonical_fields = c
                    canonical_found = 1
                    canonical_position = c_index
                    break

        if canonical_fields is None:
            canonical_fields = transcripts[0]
            canonical_position = 0
        # assert len(pieces) % n == 0
        rows.append([name, chrom, pos, ref, alt] + canonical_fields + [canonical_found])

    variant_df = pd.DataFrame(rows, columns = columns + ["canonical_found"])
    consequences_used, is_coding = annotate_consequences(variant_df["Consequence"].values)
    variant_df["consequence_use"] = consequences_used
    variant_df["is_coding"] = is_coding

    # getting  continuous values

    for score_term in score_columns:
        variant_df[score_term] = list(
            map(
                lambda x : x if x is None else parse_ampersand_value(x, position = canonical_position), 
                variant_df[score_term]
            )
        )


    variant_df["af_use"] = list(map(lambda x : x if x is None else float(x), variant_df["gnomAD_genomes_POPMAX_AF"]))
    variant_df["CADD_use"] = list(map(lambda x : x if x is None else float(x), variant_df["CADD_raw"]))
    variant_df["GERP_use"] = list(map(lambda x : x if x is None else float(x), variant_df["GERP++_RS"]))
    variant_df["phyloP_use"] = list(map(lambda x : x if x is None else float(x), variant_df["phyloP100way_vertebrate"]))
    variant_df["AA_ref_use"] = list(map(lambda x : "Stop" if x == "X" else x, variant_df["aaref"]))
    variant_df["AA_alt_use"] = list(map(lambda x : "Stop" if x == "X" else x, variant_df["aaalt"]))
    variant_df["AA_pos_use"] = list(map(lambda x : x if x is None else int(list(filter(lambda x : x != "", re.split(r'\D+', x)))[0]), variant_df["aapos"]))
    variant_df["coding_position_use"] = list(map(lambda x : x if x is None else int(list(filter(lambda x : x != "", re.split(r'\D+',x)))[0]), variant_df["CDS_position"]))
    ## need CpG context and transition + transversion filters

    columns_keeping = set(
        ['Name', 'Chrom', 'Pos', 'Ref','Alt', "SYMBOL" ,'Consequence',"VEP_canonical",
        "consequence_use","is_coding", "af_use",  "CADD_use", "GERP_use", "phyloP_use", 
        "SiPhy_29way_logOdds", "PROVEAN_score", 
         "AA_ref_use", "AA_alt_use", "AA_pos_use", "coding_position_use",
    ] + score_columns)

    all_columns = list(variant_df.columns)
    for c in columns_keeping:
        all_columns.remove(c)

    variant_df = variant_df.drop(columns = all_columns)
    variant_df.rename(columns = {
        "Consequence" : "vep_consequence",
        "consequence_use" : "consequence",
        "af_use" : "allele_frequency",
        "SYMBOL" : "Symbol",
        "CADD_use" : "CADD",
        "GERP_use" : "GERP",
        "phyloP_use" : "phyloP",
        "AA_ref_use" : "AA_ref",
        "AA_alt_use" : "AA_alt",
        "AA_pos_use" : "AA_pos",
        "coding_position_use" : "coding_position"
    }, inplace = True)

    variant_df = variant_df.set_index("Name")
    return variant_df


def place_position(breaks, pos):
    i = 0
    while pos > breaks[i] and i < len(breaks) - 1:
        i += 1
    return i

def create_precomputed_SNP_df(variant_df, gene):
    #regions = gene_to_region_breaks[gene]
    default_dict = {
         "Deleterious" : 0,
         "Missense" : 0,
         "CADD" : 0,
         "GERP" : 0,
         "phyloP" : 0,
         "allele_frequency" : -2
    }
    rows = []
    for index, row in variant_df.iterrows():
        name = row.name
        d = copy.deepcopy(default_dict)
        d["name"] = name
        if row["consequence"] == "Deleterious":
            d["CADD"] = row["CADD"]
            d["GERP"] = row["GERP"]
            d["phyloP"] = row["phyloP"]
            af = -6
            if row["allele_frequency"] > 0:
                af = np.log10(row["allele_frequency"])
            d["allele_frequency"] = af
            d["Deleterious"] = 1

        elif row["consequence"] == "Missense":
            d["CADD"] = row["CADD"]
            d["GERP"] = row["GERP"]
            d["phyloP"] = row["phyloP"]
            d["Missense"] = 1
            af = -6
            if row["allele_frequency"] > 0:
                af = np.log10(row["allele_frequency"])
            d["allele_frequency"] = af
        elif row["consequence"] == "Synonymous":
            d["CADD"] = row["CADD"]
            d["GERP"] = row["GERP"]
            d["phyloP"] = row["phyloP"]
            d["Synonymous"] = 1
            af = -6
            if row["allele_frequency"] > 0:
                af = np.log10(row["allele_frequency"])
            d["allele_frequency"] = af

        rows.append(d)

    to_return = pd.DataFrame(rows)
    return to_return

		
