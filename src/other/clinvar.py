import pandas as pd 


def load_clinvar_df(fn = "~lipid_variant_predictions/clinvar/variant_summary_2022-05.txt"):
    clinvar_df = pd.read_csv(fn, delimiter="\t")
    return clinvar_df


def get_vep_command(input_vcf, output_vcf):
    vep_command = "vep --cache --force_overwrite --offline --dir /net/data/vep --assembly GRCh38 " 
    vep_command = vep_command + "--plugin dbNSFP,~/vep_annotations_external/build_dbNSFP/dbNSFP4.0a.gz,ALL --vcf "
    vep_command = vep_command + "--no_check_variants_order --canonical --no_stats -I "
    vep_command = vep_command + input_vcf + " -o " + output_vcf
    return vep_command

def make_clinvar_vcf():
    df = load_clinvar_df()
    df = df.loc[
        (df["Assembly"] == "GRCh38") &
        (df["Type"] == "single nucleotide variant")
    ]
   
    rows = []
    header  = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    for index, row in df.iterrows():
        chrom = row["Chromosome"]
        pos = row["Start"]
        ref = row["ReferenceAlleleVCF"]
        alt = row["AlternateAlleleVCF"]
        rows.append([chrom, pos, ".", ref, alt, ".", "PASS", "."])

    df_out = pd.DataFrame(rows, columns = header)
    df_out.to_csv("clinvar_input_to_annotate.vcf", index = False, sep = "\t")

if __name__ == "__main__":
    # make_clinvar_vcf()
    command = get_vep_command("clinvar_input_to_annotate.vcf", "clinvar_annotated.vcf")
    with open("run_vep_on_clinvar.sh", "w") as f:
        f.write(command)
    f.close()