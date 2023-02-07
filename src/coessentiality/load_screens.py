import numpy as np 
import pandas as pd
from sklearn.decomposition import PCA
import statsmodels.api as sm
import os



def clean_positions_df(df):
    gene_symbols = set(df["hg38.kgXref.geneSymbol"])
    rows = []
    for symbol in gene_symbols:
        subset = df.loc[df["hg38.kgXref.geneSymbol"] == symbol]
        strand = subset.iloc[0]["hg38.knownGene.strand"]
        chrom = subset.iloc[0]["hg38.knownGene.chrom"]
        max_len = float("-inf")
        start_pos = None
        end_pos = None
        for index, row in subset.iterrows():
            if abs(row["hg38.knownGene.txStart"] - row["hg38.knownGene.txEnd"]) > max_len:
                max_len = abs(row["hg38.knownGene.txStart"] - row["hg38.knownGene.txEnd"])
                start_pos = row["hg38.knownGene.txStart"]
                end_pos = row["hg38.knownGene.txEnd"]
                
        rows.append([symbol, chrom, start_pos, end_pos, strand])
    to_return = pd.DataFrame(rows, columns = ["gene", "chr", "start", "end", "strand"]).set_index("gene")
    return to_return 

### functions used for corerction as described in Mendillo Lab's FIREWORKs publication
### https://www.life-science-alliance.org/content/4/2/e202000882
### https://github.com/mendillolab/fireworks/blob/master/generate_corr_matrices/generate_corrected_coessential.py
def correct_for_proximity(achilles_fn):
    df_wide = pd.read_csv(achilles_fn, index_col=0) #wide format, cell lines index, gene symbols as columns
    df_wide.columns = list(map(lambda x : x.split(" ")[0], df_wide.columns))
    
    df1 = clean_positions_df(pd.read_csv("./data/UCSC_genome_positions_achilles.txt", delim_whitespace = True))
    df1["Gene"] = df1.index
    df1.rename(columns = {"chr" : "Chromosome", "start" : "Start", "end" : "End"}, inplace = True)
    position_genes = set(list(df1.index))
    
    dgd = pd.read_csv('./data/dgd_Hsa_all_v71.tsv', sep='\t', index_col='Name')
    dup_genes = set(list(dgd.index))
    
    rank_threshold = 20 #on each side
    a = pd.DataFrame()
    for gene in list(df_wide.columns):
        gene = gene.split(" ")[0]
        if gene in dup_genes or gene not in position_genes:
            a[gene] = df_wide[gene] #just copy the original essentiality score if a duplicate gene, or no position UCSC data
        else:
            chromosome = df1[df1['Gene'] == gene].Chromosome.values[0]
            df2 = df1[df1['Chromosome'] == chromosome].sort_values(by='End').reset_index()
            pos = df2[df2['Gene'] == gene].index[0]
            locus_genes = list(df2.loc[pos - rank_threshold:pos+rank_threshold, :].Gene)
            locus_genes.remove(gene)
            a[gene] = df_wide[gene] - (df_wide[[i for i in locus_genes if i in df_wide]].median(axis=1)/2)
    return a


def load_screens(use_changes = True):
    if use_changes:
        screens = correct_for_proximity('./data/Achilles_gene_effect.csv').T
        screens.index = screens.index.str.split(' ').str.get(0)
        cell_lines = pd.read_csv('./data/new_sample_info.csv', index_col='DepMap_ID',
                                     usecols=['DepMap_ID', 'CCLE_name'], squeeze=True)
        annotated_cell_lines = set(cell_lines.index)
        achilles_observed_cell_lines = list(screens.columns)
        achilles_cols_to_drop = list(filter(lambda x : x not in annotated_cell_lines, list(screens.columns)))
        cell_line_cols_to_drop = list(filter(lambda x : x not in achilles_observed_cell_lines, list(cell_lines.index)))
        screens.drop(achilles_cols_to_drop, axis = 1, inplace = True)
        cell_lines.drop(cell_line_cols_to_drop, inplace = True)
        screens.columns = screens.columns.reindex(cell_lines)[0]
        screens = screens.dropna()
        
    else:
        # Load screens
        screens = pd.read_csv('gene_effect.csv', index_col=0).T
        screens.index = screens.index.str.split(' ').str.get(0)
        # Map Broad ID to CCLE name
        cell_lines = pd.read_csv('sample_info.csv', index_col='Broad_ID',
                                 usecols=['Broad_ID', 'CCLE_name'], squeeze=True)
        screens.columns = screens.columns.reindex(cell_lines)[0]
        
    # Bias-correct using "molecular function:olfactory receptor activity" genes
    olfactory_genes = pd.read_csv('olfactory_genes.txt', header=None, squeeze=True)
    olfactory_data = screens.reindex(olfactory_genes).dropna()
    transformation = PCA(n_components=4)
    transformation.fit(olfactory_data)
    top_PC_effects = transformation.inverse_transform(
        transformation.transform(screens))
    screens -= top_PC_effects
    screens = screens.iloc[:, :-4]
    return screens


if __name__ == "__main__":
    #screens = correct_for_proximity('./data/Achilles_gene_effect.csv').T
    #screens.to_csv("Achilles_proximity_corrected.csv")
    ls = load_screens(use_changes = True)
    
