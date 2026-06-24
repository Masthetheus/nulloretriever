"""Multiple nullomer data analysis."""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from nulloretriever.analysis.processing import json_to_csv_mapping
# from nulloretriever.analysis.results import *


def setup_argparser() -> argparse.ArgumentParser:
    """Argument passing function for the kmer extraction script."""
    parser = argparse.ArgumentParser(description="""Script for k-mer processing
                                     and nullomer extraction on the snakemake
                                     pipeline""")
    parser.add_argument(
        '-k',
        '--kvalues',
        help="Values of k for k-mer extraction. More than one value can be"
        "informed. Default = 10",
        nargs="*",
        default=[10]
    )
    parser.add_argument(
        '-m',
        '--mode',
        help="Output mode, if in compact txt format or binary format."
        " Default = binary.",
        choices=['binary', 'txt'],
        default='binary'
    )
    parser.add_argument(
        '-o',
        '--output',
        help="Output for the found nullomers file. Default = workflow/data/",
        default='workflow/data/'
    )
    return parser


def main():

    orgdb = pd.read_json('organisms.json')
    print(orgdb.head())
    orgdbt = orgdb.transpose()
    print(orgdbt)
    nulldata = pd.read_csv(
        'nullomer_statistics_summarized.csv')
    print(nulldata.head())
    print(nulldata.columns)
    print(nulldata['counter'])
    df = json_to_csv_mapping(orgdbt, nulldata)
    print(df.columns)
    media = nulldata.groupby('k')['counter'].sum()
    values = nulldata.groupby('k')['counter'].apply(lambda x: (x > 0).sum())
    comp_mean = nulldata.groupby('k')['composition'].mean()
    print(df['speciesname'])
    df['genome_length'] = df['genome_length'].astype(int)
    df.loc[:,'max_poss_kmer'] = df.loc[:,'genome_length'] - df.loc[:,'k'] + 1
    df.loc[:,'null_poss_perc'] = (df.loc[:,'counter']/df.loc[:,'max_poss_kmer'])*100
    df.loc[:,'non_trivial'] = df.loc[:,'counter']-df.loc[:,'trivial']
    df.loc[:,'trivial_porc'] = (df.loc[:,'trivial']/df.loc[:,'counter'])*100
    df.loc[:,'unique_kmers'] = df.loc[:,'max_poss_kmer']-((4**(df.loc[:,'k']))-df.loc[:,'counter'])
    df.loc[:,'unique_per_null'] = df.loc[:,'unique_kmers']/df.loc[:,'counter']
    df.loc[:,'trivial_coverage'] = (df.loc[:,'trivial']/df.loc[:,'genome_length'])*100
    #df.loc[:,'gc_relation'] = (df.loc[:,'gc_perc']-df.loc[:,'composition'])
    df.loc[:,'trivial_non_trivial'] = (df.loc[:,'trivial']/df.loc[:,'non_trivial'])*100
    print(df['non_trivial'])
    df.loc[:,'possible_trivial'] = df.loc[:,'counter']*8
    maxpo = nulldata.groupby('k')['composition'].mean()
    composition_max_count = nulldata[nulldata['composition'] > 1].groupby(
        'k')['composition'].agg(['min', 'max', 'mean'])
    try:
        cpg_mean = nulldata.groupby('k')['motifs_cpg_global_mean'].mean()
        with_cpg_mean = nulldata.groupby(
            'k')['motifs_cpg_mean_nullomers_with_cpg'].mean()
        max_count = nulldata[nulldata['counter'] > 1].groupby(
            'k')['counter'].agg(['min', 'max', 'mean'])
        cpg_max_count = nulldata[nulldata['motifs_cpg_total'] > 1].groupby(
            'k')['motifs_cpg_total'].agg(['min', 'max', 'mean'])
        palindromy_max_count = nulldata[nulldata['motifs_palindromy_count'] > 1].groupby(
            'k')['motifs_palindromy_count'].agg(['min', 'max', 'mean'])
    except Exception as e:
        print("Some stats not found lol")
    grouping_columns = {'k', 'organism', 'organism_name',
                            'tax_id', 'class_id',
                            'assemblystatus', 'taxid', 'speciestaxid'}
    df.to_csv("treated_final_results.csv")
    print(df)
    data_columns = list(set(list(nulldata)) - grouping_columns)
    graph_output = "workflow/results/graphs"
    groups = ["phylum_id","family_id","genus_id","order_id"]
    df_k13 = df[df['k'] == 13]
    sns.scatterplot(data=df_k13, x="genome_length",
                    y="trivial_porc", hue="phylum_id")
    plt.xscale("log")
    plt.xlabel("Genome length (bp, log scale)")
    plt.ylabel("Trivial nullomers (%)")
    plt.savefig(f"{graph_output}/genome_trivial")
    plt.close()
    sns.lineplot(data=df, x="k", y="trivial_porc",
             hue="phylum_id", estimator="median", errorbar="sd")
    plt.savefig(f"{graph_output}/trivial_per_k")
    plt.close()
    sns.boxplot(data=df_k13, x="phylum_id",
            y="trivial_porc")
    plt.xticks(rotation=45)
    plt.savefig(f"{graph_output}/boxplot_phylum")
    plt.close()
    sns.boxplot(data=df_k13, x="phylum_id",
                y="unique_per_null")
    plt.xticks(rotation=45)
    plt.savefig(f"{graph_output}/unique_per_null")
    plt.close()
    sns.scatterplot(data=df_k13, x="unique_kmers",
                    y="counter", hue="phylum_id")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Genome length (bp, log scale)")
    plt.ylabel("Trivial nullomers (%)")
    plt.savefig(f"{graph_output}/unique_per_null_scatter")
    plt.close()

    #for group in groups:
        #sns.scatterplot(data=df, x="genome_length",y="trivial_porc", hue=f"{group}")
        #plt.xscale("log")
        #plt.savefig(f"{graph_output}/{group}_scatter_trivial")
        #plt.close()
        #sns.scatterplot(data=df, x="unique_kmers",y="counter", hue= f"{group}")
        #plt.xscale("log")
        #plt.yscale("log")
        #plt.savefig(f"{graph_output}/{group}_unique_to_counter_ratio")
        #plt.close()
        #sns.scatterplot(data=df, x="genome_length",y="trivial_non_trivial", hue= f"{group}")
        #plt.xscale("log")
        #plt.savefig(f"{graph_output}/{group}_scatter_trivial_non_trivial")
        #plt.close()
        #sns.violinplot(data=df, x="k",y="counter", hue=f"{group}")
        #plt.savefig(f"{graph_output}/{group}_violin_counter")
        #plt.close()
        #sns.boxplot(data=df, x="phylum_id",y="composition", hue=f"{group}")
        #plt.savefig(f"{graph_output}/{group}_composition_box")
        #plt.close()
        # for column in data_columns:
        #     plt.figure(figsize=(12, 6))
        #     sns.lineplot(data=nulldata, x='k', y=column, hue=f"{group}")
        #     plt.savefig(f'{graph_output}/{group}_{column}')
        #     plt.close()
    pass
#    print(media)
#    print(values)
#    print(cpg_mean)
#    print(comp_mean)
#    print(with_cpg_mean)
#    print(max_count)
#    print(cpg_max_count)
#    print(composition_max_count)
#    print(palindromy_max_count)
#    print(nulldata[['k', 'counter']])
#

if __name__ == "__main__":
    main()
