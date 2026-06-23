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

    orgdb = pd.read_json('workflow/data/organisms.json')
    print(orgdb.head())
    orgdbt = orgdb.transpose()
    print(orgdbt)
    nulldata = pd.read_csv(
        'workflow/results/nullomer_statistics_summarized.csv')
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
    df.loc['null_poss_perc'] = (df.loc[:,'counter']/df.loc[:,'max_poss_kmer'])*100
    df.loc[:,'non_trivial'] = df.loc[:,'counter']-df.loc[:,'trivial']
    df.loc[:,'trivial_porc'] = (df.loc[:,'trivial']/df.loc[:,'counter'])*100
    df.loc[:,'unique_kmers'] = df.loc[:,'max_poss_kmer']-((4**(df.loc[:,'k']))-df.loc[:,'counter'])
    df.loc[:,'trivial_coverage'] = (df.loc[:,'trivial']/df.loc[:,'genome_length'])*100
    df.loc[:,'trivial_non_trivial'] = (df.loc[:,'trivial']/df.loc[:,'non_trivial'])*100
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
    except:
        print("Some stats not found lol")
    grouping_columns = {'k', 'organism', 'organism_name',
                            'tax_id', 'class_id',
                            'assemblystatus', 'taxid', 'speciestaxid'}
    df.to_csv("treated_final_results.csv")
    print(df)
    data_columns = list(set(list(nulldata)) - grouping_columns)
    graph_output = "workflow/results/graphs"
    sns.scatterplot(data=df, x="genome_length",y="trivial_porc", hue="phylum_id")
    plt.xscale("log")
    plt.savefig(f"{graph_output}/scatter_trivial")
    plt.close()
    sns.scatterplot(data=df, x="unique_kmers",y="counter", hue="phylum_id")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(f"{graph_output}/unique_to_counter_ratio")
    plt.close()
    sns.scatterplot(data=df, x="genome_length",y="trivial_non_trivial", hue="phylum_id")
    plt.xscale("log")
    plt.savefig(f"{graph_output}/scatter_trivial_non_trivial")
    plt.close()
    sns.violinplot(data=df, x="k",y="counter", hue="phylum_id")
    plt.savefig(f"{graph_output}/violin_counter")
    plt.close()
    sns.boxplot(data=df, x="phylum_id",y="composition", hue="phylum_id")
    plt.savefig(f"{graph_output}/composition_box")
    plt.close()
    for column in data_columns:
        plt.figure(figsize=(12, 6))
        sns.lineplot(data=nulldata, x='k', y=column, hue="phylum_id")
        plt.savefig(f'{graph_output}/{column}')
        plt.close()
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
