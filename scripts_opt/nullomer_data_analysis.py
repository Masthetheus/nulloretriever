"""Multiple nullomer data analysis."""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
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
        'workflow/results/nullomer_statistics_summarized_9_to_14.csv')
    print(nulldata.head())
    print(nulldata.columns)
    print(nulldata['counter'])
    media = nulldata.groupby('k')['counter'].sum()
    values = nulldata.groupby('k')['counter'].apply(lambda x: (x > 0).sum())
    cpg_mean = nulldata.groupby('k')['motifs_cpg_global_mean'].mean()
    with_cpg_mean = nulldata.groupby(
        'k')['motifs_cpg_mean_nullomers_with_cpg'].mean()
    comp_mean = nulldata.groupby('k')['composition'].mean()
    max_count = nulldata[nulldata['counter'] > 1].groupby(
        'k')['counter'].agg(['min', 'max', 'mean'])
    cpg_max_count = nulldata[nulldata['motifs_cpg_total'] > 1].groupby(
        'k')['motifs_cpg_total'].agg(['min', 'max', 'mean'])
    composition_max_count = nulldata[nulldata['composition'] > 1].groupby(
        'k')['composition'].agg(['min', 'max', 'mean'])
    palindromy_max_count = nulldata[nulldata['motifs_palindromy_count'] > 1].groupby(
        'k')['motifs_palindromy_count'].agg(['min', 'max', 'mean'])
    print(media)
    print(values)
    print(cpg_mean)
    print(comp_mean)
    print(with_cpg_mean)
    print(max_count)
    print(cpg_max_count)
    print(composition_max_count)
    print(palindromy_max_count)


if __name__ == "__main__":
    main()
