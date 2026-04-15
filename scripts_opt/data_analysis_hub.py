"""Generate graphs and csv's of nullomer data analysis."""

import argparse
import pandas as pd
import csv
import matplotlib.pyplot as plt
import seaborn as sns
from nulloretriever.analysis.processing import json_to_csv_mapping


def setup_argparser() -> argparse.ArgumentParser:
    """Parse function for the current script."""
    parser = argparse.ArgumentParser(description="""
    Script that focus on nullomer csv data analysis and manipulation.""")
    parser.add_argument(
        '-db', '--json_db', help='Relative path to the json'
        ' file containing organism data. Default ='
        'workflow/data/organisms.json',
        default="workflow/data/organisms.json"
    )
    parser.add_argument(
        '-d',
        '--data',
        help="Relative path to the csv containing the data to be analyzed."
        "The default takes the range of k present on 'workflow/results/'",
        default="workflow/results/nullomer_statistics_summarized_"
    )
    parser.add_argument(
        '-k',
        '--k_range',
        help="k value range of obtained data. Receives the min and max value"
        ", in this order.",
        nargs=2,
        type=int
    )
    parser.add_argument(
        '-p',
        '--parameters',
        help="Which parameters shall be analyzed. Default is all."
    )
    parser.add_argument(
        '-g',
        '--grouping',
        help="Grouping method for further analysis, up to two parameters."
        "Default = k."
    )
    parser.add_argument(
        '--export',
        help="Exports the grouped and processed data to a csv.",
        action='store_true'
    )
    parser.add_argument(
        '--columns',
        help="States columns used only for grouping or comparison.",
    )
    parser.add_argument(
        '-m',
        '--mode',
        help="Defines if the goal is graph generation(0), csv processing(1)"
        " or both(2).",
        type=int,
        choices=[0, 1, 2],
        default=2
    )
    parser.add_argument(
        '-og',
        '--output_graphs',
        help="Destination for graph storage."
        "Default = workflow/analysis",
        default='workflow/analysis/graphs'
    )
    parser.add_argument(
        '-oc',
        '--output_csvs',
        help="Destination for processed csv storage."
        "Default = workflow/analysis",
        default='workflow/analysis/csv'
    )

    return parser


def main():
    """Execute certain analysis over a csv file and a json organism db."""
    parser = setup_argparser()
    args = parser.parse_args()
    df_db = pd.read_json(args.json_db).transpose()
    k_range = args.k_range
    mode = args.mode
    csv_path = args.data + f"{k_range[0]}_to_{k_range[1]}.csv"
    graph_output = args.output_graphs
    csv_output = args.output_csvs
    if not args.grouping:
        grouping = 'k'
    else:
        grouping = args.grouping
    # Removes all rows where 'column_name' has the value 10
    df = pd.read_csv(csv_path)
    df = df[df['k'] != 14]
    df = df[df['k'] != 13]
    json_to_csv_mapping(df_db, df)
    df_grouped = df.groupby([f'{grouping}'])
    if not args.columns:
        grouping_columns = {'k', 'organism', 'organism_name',
                            'tax_id', 'motifs_homopolymers', 'class_id',
                            'assemblystatus', 'taxid', 'speciestaxid'}
    data_columns = list(set(list(df)) - grouping_columns)

    if mode == 1:
        for column in data_columns:
            plt.figure(figsize=(12, 6))
            sns.lineplot(data=df, x='k', y=column, hue="phylum_id")
            plt.savefig(f'{graph_output}/phylum_id/{column}')
        pass

    elif mode == 2:
        for name, group in df_grouped:
            count = 0
            for column in data_columns:
                plt.figure(figsize=(12, 6))
                current = sns.lineplot(x='k', y=column, hue="phylum_id",
                                       data=group, palette='Paired')
                sns.move_legend(current, "upper right",
                                bbox_to_anchor=(1.05, 1))
                try:
                    if int(name[0]) > 20:
                        plt.savefig(
                            f'{graph_output}/phylum_id/{int(name[0])}_{column}', dpi=200)
                    else:
                        plt.savefig(
                            f'{graph_output}/{name[0]}/{column}', dpi=200)
                except Exception:
                    plt.savefig(
                        f'{graph_output}/phylum_id/{count}_{column}', dpi=200)
                    plt.close()
                    pass
                plt.close()
        df_grouped.describe().to_csv(
            f'{graph_output}/dataset_grouped_by_{grouping}.csv')


if __name__ == "__main__":
    main()
