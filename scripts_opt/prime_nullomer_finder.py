"""Compares multiple organisms nullomer tries for primes."""

import argparse
import yaml
import csv
import json
from collections import defaultdict

from nulloretriever.analysis.processing import mount_trie_from_bitfile


def setup_argparser() -> argparse.ArgumentParser:
    """Parsing function for the prime finder script."""
    parser = argparse.ArgumentParser(description="""
                                     Script aimed at searching different
                                     organisms nullomer information for
                                     primes.""")
    parser.add_argument(
        "-c",
        "--config",
        help="Relative path to the YAML config file containing the organisms"
        " names. Default = workflow/config/config.yaml",
        default="workflow/config/config.yaml",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Directory to store the prime sequences file, if"
        " any are found."
        "Default = workflow/results",
        default="workflow/results",
    )
    parser.add_argument(
        "-k",
        "--k_range",
        help="Range of k values to search prime nullomers. Default = 10-12",
        nargs=2,
        type=int,
        default=[10, 12],
    )

    parser.add_argument(
        "-g",
        "--group_by",
        help="Mode of grouping organisms while searching for primes.",
        choices=["none", "phylum_id", "order_id", "family_id", "genus_id"],
        default="none",
    )

    return parser


def dict_flattener(full_dict, final_dict=None, parent_key=""):
    if final_dict is None:
        final_dict = {}

    for key, value in full_dict.items():
        new_key = f"{parent_key}_{key}" if parent_key else key
        if isinstance(value, dict):
            dict_flattener(value, final_dict, parent_key=new_key)
        else:
            final_dict[new_key] = value

    return final_dict


def normalize_accession(acc):
    return acc.replace(".", "_").replace("-", "_")


def main():
    """Search for nullomer primes.
    Outputs:
        {output}prime_null_{group}_{k}: file containing
            the sequences of any found prime nullomer
            for given group in given k value. The root
            directory can be changed via arg parsing.
        prime_range_{k[0]}_to_{k[last]}: file containing
            a summary of found primers statistics grouping
            them by the correspondent k value.
    """
    parser = setup_argparser()
    args = parser.parse_args()
    output = args.output
    with open("workflow/data/organisms.json", "r") as f:
        org_data = json.load(f)
    group = args.group_by
    k_range = args.k_range
    config_path = args.config
    all_stats = []
    org_lookup = {normalize_accession(k): v for k, v in org_data.items()}

    with open(config_path, "r") as f:
        try:
            data = yaml.safe_load(f)
            genomes = data.get("organisms")
        except yaml.YAMLError as exc:
            print(f"An error was found parsing the yaml file: {exc}.")

    grouping = defaultdict(list)
    if group == "none":
        grouping["all"] = genomes
    else:
        for accession in genomes:
            norm = normalize_accession(accession)
            if norm in org_lookup:
                org_group = org_lookup[norm].get(f"{group}", "unknown")
                grouping[org_group].append(accession)

    k = int(k_range[0])
    for group_id, group_genomes in grouping.items():
        k = k_range[0]
        while k <= k_range[1]:
            print(f"======{k}======")
            anchor_trie = mount_trie_from_bitfile(
                f"workflow/results/k{k}/{group_genomes[0]}/null_bit_format"
            )
            orgs_analyzed = len(group_genomes)

            for genome in group_genomes[1:]:
                try:
                    trie2 = mount_trie_from_bitfile(
                        f"workflow/results/k{k}/{genome}/null_bit_format"
                    )
                    anchor_trie.retrieve_prime_null(trie2)
                    print(anchor_trie.count_kmers())
                    if (anchor_trie.count_kmers()) == 0:
                        break
                except Exception as e:
                    print(e)
                    continue

            count = anchor_trie.count_kmers()
            if count > 0:
                retrieved_stats = {
                    "counter": count,
                    "org_count": orgs_analyzed,
                    "gc": anchor_trie.count_gc(),
                    "cpg_stats": anchor_trie.retrieve_nullomers_cpg_stats(),
                    "palindrome_stats": anchor_trie.retrieve_palindrome_stats(),
                }
                base_dict = {
                    "k": k,
                    "group_id": group_id,
                    "group_by": args.group_by,
                }

                row_data = dict_flattener(retrieved_stats, base_dict)
                all_stats.append(row_data)
                if group_id == "N/A":
                    group_id = "not_found"

                identifier = 1
                anchor_trie.write_sequences(f"{output}/prime_null_{group_id}_{k}",identifier)
            k += 1
    if all_stats:
        file_name = f"prime_{k_range[0]}_to_{k-1}.csv"
        header = list(
            {column: True for line in all_stats for column in line.keys()}.keys()
        )
        with open(file_name, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(all_stats)


if __name__ == "__main__":
    main()
