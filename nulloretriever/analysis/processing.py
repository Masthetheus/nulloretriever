"""Functions for nullomer data processing aiming further analysis."""

import pandas as pd


def json_to_csv_mapping(df_json, df_csv):
    """Add all json columns to a csv df."""
    columns_to_add = list(df_json)
    print(columns_to_add)
    print(df_json.index())
    for column in columns_to_add:
        try:
            mapping = df_json.set_index('e')[column].replace(".","_")
            print("\n")
            print(mapping)
            column = column.lower().replace(" ", "_")
            df_csv[column] = df_csv['organism'].map(mapping)
        except Exception:
            print(f"Colum {column} Nao foi")
            continue
    return df_csv
