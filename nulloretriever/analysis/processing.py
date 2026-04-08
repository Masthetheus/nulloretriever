"""Functions for nullomer data processing aiming further analysis."""

import pandas as pd


def json_to_csv_mapping(df_json, df_csv):
    """Add all json columns to a csv df."""
    columns_to_add = list(df_json)
    for column in columns_to_add:
        try:
            mapping = df_json.set_index('SpeciesName')[column]
            column = column.lower().replace(" ", "_")
            df_csv[column] = df_csv['organism_name'].map(mapping)
        except Exception:
            continue
    return df_csv
