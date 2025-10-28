"""Nullomer CSV summary generation."""
import pandas as pd

from snakemake.script import snakemake


def main():
    """Return single csv with all nullomer statistics data."""
    out = snakemake.output[0]
    input_files = snakemake.input
    df_list = []
    for file_path in input_files:
        try:
            df = pd.read_csv(file_path, index_col=None, header=0)
            df_list.append(df)
        except FileNotFoundError:
            print(f"Warning! File {file_path} not found! Skipping it")
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    df = pd.concat(df_list, axis=0, ignore_index=True)
    df.to_csv(out)


if __name__ == "__main__":
    """Main function call."""
    main()
