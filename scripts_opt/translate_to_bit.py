"""From a nullomer sequence file obtains the correspondent compact bit file"""

import argparse
import yaml

from pathlib import Path
from nulloretriever.utils.bit_encoding import write_bit_file

def setup_argparser() -> argparse.ArgumentParser:
    """Argument passing function for the encoder script."""
    parser = argparse.ArgumentParser(description="""Script for bit file encoding from nullomeric sequences.""")
    parser.add_argument(
        "-k",
        "--kvalues",
        help="Values of k from the source files. If the files aren't sorted in foldes by k value, more than one k value will result in error."
        "Default = 10",
        nargs="*",
        default=[10],
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output directory for the bit files. Default = workflow/results",
        default="workflow/results",
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input directory with sequence files to convert to compact bit format. Default = workflow/results",
        default = "workflow/results",
    )
    parser.add_argument(
        "-f",
        "--filename",
        help="Full path of the nullomer file to be encoded. Used in the case of single file conversion."
    )
    parser.add_argument(
        "--multiple",
        action = "store_true",
        help = "Allows for the processing of multiple sequence files. When this flag is used. all files inside the input directory are read and outputted in compact binary format."
    )
    return parser


def main():
    """Converts pure sequence nullomer file to the compact bit version from nulloretriever."""
    parser = setup_argparser()
    args = parser.parse_args()
    dir_input = args.input
    dir_output = args.output
    filename = args.filename
    k_values = args.kvalues
    multiple = args.multiple

    if not dir_input and not filename:
        print("Please inform a config file, filename or directory with the files to be analyzed")
        sys.exit(1)
    if filename and k_values:
        if len(k_values) == 1:
            full_out_path = f"{filename}_bit_k{k_values[0]}"
            write_bit_file(filename, full_out_path, int(k_values[0]))
            print("Sequence file converted to its binary compact format.")
        else:
            for k in k_values:
                try:
                    dir_k = f"{dir_input}/k{k}/{filename}"
                    full_out_path = f"{dir_output}/{filename}_bit_k{k}"
                    write_bit_file(dir_k, full_out_path, int(k))
                except Exception as e:
                    print(f"Something went wrong. Check the pathing inputs. Error: {e}")
            print("All files correctly outputted in compact binary format.")
    elif multiple:
        for k in k_values:
            try:
                dir_k = Path(f"{dir_input}/k{k}")
                for file in dir_k.iterdir():
                    if file.is_file():
                        print(f"Writing file {file.stem} for k {k}.")
                        full_out_path = f"{dir_output}/{file.stem}_bit_format_k{k}"
                        write_bit_file(file, full_out_path, int(k))
            except Exception as e:
                print(f"Failure writing the binary output. Please check if the informed path exists. Error: {e}.")
                break

        print("All informed sequences correctly converted to the compact bit format!")
    else:
        print("Please select the operating format. Filename for one organism conversion or use the --multiple flag for converting all sequences file inside the input directory.")

if __name__ == "__main__":
    main()
