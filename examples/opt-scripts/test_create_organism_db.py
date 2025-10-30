"""Provides an example to the create_organism_db optional script."""
import json
import yaml

from Bio import Entrez
from nulloretriever.utils.validation import get_valid_email, get_valid_tool
from nulloretriever.data.ncbiapidata import get_genome_metadata

config_file = "../data/config.yaml"
output = "../output/example_organisms_db.json"

print(f"""
Initializing JSON organism db creation via NCBI API using the config file at
{config_file} and outputting the data at {output}.""")

Entrez.email = get_valid_email()
Entrez.tool = get_valid_tool()

with open(config_file, 'r') as f:
    content = yaml.safe_load(f)
    organisms = content['organisms']
metadata = get_genome_metadata(organisms)
print(metadata)
with open(output, 'w') as f:
    json.dump(metadata, f)

print(f"Final example file written, you can check it at {output}.")
