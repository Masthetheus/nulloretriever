"""Compares multiple organisms nullomer tries for primes."""

from nulloretriever.analysis.prime_nullomer import *
from nulloretriever.analysis.parsings import *
import argparse

def setup_argparser() -> argparse.ArgumentParser:
    """Parsing function for the prime finder script."""
    parser = argparse.ArgumentParser(description="""
                                     Script aimed at searching different
                                     organisms nullomer information for
                                     primes.""")
    return parser

def main():
    """Search for nullomer primes."""
    parser = setup_argparser()
    args = parser.parse_args()
    first_org = "file_placeholder"
    with open(first_org,'r') as f:
        
