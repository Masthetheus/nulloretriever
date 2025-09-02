"""Functions related to path management and manipulation"""
import os
from pathlib import Path

def gather_files_names(location):
    """Gather the name of all files inside a specified location
    Args:
        location(str): path to search for the files names
    Returns:
        files(arr): array containing all file names
    """
    try:
        location = Path(location)
        files = [f.name for f in location.iterdir() if f.is_file()]
        return files
    except OSError as e:
        print(f"Unable to gather file names, error {e}")
        return []

def gather_files_paths(location):
    """Gather the complete path of all files inside a specified location
    Args:
        location(str): path to search for files
    Returns:
        files_paths(dict): dict containing the full path to all files
    """
    try:
        loc_iter = Path(location)
        files_paths = {}
        files_paths[location] = [f.name for f in loc_iter.iterdir() if f.is_file()]
        return files_paths
    except OSError as e:
        print(f"Unable to gather file path, error {e}")
        return []