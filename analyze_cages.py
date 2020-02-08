#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build simple imine organic cages.

Author: Andrew Tarzia

Date Created: 07 Feb 2020

"""

import sys
import stk
from os import getcwd, mkdir
from os.path import exists, join

def main():
    """
    Run script.

    """
    usage_string = (
        'analyze_cages.py settings_file '
        'cage_list_file'
    )
    if (not len(sys.argv) == 6):
        print(f"""
Usage: {usage_string}
    settings_file : str
        File with settings for code.

    cage_list_file :
        .json file to populate with cage definitions of form:
        `name`,`alde`,`ami`,`topology`,`alde`

""")
        sys.exit()
    else:
        settings = read_settings(sys.argv[1])
        cage_list_file = sys.argv[2]

    # Define list of cages to be built.
    # Load this file if it exists to avoid rebuilding cages.
    if exists(cage_list_file):
        # Write a cage list file.
        cage_dct = load_cage_list(cage_list_file)
    else:
        raise FileNotFoundError(
            f'{cage_list_file} not found. Build cages!'
        )

if __name__ == "__main__":
    main()
