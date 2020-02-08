#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyze simple imine organic cages.

Author: Andrew Tarzia

Date Created: 07 Feb 2020

"""

import sys
import stk
from os.path import exists

from utilities import load_cage_list, update_cage_json
import cage_analysis


def main():
    """
    Run script.

    """
    usage_string = (
        'analyze_cages.py cage_list_file'
    )
    if (not len(sys.argv) == 2):
        print(f"""
Usage: {usage_string}
    cage_list_file :
        .json file to populate with cage definitions of form:
        `name`,`alde`,`ami`,`topology`,`alde`

""")
        sys.exit()
    else:
        cage_list_file = sys.argv[1]

    # Define list of cages to be built.
    # Load this file if it exists to avoid rebuilding cages.
    if exists(cage_list_file):
        # Write a cage list file.
        cage_dct = load_cage_list(cage_list_file)
    else:
        raise FileNotFoundError(
            f'{cage_list_file} not found. Build cages!'
        )

    print(f'.............{len(cage_dct)} cages to analyse...........')

    # Build and optimise cages.
    for cage_name in cage_dct:
        print(cage_name)
        cd = cage_dct[cage_name]
        opt_mol = f'{cage_name}_opt.mol'
        opt_json = f'{cage_name}_opt.json'
        pw_json = f'{cage_name}_opt.json'
        pw_xyz = f'{cage_name}_opt_pwout.xyz'
        # If analysed, then skip.
        if 'analysed' in cd.keys() and cd['analysed']:
            continue
        else:
            cd['analysed'] = False
        print(cd)

        # If output file of optimization is present, then analyse.
        if exists(opt_json):
            print(f'.............pyWindow analysis of: {cage_name}...')
            cage_analysis.pw_analyze_cage_from_MOL(
                file=opt_mol,
                prop_out=pw_json,
                xyz_out=pw_xyz,
                include_coms=True
            )
        # If pywindow output file is present then save analysis to dct.
        if exists(pw_json):
            # Read in cage structure.
            cage = stk.ConstructedMolecule.load(opt_json)
            cd = cage_analysis.assign_cage_properties(
                dct=cd,
                cage=cage,
                pw_file=pw_json
            )
            cd['analysed'] = True

        print(cd)
        # Update cage dict file.
        cage_dct[cage_name] = cd
        update_cage_json(cage_dct, cage_list_file)
        print(cage_dct)
        sys.exit('check that this all matches files')


if __name__ == "__main__":
    main()
