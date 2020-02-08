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

from utilities import (
    read_settings,
    read_mol_dict,
    viz_precursors,
    load_cage_list,
    update_cage_json
)
import cage_building


def main():
    """
    Run script.

    """
    usage_string = (
        'build_cages.py settings_file amine_list alde_list '
        'cage_list_file opt_type'
    )
    if (not len(sys.argv) == 6):
        print(f"""
Usage: {usage_string}
    settings_file : str
        File with settings for code.

    amine_list :
        File with amine definitions of form:
        `name`,init (`file` or `smiles`),`-` if file, else SMILES

    alde_list :
        File with aldehyde definitions of form:
        `name`,init (`file` or `smiles`),`-` if file, else SMILES

    cage_list_file :
        .json file to populate with cage definitions of form:
        `name`,`alde`,`ami`,`topology`,`alde`

    opt_type :
        String to determine optimsation settings.
        Options:
            `long`: Recommended (see cage_building.py for settings).
            `short`: (see cage_building.py for settings).

""")
        sys.exit()
    else:
        settings = read_settings(sys.argv[1])
        amine_dct = read_mol_dict(sys.argv[2])
        alde_dct = read_mol_dict(sys.argv[3])
        cage_list_file = sys.argv[4]
        opt_type = sys.argv[5]

    macromod_ = settings['macromod_']
    macromod_output = join(getcwd(), 'data/')
    if not exists(macromod_output):
        mkdir(macromod_output)

    print(settings)

    amine_dct = cage_building.load_amines(amine_dct)
    alde_dct = cage_building.load_aldehydes(alde_dct)

    # Visualise the precursors with rdkit image.
    viz_precursors(dct=amine_dct, out='amine_precursors')
    viz_precursors(dct=alde_dct, out='aldehyde_precusors')

    # Define list of cages to be built.
    # Load this file if it exists to avoid rebuilding cages.
    if exists(cage_list_file):
        # Write a cage list file.
        cage_dct = load_cage_list(cage_list_file)
    else:
        # Determine list of cages to build based on precursor FGs and
        # available toplogies.
        cage_dct = cage_building.determine_cage_list(
            amine_dct=amine_dct,
            alde_dct=alde_dct
        )
        update_cage_json(cage_dct, cage_list_file)

    print(f'.............{len(cage_dct)} cages to build.............')

    # Build and optimise cages.


if __name__ == "__main__":
    main()
