#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build open-form of simple, collapsed imine organic cages.

Author: Andrew Tarzia

Date Created: 08 Feb 2020

"""
def main():
    """
    Run script.

    """
    usage_string = (
        'build_open_cages.py settings_file cage_list_file guest'
    )
    if (not len(sys.argv) == 4):
        print(f"""
Usage: {usage_string}
    settings_file : str
        File with settings for code.

    cage_list_file :
        .json file to populate with cage definitions of form:
        `name`,`alde`,`ami`,`topology`,`alde`

    guest :
        String to guest file name. This guets will be used to hold the
        cages open.

""")
        sys.exit()
    else:
        settings = read_settings(sys.argv[1])
        cage_list_file = sys.argv[2]
        guest = stk.BuildingBlock.init_from_file(sys.argv[3])

    macromod_ = settings['macromod_']
    macromod_output = join(getcwd(), 'open_data/')
    if not exists(macromod_output):
        mkdir(macromod_output)

    print(settings)

    # Define list of cages to be built.
    # Load this file if it exists to avoid rebuilding cages.
    if exists(cage_list_file):
        # Write a cage list file.
        cage_dct = load_cage_list(cage_list_file)
    else:
        raise FileNotFoundError(
            f'{cage_list_file} not found. Build cages!'
        )

    collapsed_cages = [
        i for i in cage_dct if cage_dct[i]['collapsed'] is True
    ]

    print(
        f'.............{len(collapsed_cages)} '
        ' cages to build........'
    )

    # Build and optimise cages.
    for cage_name in cage_dct:
        cd = cage_dct[cage_name]
        unopt_open_mol = f'{cage_name}_unoptopen.mol'
        unopt_open_json = f'{cage_name}_unoptopen.json'
        opt_compl_mol = f'{cage_name}_optcompl.mol'
        opt_compl_json = f'{cage_name}_optcompl.json'
        opt_open_mol = f'{cage_name}_optopen.mol'
        opt_open_json = f'{cage_name}_optopen.json'

        # Skip if done or not collapsed.
        if cage_name not in collapsed_cages:
            continue

        sys.exit()


if __name__ == "__main__":
    main()
