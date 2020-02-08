#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build open-form of simple, collapsed imine organic cages.

Author: Andrew Tarzia

Date Created: 08 Feb 2020

"""

import sys
import stk
from os import getcwd, mkdir
from os.path import exists, join

from utilities import read_settings, load_cage_list
import cage_building


def remove_guest(complex, cage, guest):
    """
    Remove guest from cage complex.

    """
    guest_ident = guest.get_identity_key()
    # Determine building block id of host.
    for bb in complex.get_building_blocks():
        ident = bb.get_identity_key()
        if ident != guest_ident:
            print(bb)
            bb_id = bb.id

    # Get host atom ids.
    is_host_atom = []
    for atom in complex.atoms:
        id = atom.id
        print(atom, id)
        if atom.building_block_id == bb_id:
            is_host_atom.append(True)
        else:
            is_host_atom.append(False)

    # Get positon matrix of host in complex.
    pos_mat = complex.get_position_matrix()
    print(pos_mat)
    # Mask pos_mat with host atom ids.
    cage_pos_mat = pos_mat[is_host_atom is True]
    print(cage_pos_mat)
    sys.exit()
    # Update cage position.
    cage.set_position_matrix(cage_pos_mat)
    return cage


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

        print(cage_name)
        print(cd)

        # Read in unoptimised structure.
        unopt_json = f'{cage_name}_unopt.json'
        cage = stk.ConstructedMolecule.load(unopt_json)
        # If output file not present: build cage.
        if not exists(unopt_open_json):
            print(
                '.............building open form of '
                f'{cage_name}..............'
            )
            # No need for guests for certain topologies.
            # In those cages, we just apply a shorter optimisation.
            if cd['topo'] in ['1p1', '2p3']:
                cage_guest = cage.deepcopy()
            else:
                cage_guest = cage_building.build_cage_complex(
                    host=cage,
                    guest=guest
                )
            cage_guest.write(unopt_open_mol)
            cage_guest.dump(unopt_open_json)
        else:
            cage_guest = stk.ConstructedMolecule.load(unopt_open_json)

        sys.exit('view cage guest')

        # If output file not present: optimize cage.
        if not exists(opt_open_json):
            print(
                '.............optimizing open form of '
                f'{cage_name}..............'
            )
            # Optimze cage with guest.
            cage_guest = cage_building.FF_optimize_open_cage(
                name=cage_name,
                cage=cage,
                output_dir=macromod_output,
                macromodel_path=macromod_
            )

            # Write complexed structure out.
            cage_guest.write(opt_compl_mol)
            cage_guest.dump(opt_compl_json)
            sys.exit('view cage guest opt')

            # Remove guest if present.
            if cd['topo'] in ['1p1', '2p3']:
                cage = cage_guest.deepcopy()
            else:
                cage = remove_guest(
                    comlex=cage_guest,
                    cage=cage,
                    guest=guest
                )

            cage.write(opt_open_mol)
            cage.dump(opt_open_json)
            sys.exit('view cage opened opt')

        sys.exit()


if __name__ == "__main__":
    main()
