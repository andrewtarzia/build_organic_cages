#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Fucntions for cage building.

Author: Andrew Tarzia

Date Created: 07 Feb 2020

"""

import stk
from itertools import product


def load_amines(amine_dct):
    # Load in amines.
    for m in amine_dct:
        if amine_dct[m]['method'] == 'file':
            molecule = stk.BuildingBlock.init_from_file(
                amine_dct[m]['file'], ['amine']
            )
        elif amine_dct[m]['method'] == 'smiles':
            molecule = stk.BuildingBlock(
                amine_dct[m]['smiles'], ['amine']
            )
        amine_dct[m]['stk_molecule'] = molecule

    return amine_dct


def load_aldehydes(alde_dct):
    # Load in aldehydes.
    for m in alde_dct:
        if alde_dct[m]['method'] == 'file':
            molecule = stk.BuildingBlock.init_from_file(
                alde_dct[m]['file'], ['aldehyde']
            )
        elif alde_dct[m]['method'] == 'smiles':
            molecule = stk.BuildingBlock(
                alde_dct[m]['smiles'], ['aldehyde']
            )
        alde_dct[m]['stk_molecule'] = molecule

    return alde_dct


def topo_2_property(topology, property):
    """Returns properties of a topology for a given topology name.

    Properties:
        'stk_func' - gives the stk topology function for building cages
        'stoich' - gives the stoichiometries of both building blocks
            assuming that the first building block has the larger
            number of functional groups.
        'noimines' - gives the number of imines formed to build that
            topology
        'expected_wind' - gives the number of windows expected

    Currently defined topologies:
        TwoPlusThree topologies
        ThreePlusThree topologies

    """
    properties = ['stk_func', 'stoich', 'noimines', 'expected_wind']
    if property not in properties:
        raise ValueError(
            f'{property} not defined'
            f'possible properties: {properties}'
        )

    dict = {
        '2p3': {
            'stk_func': stk.TwoPlusThree(),
            'stoich': (2, 3),
            'noimines': 6,
            'expected_wind': 3,
        },
        '4p6': {
            'stk_func': stk.FourPlusSix(),
            'stoich': (4, 6),
            'noimines': 12,
            'expected_wind': 4,
        },
        '4p62': {
            'stk_func': stk.FourPlusSix2(),
            'stoich': (4, 6),
            'noimines': 12,
            'expected_wind': 4,
        },
        '6p9': {
            'stk_func': stk.SixPlusNine(),
            'stoich': (6, 9),
            'noimines': 18,
            'expected_wind': 5,
        },
        'dodec': {
            'stk_func': stk.Dodecahedron(),
            'stoich': (20, 30),
            'noimines': 60,
            'expected_wind': 12,
        },
        '8p12': {
            'stk_func': stk.EightPlusTwelve(),
            'stoich': (8, 12),
            'noimines': 24,
            'expected_wind': 6,
        },
        '1p1': {
            'stk_func': stk.OnePlusOne(
                # place bb1 on vertex (0), bb2 on vertex (1)
                bb_positions={0: [0], 1: [1]}),
            'stoich': (1, 1),
            'noimines': 3,
            'expected_wind': 3,
        },
        '2p2': {
            'stk_func': stk.TwoPlusTwo(
                # place bb1 on vertex (0, 2), bb2 on vertex (1, 3)
                bb_positions={0: [0, 3], 1: [1, 2]}),
            'stoich': (2, 2),
            'noimines': 6,
            'expected_wind': 4,
        },
        '4p4': {
            'stk_func': stk.FourPlusFour(
                # place bb1 on vertex (0, 2), bb2 on vertex (1, 3)
                bb_positions={0: [0, 3, 5, 6], 1: [1, 2, 4, 7]}),
            'stoich': (4, 4),
            'noimines': 12,
            'expected_wind': 6,
        },
    }
    if topology not in dict:
        raise ValueError(f'properties not defined for {topology}')
    return dict[topology][property]


def get_possible_topologies(no_ami_fgs, no_alde_fgs):
    """
    Determine possible topologies based on the FGs in BBs.

    """

    if no_ami_fgs == 2 and no_alde_fgs == 3:
        pos_topo = ['2p3', '4p6', '4p62', '6p9']
    elif no_ami_fgs == 3 and no_alde_fgs == 3:
        pos_topo = ['1p1', '4p4']

    return pos_topo


def determine_cage_list(amine_dct, alde_dct):
    """
    Determines complete cage list to be built from precursors.

    Defines possible topologies based on FGs of each precursor.

    """

    cage_dct = {}

    # Iterate over amine-aldehyde pairs.
    iteration = product(amine_dct, alde_dct)
    for ami_k, alde_k in iteration:
        ami = amine_dct[ami_k]
        alde = alde_dct[alde_k]
        no_ami_fgs = len(ami['stk_molecule'].func_groups)
        no_alde_fgs = len(alde['stk_molecule'].func_groups)
        # Get possible topologies.
        pos_topo = get_possible_topologies(no_ami_fgs, no_alde_fgs)
        for topo in pos_topo:
            cage_name = f"{alde['name']}_{ami['name']}_{topo}"
            cage_dct[cage_name] = {}
            cage_dct[cage_name]['ami'] = ami['name']
            cage_dct[cage_name]['alde'] = alde['name']
            cage_dct[cage_name]['topo'] = topo

    return cage_dct

    return cage_list
