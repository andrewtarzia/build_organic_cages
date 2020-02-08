#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Fucntions for cage analysis.

Author: Andrew Tarzia

Date Created: 07 Feb 2020

"""

import json
import numpy as np
import pywindow as pw
from rdkit.Chem import AllChem as Chem


def pw_analyze_cage(
    cage,
    propfile=None,
    structfile=None,
    include_coms=True
):
    """
    Analyze cage already loaded into pywindow.

    """
    # Perform full pyWindow analysis
    cage.full_analysis()
    # Dump pyWindow properties into JSON and cage into xyz
    if propfile is not None:
        cage.dump_properties_json(propfile)
    if structfile is not None:
        cage.dump_molecule(structfile, include_coms=include_coms)


def pw_analyze_cage_from_MOL(
    file,
    prop_out,
    xyz_out,
    include_coms=True
):
    """
    Run all desired analysis on a single built cage molecule.

    Output cage with COM atoms to XYZ and properties to JSON.

    """
    # Import optimised cage into pyWindow, via RDkit mol file
    cage_rd = Chem.MolFromMolFile(file)
    cage_sys = pw.MolecularSystem.load_rdkit_mol(cage_rd)
    cage_mol = cage_sys.system_to_molecule()
    pw_analyze_cage(
        cage=cage_mol,
        propfile=prop_out,
        structfile=xyz_out,
        include_coms=include_coms
    )


def topo_2_property(topology, property):
    """
    Returns properties of a topology for a given topology name.

    Properties:
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
    properties = ['stoich', 'noimines', 'expected_wind']
    if property not in properties:
        raise ValueError(
            f'{property} not defined'
            f'possible properties: {properties}'
        )

    dict = {
        '2p3': {
            'stoich': (2, 3),
            'noimines': 6,
            'expected_wind': 3,
        },
        '4p6': {
            'stoich': (4, 6),
            'noimines': 12,
            'expected_wind': 4,
        },
        '4p62': {
            'stoich': (4, 6),
            'noimines': 12,
            'expected_wind': 4,
        },
        '6p9': {
            'stoich': (6, 9),
            'noimines': 18,
            'expected_wind': 5,
        },
        '1p1': {
            'stoich': (1, 1),
            'noimines': 3,
            'expected_wind': 3,
        },
        '4p4': {
            'stoich': (4, 4),
            'noimines': 12,
            'expected_wind': 6,
        },
    }
    if topology not in dict:
        raise ValueError(f'properties not defined for {topology}')
    return dict[topology][property]


def is_collapsed(topo, pore_diameter, no_window):
    """
    Returns True if a cage is deemed to be collapsed.

    A collapsed cage is defined as having:
        - pore_diam_opt < 2.8 Angstrom (H2 kinetic diameter)
        - number of windows != expected number based on topology.

    """
    expected_wind = topo_2_property(topo, property='expected_wind')
    if expected_wind != no_window:
        return True
    elif pore_diameter < 2.8:
        return True
    else:
        return False


def assign_cage_properties(dct, cage, pw_file):
    '''
    Assigns and outputs arbitrary cage properties to dct.

    '''

    with open(pw_file, 'r') as f:
        pwdata = json.load(f)

    # Get pywindow data.
    dct['max_diam'] = pwdata['maximum_diameter']['diameter']
    dct['p_diam'] = pwdata['pore_diameter']['diameter']
    dct['p_diam_opt'] = pwdata['pore_diameter_opt']['diameter']
    dct['p_vol'] = pwdata['pore_volume']
    dct['p_vol_opt'] = pwdata['pore_volume_opt']
    if pwdata['windows']['diameters'] is None:
        window_data = None
    elif len(pwdata['windows']['diameters']) == 0:
        window_data = None
    else:
        window_data = {
            'w_no': len(pwdata['windows']['diameters']),
            'w_max': max(pwdata['windows']['diameters']),
            'w_min': min(pwdata['windows']['diameters']),
            'w_avg': np.average(pwdata['windows']['diameters'])
        }
    dct['windows'] = window_data

    # Following structural analysis fails if pyWindow bug is
    # encountered. The bug produces overly large windows.
    # Do not do analysis and assume collapsed if big occurs.
    if max(pwdata['windows']['diameters']) < 200:
        # Get window differences if topology is not 4p62.
        if dct['topo'] != '4p62':
            dct['w_diff'] = cage.window_difference()
        else:
            dct['w_diff'] = None

        # Get collapsed flag.
        dct['coll_flag'] = is_collapsed(
            topo=dct['topo'],
            pore_diameter=dct['p_diam_opt'],
            no_window=dct['windows']['w_no']
        )
    else:
        dct['w_diff'] = None
        dct['coll_flag'] = True

    return dct
