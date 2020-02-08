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
from os.path import join


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


def define_imine_topology(bb1, bb2, topology, property):
    """
    Define known topology graphs.

    This requires modification for new cage types.

    Currently defined topologies:
        TwoPlusThree topologies
        ThreePlusThree topologies

    """
    # Define the topologies.
    topo_func = None
    bb_vect = None
    if property == 'stk_func':
        if topology == '2p3':
            topo_func = stk.cage.TwoPlusThree()
            bb_vect = None
        elif topology == '4p6':
            topo_func = stk.cage.FourPlusSix()
            bb_vect = None
        elif topology == '4p62':
            topo_func = stk.cage.FourPlusSix2()
            bb_vect = None
        elif topology == '6p9':
            topo_func = stk.cage.SixPlusNine()
            bb_vect = None
        elif topology == '1p1':
            topo_func = stk.cage.OnePlusOne()
            bb_vect = {
                bb1: [topo_func.vertices[0]],
                bb2: [topo_func.vertices[1]],
            }
        elif topology == '4p4':
            topo_func = stk.cage.FourPlusFour()
            bb_vect = {
                bb1: topo_func.vertices[0, 2, 4, 6],
                bb2: topo_func.vertices[1, 3, 5, 7],
            }

    return topo_func, bb_vect


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


def build_cage(bb1, bb2, topology, bb_vert=None):
    """
    Build cage with stk.

    """
    cage = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2],
        topology_graph=topology,
        building_block_vertices=bb_vert
    )
    return cage


def build_cage_complex(host, guest):
    """
    Build cage complex with stk.

    """
    complex = stk.ConstructedMolecule(
        building_blocks=[host, guest],
        topology_graph=stk.host_guest.Complex()
    )
    return complex


def FF_optimize_cage(
    name,
    cage,
    settings,
    output_dir,
    macromodel_path
):
    """
    Optimize cage with stk.

    """

    ff = stk.MacroModelForceField(
        macromodel_path=macromodel_path,
        output_dir=join(output_dir, f'{name}_FFout'),
        restricted=True
    )
    # MD process - run MD, collect N conformers, optimize each,
    # return lowest energy conformer.
    md = stk.MacroModelMD(
        macromodel_path=macromodel_path,
        output_dir=join(output_dir, f'{name}_MDout'),
        timeout=settings['timeout'],
        force_field=settings['force_field'],
        temperature=settings['temperature'],
        conformers=settings['conformers'],
        time_step=settings['time_step'],
        eq_time=settings['eq_time'],
        simulation_time=settings['simulation_time'],
        maximum_iterations=settings['maximum_iterations'],
        minimum_gradient=settings['minimum_gradient'],
        use_cache=settings['use_cache']
    )
    seq = stk.Sequence(ff, md)
    seq.optimize(mol=cage)
    return cage


def FF_optimize_open_cage(
    name,
    complex,
    output_dir,
    macromodel_path
):
    """
    Crudely optimize cage complex with stk.

    """

    ff_res = stk.MacroModelForceField(
        macromodel_path=macromodel_path,
        output_dir=join(output_dir, f'{name}_FFResout'),
        restricted=True
    )
    ff_unres = stk.MacroModelForceField(
        macromodel_path=macromodel_path,
        output_dir=join(output_dir, f'{name}_FFUnresout'),
        restricted=False
    )
    seq = stk.Sequence(ff_res, ff_unres)
    seq.optimize(mol=complex)
    return complex


def get_settings(opt_type):
    """
    Define settings for MacroModel optimization in stk.

    Options:
        `default`:
            Default settings from stk source code as of 26/04/19.

        `long`:
            My settings for rigorous cage optimizations in stk.
            Mimics: Computationally-inspired discovery of an
            unsymmetrical porous organic cage - DOI:10.1039/C8NR06868B
            Modified on 26/04/19.
            Modified on 06/06/19.

        `short`:
            My default settings for short cage optimizations in stk.
            Modified on 26/04/19.

    """

    setting_options = {
        'default': {
            'output_dir': None, 'timeout': None,
            'force_field': 16, 'temperature': 300,  # K
            'conformers': 50, 'time_step': 1.0,  # fs
            'eq_time': 10, 'simulation_time': 200,  # ps
            'maximum_iterations': 2500, 'minimum_gradient': 0.05,
            'use_cache': False
        },
        'long': {
            'output_dir': None, 'timeout': None,
            'force_field': 16, 'temperature': 700,  # K
            # change from 10000
            'conformers': 5000, 'time_step': 0.5,  # fs
            'eq_time': 100,  # ps
            # ps -- 50 ns changed from 100 ns
            'simulation_time': -500,
            'maximum_iterations': 2500, 'minimum_gradient': 0.05,
            'use_cache': False
        },
        'short': {
            'output_dir': None, 'timeout': None,
            'force_field': 16, 'temperature': 700,  # K
            'conformers': 50, 'time_step': 1,  # fs
            'eq_time': 50, 'simulation_time': 1000,  # ps -- 1 ns
            'maximum_iterations': 2500, 'minimum_gradient': 0.05,
            'use_cache': False
        }
    }

    if opt_type not in list(setting_options.keys()):
        raise ValueError(
            f'{opt_type} is not a defined option'
            f'Options incl: {setting_options.keys}'
        )

    return setting_options[opt_type]
