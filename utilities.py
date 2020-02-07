#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utility functions.

Author: Andrew Tarzia

Date Created: 07 Feb 2020

"""


from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import json


def load_cage_list(file):
    """
    Read .json file containing cage list.

    """

    with open(file, 'r') as f:
        cage_list = json.load(f)

    return cage_list


def update_cage_json(dct, file):
    """
    Update .json file containing cage list.

    """

    with open(file, 'w') as f:
        json.dump(dct, f)


def read_settings(file):
    """
    Read settings.

    Returns
    -------
    params : :class:`dict`
        Dictionary of parameters.

    """
    params = {}
    with open(file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        key, val = line.rstrip().split(':')
        if val == 'True' or val == 'False':
            val = True if val == 'True' else False
        else:
            try:
                val = float(val)
            except ValueError:
                pass
        params[key] = val

    return params


def read_mol_list(file):
    """
    Read list of molecules from file.

    Returns
    -------
    molecules : :class:`list` of :class:`dict`
        List of molecule definitions.

    """
    with open(file, 'r') as f:
        lines = f.readlines()

    molecules = []
    for line in lines:
        molecule = {}
        name, method, smiles = line.rstrip().split(',')
        molecule['name'] = name
        molecule['method'] = method
        molecule['smiles'] = smiles
        molecule['file'] = f'{name}.mol'
        molecules.append(molecule)

    return molecules


def mol_list2grid(
    molecules,
    filename,
    mol_per_row,
    maxrows,
    subImgSize=(200, 200),
    names=None
):
    """
    Produce a grid of molecules in mol_list.

    molecules (list) - list of molecule SMILEs

    """

    if len(molecules) > mol_per_row * maxrows:
        print('here')
        # have to make multiple images
        new_mol_list = []
        new_names = []
        count = 1
        for i, mol in enumerate(molecules):
            print('here', names[i], mol)
            new_mol_list.append(mol)
            if names is None:
                new_names = None
            else:
                new_names.append(names[i])
            # make image
            if len(new_mol_list) == mol_per_row * maxrows:
                img = Draw.MolsToGridImage(
                    new_mol_list,
                    molsPerRow=mol_per_row,
                    subImgSize=subImgSize,
                    legends=new_names,
                    useSVG=True
                )
                save_svg(
                    filename=f'{filename}_{count}.svg',
                    string=img
                )
                # img.save(filename + '_' + str(count) + '.png')
                new_mol_list = []
                new_names = []
                count += 1
    else:
        img = Draw.MolsToGridImage(
            molecules,
            molsPerRow=mol_per_row,
            subImgSize=subImgSize,
            legends=names,
            useSVG=True
        )
        save_svg(filename=f'{filename}.svg', string=img)
        # img.save(filename + '.png')


def save_svg(filename, string):
    """
    Save svg text to a file.

    """

    with open(filename, 'w') as f:
        f.write(string)


def viz_precursors(lst, out):
    """
    Make RDKit image of lst of molecules.

    """

    mol_list = [
        Chem.MolFromSmiles(
            Chem.MolToSmiles(Chem.MolFromMolFile(i['file']))
        )
        for i in lst
    ]
    mol_names = [i['name'] for i in lst]
    mol_list2grid(
        molecules=mol_list,
        names=mol_names,
        filename=out,
        mol_per_row=4,
        maxrows=5
    )
