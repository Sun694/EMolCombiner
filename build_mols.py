import argparse
import os
import sys

import networkx as nx
import numpy as np
import rdkit.Chem.Draw
from rdkit import Chem
from rdkit import RDConfig

from GeneratedMol import GeneratedMol
from utils import read_all, process_binding_data, remol_smiles

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))

parser = argparse.ArgumentParser(
    description="Probabilistically generate new molecules based off fragments."
)

parser.add_argument(
    "--input_directory", type=str, help="emolfrag output dir", required=True
)
parser.add_argument(
    "--output_directory",
    type=str,
    help="emolcombiner output dir",
    required=False,
    default="emolcombinerout/",
)
parser.add_argument(
    "--num_mols",
    type=int,
    help="number of mols to generated",
    default=100,
    required=False,
)
parser.add_argument(
    "--draw",
    help="whether to save output mols",
    default=False,
    action="store_true",
    required=False,
)

args = parser.parse_args()

base_path = args.input_directory + "\\"
path_to_bricks = base_path + "output-brick/"
path_to_linkers = base_path + "output-linker/"

list_of_bricks = os.listdir(path_to_bricks)
list_of_linkers = os.listdir(path_to_linkers)

mol_bricks = [Chem.MolFromMolFile(path_to_bricks + x) for x in list_of_bricks]
mol_linkers = [Chem.MolFromMolFile(path_to_linkers + x) for x in list_of_linkers]

raw_bricks = read_all(list_of_bricks, path_to_bricks)
raw_linkers = read_all(list_of_linkers, path_to_linkers)

brick_linker_pairs, linker_brick_pairs = process_binding_data(raw_bricks, raw_linkers)


def gen(i):
    """
    Generates a single molecule under 500MW
    Args:
        i: current iterator, determines file-name

    Returns: i+1 if generation was successful, i otherwise

    """
    G = GeneratedMol(brick_linker_pairs, linker_brick_pairs, mol_bricks, mol_linkers)
    G.enable_recursion(np.random.randint(0, len(mol_bricks) - 1))
    try:
        G.generate()
        nx.draw(G.graph)
        out_mol = G.finalize_bonds()
        i += 1
        Chem.MolToMolFile(
            out_mol, args.output_directory + "{}.mol".format(i), kekulize=False
        )
        if args.draw:
            rdkit.Chem.Draw.MolToFile(
                remol_smiles(out_mol), args.output_directory + "ims/{}.png".format(i)
            )
    except AttributeError as e:
        pass
    return i


try:
    os.mkdir(args.output_directory)
    os.mkdir(args.output_directory + "ims/")
except:
    print("directory already exists, continuing")

i = 0
while i < args.num_mols:
    i = gen(i)
