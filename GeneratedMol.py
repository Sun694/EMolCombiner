import networkx as nx
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

from MolBuilder import MolBuilder


class GeneratedMol:
    def __init__(self, mol_binder_pairs, binder_mol_pairs, base_mols, binder_mols):
        """
        A huge wrapper for molecule-parts, has some helper functions for the required recursion
        This class is agnostic to whether you input bricks or linkers, but comments and
        variable names assume you input bricks
        Args:
            mol_binder_pairs: Set of indices that tells the object how bricks can bind with linkers
            binder_mol_pairs:Set of indices that tells the object how linkers can bind with bricks
            base_mols: The bricks themselves
            binder_mols: The linkers themselves
        """
        self.open_nodes = []
        self.bound_nodes = []
        self.recusion_ready = False

        self.mol_binder_pairs = mol_binder_pairs
        self.binder_mol_pairs = binder_mol_pairs
        self.base_mols = base_mols
        self.binder_mols = binder_mols

    def enable_recursion(self, starter_indice):
        """
        Makes the object recursive.
        Args:
            starter_indice: what brick to use as the seed

        Returns: None

        """
        mol_objs = self.make_base_objs(
            self.mol_binder_pairs, self.base_mols, self.binder_mols
        )
        binder_objs = self.make_base_objs(
            self.binder_mol_pairs, self.binder_mols, self.base_mols
        )

        self.setup_recursive_objs(mol_objs, binder_objs, self.mol_binder_pairs)
        self.setup_recursive_objs(binder_objs, mol_objs, self.binder_mol_pairs)

        self.base_mols = mol_objs
        self.binder_mols = binder_objs

        initial_mol = mol_objs[starter_indice]

        self.mol = initial_mol.mol
        self.open_nodes.append([0, initial_mol])
        self.graph = nx.DiGraph()
        self.graph.add_node(0, mol=initial_mol.mol)
        self.node_index = 1

        self.recursion_ready = True

    def make_base_objs(self, matching_list, base_mol_list, binder_mol_list):
        """
        Binds mols and fragments they can connect to into one object
        Args:
            matching_list: list detailing how the bricks bind with the linkers
            base_mol_list: bricks
            binder_mol_list: linkers

        Returns:

        """
        mol_objs = []
        for i, mol in enumerate(base_mol_list):
            if matching_list[i]:
                b = MolBuilder(mol)
                mol_binding_data = matching_list[i]
                linkers_mols = [
                    binder_mol_list[x] for x in np.array(mol_binding_data).T[0]
                ]
                mol_atom_number = list(np.array(mol_binding_data).T[1])
                binder_atoms_numbers = list(np.array(mol_binding_data).T[2])
                b.initialize_binders(
                    linkers_mols, mol_atom_number, binder_atoms_numbers
                )
                mol_objs.append(b)
            else:
                b = MolBuilder(None)
                mol_objs.append(b)
        return mol_objs

    def setup_recursive_objs(self, base_objs_list, binding_objs_list, matching_list):
        """
        Makes the objects recurse, since a brick can bind to a linker, a linker to a brick, until convergance
        Args:
            base_objs_list: list of brick objects with pairs baked in
            binding_objs_list: list of linker objects with pairs baked in
            matching_list: list detailing how the bricks and linkers bind

        Returns: None

        """
        for i, base_obj in enumerate(base_objs_list):
            if base_obj.mol is not None:
                list_of_binder_objs = []
                binder_indices = np.array(matching_list[i]).T[0]
                base_obj.setup_recursive_binders(
                    list(np.array(binding_objs_list)[binder_indices])
                )

    def bind(self, mol_a, mol_b, mol_a_atom, mol_b_atom):
        """
        Binds the given atoms together
        Args:
            mol_a: first mol
            mol_b: second mol
            mol_a_atom: atom of a to be bound
            mol_b_atom: atom of b to be bound

        Returns: mol with bind

        """
        combined = Chem.CombineMols(mol_a, mol_b)
        editable_combined = Chem.EditableMol(combined)
        editable_combined.AddBond(
            int(mol_a_atom),
            int(mol_a.GetNumAtoms() + mol_b_atom),
            rdkit.Chem.rdchem.BondType.SINGLE,
        )
        bound = editable_combined.GetMol()
        return bound

    def skeletonize_binds(self, index, binder, base_indice, binder_indice):
        """
        Since we abstract the molecule into a graph, actually binding them would make the data structure hard to work with
        Instead, store how they should be bound in a networkx graph, and bind them at the very end
        Args:
            index: index of mol
            binder: index of linker
            base_indice: index of mol atom to bind to linker
            binder_indice: where this atom should bind

        Returns: None

        """
        node_count = self.graph.number_of_nodes()
        self.graph.add_node(
            node_count, mol=binder, parent_atom=base_indice, binder_atom=binder_indice
        )

        self.graph.add_edge(index, node_count)

    def finalize_bonds(self):
        """
        Step through the entire graph, combining all nodes without any children that share a parent
        Returns: generated molecule

        """
        finished = False
        while not finished:
            parent_of_leafs_list = [
                x
                for x in self.graph.nodes()
                if self.graph.out_degree(x) > 0
                and all(
                    [self.graph.out_degree(y) == 0 for y in self.graph.successors(x)]
                )
            ]
            if parent_of_leafs_list:
                for parent_of_leaf in parent_of_leafs_list:
                    leaf_list = list(self.graph.successors(parent_of_leaf))
                    for leaf_indice in leaf_list:
                        leaf = self.graph.nodes.data()[leaf_indice]
                        parent_mol = self.graph.nodes.data()[parent_of_leaf]["mol"]
                        leaf_mol = leaf["mol"]
                        parent_atom_indice = leaf["parent_atom"]
                        leaf_atom_indice = leaf["binder_atom"]

                        bound_mol = self.bind(
                            parent_mol, leaf_mol, parent_atom_indice, leaf_atom_indice
                        )
                        self.graph.add_node(parent_of_leaf, mol=bound_mol)
                        self.graph.remove_node(leaf_indice)
            else:
                finished = True

        cleaned_mol = self.graph.nodes.data()[0]["mol"]
        return cleaned_mol

    def step(self):
        """
        Take a single population step.
        Bind all open atoms to a random uniformly sampled atom that can fit there.
        Returns: True

        """
        assert self.recursion_ready is True

        open_linkers = []
        all_pass = True
        for ind, node in self.open_nodes:
            binders = node.sample_binder()
            if binders is None:
                pass
            else:
                for b in binders:
                    if (
                        sum(
                            [
                                Descriptors.ExactMolWt(x[1]["mol"])
                                for x in self.graph.nodes.data()
                            ]
                        )
                        + Descriptors.ExactMolWt(b[0])
                        < 500
                    ):
                        all_pass = False
                        self.skeletonize_binds(ind, b[0], b[1], b[2])
                        open_linkers.append([self.node_index, b[3]])
                        self.node_index += 1
                    else:
                        pass

        self.open_nodes = []
        [self.open_nodes.append(x) for x in open_linkers]
        if all_pass:
            return True

    def generate(self):
        """
        Skeletonize an entire molecule, but does not finalize graph
        Returns: None

        """
        converged = False
        while not converged:
            converged = self.step()
