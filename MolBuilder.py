import numpy as np


class MolBuilder:
    def __init__(self, molecule):
        """
        Helper class to decide what can bind where
        Args:
            molecule: rdkit Mol object that will be the seed
        """
        self.mol = molecule
        self.bound = False

    def initialize_binders(
        self, list_of_molecules, molecule_atom_numbers, binder_atoms_numbers
    ):
        """
        Store what molecules can bind where and how
        Args:
            list_of_molecules: list of mols
            molecule_atom_numbers: what atom of the seed they can bind to
            binder_atoms_numbers: what atom of the seed binds to the mol

        Returns: None

        """
        assert len(list_of_molecules) == len(molecule_atom_numbers)
        assert len(list_of_molecules) == len(binder_atoms_numbers)
        self.binders = list_of_molecules
        self.molecule_atom_numbers = molecule_atom_numbers
        self.binder_atoms_numbers = binder_atoms_numbers

    def setup_recursive_binders(self, binder_objects):
        """
        Enables recursion for easier generations
        Args:
            binder_objects: other Mol wrappers

        Returns: None

        """
        self.binder_objects = binder_objects

    def sample_binder(self):
        """
        Randomly sample an atom to bind to each open indice of seed
        Returns: List that details how to bind if an atom is available, None if no open spots

        """
        if not self.bound:
            (
                chosen_binders,
                chosen_mol_nums,
                chosen_binder_atoms,
                chosen_binder_objects,
            ) = ([], [], [], [])
            self.bound = True
            for unique_node_atom in np.unique(np.array(self.molecule_atom_numbers).T):
                possible_linker_indices = np.arange(len(self.molecule_atom_numbers))[
                    np.where(
                        np.array(self.molecule_atom_numbers).T == unique_node_atom,
                        True,
                        False,
                    )
                ]
                chosen_index = np.random.choice(possible_linker_indices)

                if self.binder_atoms_numbers[chosen_index] in np.unique(
                    self.binder_objects[chosen_index].molecule_atom_numbers
                ):
                    chosen_binders.append(self.binders.pop(chosen_index))
                    chosen_mol_nums.append(self.molecule_atom_numbers.pop(chosen_index))

                    chosen_binder_atom_indice = self.binder_atoms_numbers.pop(
                        chosen_index
                    )
                    chosen_binder_atoms.append(chosen_binder_atom_indice)
                    chosen_binder_object = self.binder_objects.pop(chosen_index)
                    chosen_binder_object.remove_bindable_atom(chosen_binder_atom_indice)

                    chosen_binder_objects.append(chosen_binder_object)

            return list(
                np.array(
                    [
                        chosen_binders,
                        chosen_mol_nums,
                        chosen_binder_atoms,
                        chosen_binder_objects,
                    ]
                ).T
            )
        else:
            return None

    def remove_bindable_atom(self, atom_indice):
        """
        Marks an atom of a given Mol as "taken", to prevent multiple bonds on a single atom where not allowed
        Args:
            atom_indice: what atom to mark as taken

        Returns: None

        """

        assert atom_indice in self.molecule_atom_numbers

        indices_to_keep = np.arange(len(self.binders))[
            np.where(np.array(self.molecule_atom_numbers) == atom_indice, False, True)
        ]

        self.binders = list(np.array(self.binders)[indices_to_keep])
        self.molecule_atom_numbers = list(
            np.array(self.molecule_atom_numbers)[indices_to_keep]
        )
        self.binder_atoms_numbers = list(
            np.array(self.binder_atoms_numbers)[indices_to_keep]
        )
        self.binder_objects = list(np.array(self.binder_objects)[indices_to_keep])
