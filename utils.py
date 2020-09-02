from rdkit import Chem


def process_binding_data(raw_bricks, raw_linkers):
    """
    Takes raw SDF files, processes them into a usable format
    Args:
        raw_bricks: brick SDFs
        raw_linkers: linker SDFs

    Returns: Two lists detailing what mols can bind to what, where, and how

    """
    binding_sites_for_bricks = [
        list(
            filter(
                None,
                x[
                    x.index("> <BRANCH @atom-number eligible-atmtype-to-connect> ")
                    + 1 : x.index("> <fragments similar> ")
                ],
            )
        )
        for x in raw_bricks
    ]
    atom_types_for_linkers = [
        list(
            filter(
                None,
                x[
                    x.index("> <MAX-NUMBER-Of-CONTACTS ATOMTYPES> ")
                    + 1 : x.index("$$$$")
                ],
            )
        )
        for x in raw_linkers
    ]

    atom_types_for_linkers = [
        [x.split() for x in sublist] for sublist in atom_types_for_linkers
    ]
    binding_sites_for_bricks = [
        [x.split() for x in sublist] for sublist in binding_sites_for_bricks
    ]

    cleaned_linkers = []
    for linker in atom_types_for_linkers:
        temp_atom_binder_list = []
        for i, atom in enumerate(linker):
            if int(atom[0]) == 0:
                pass
            elif int(atom[0]) > 0:
                temp_atom_binder_list.append(
                    [i, int(atom[0]), atom[1]]
                )  # index, max number of bonds, atom type
            else:
                assert False
        cleaned_linkers.append(temp_atom_binder_list)
    assert len(cleaned_linkers) == len(atom_types_for_linkers)

    brick_linker_pairs = []
    for i, brick in enumerate(binding_sites_for_bricks):
        compatable_linkers_for_brick = []
        for binding_site in brick:
            for j, linker in enumerate(cleaned_linkers):
                for atom in linker:
                    if binding_site[1] == atom[2]:
                        compatable_linkers_for_brick.append(
                            [j, int(binding_site[0]) - 1, atom[0]]
                        )
        brick_linker_pairs.append(compatable_linkers_for_brick)
    # brick_linker_pairs[brick_number] = [linker_index, atom_of_brick_binding_site, atom of linker binding site,]

    linker_brick_pairs = []
    for linker in cleaned_linkers:
        compatable_bricks_for_linker = []
        for bindable_atom in linker:
            for j, brick in enumerate(binding_sites_for_bricks):
                for binding_site in brick:
                    if bindable_atom[2] == binding_site[1]:
                        compatable_bricks_for_linker.append(
                            [j, int(bindable_atom[0]), int(binding_site[0]) - 1]
                        )
        linker_brick_pairs.append(compatable_bricks_for_linker)
    # linker_brick_pairs[linker_number] = [brick_index, atom of linker binding site, atom_of_brick_binding_site]
    return brick_linker_pairs, linker_brick_pairs


def read_all(names, path):
    """
    Read SDFs and return lines
    Args:
        names: all filenames
        path: path to filenames

    Returns: file lines in a list

    """
    raw_sdf = []
    for file in names:
        with open(path + file, "r") as f:
            raw_sdf.append(f.read().splitlines())
    return raw_sdf


def remol_smarts(m):
    """
    Takes an RDKIT Mol, converts it to SMARTS, then back to Mol.
    This fixes a weird bug with vizualization making the object look weird.
    Args:
        m: rdkit Mol

    Returns: rdkit Mol

    """
    return Chem.MolFromSmarts(Chem.MolToSmarts(m))


def remol_smiles(m):
    """
    Takes an RDKIT Mol, converts it to SMARTS, then back to Mol.
    This fixes a weird bug with vizualization making the object look weird.
    Args:
        m: rdkit Mol

    Returns: rdkit Mol

    """
    return Chem.MolFromSmiles(Chem.MolToSmiles(m))
