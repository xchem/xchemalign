import gemmi
from loguru import logger

# import loguru
from xchemalign.data import (
    Atom,
    AtomID,
    Dataset,
    LigandID,
    LigandNeighbourhood,
    Structure,
    SystemData,
)


def get_structure_fragments(
    dataset: Dataset, structure: Structure
) -> dict[LigandID, gemmi.Residue]:
    fragments: dict[LigandID, gemmi.Residue] = {}
    lig_number: int = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.name == "LIG":
                    ligand_id: LigandID = LigandID(
                        dtag=dataset.dtag, id=lig_number
                    )
                    fragments[ligand_id] = residue
                    lig_number = lig_number + 1

    return fragments


def get_model_and_artefact_atoms(
    residue_neighbours: dict[
        tuple[float, float, float], gemmi.NeighborSearch.Mark
    ],
    structure: Structure,
) -> tuple[list[gemmi.NeighborSearch.Mark], list[gemmi.NeighborSearch.Mark]]:
    # Check each mark for its image and partition them on this
    model_atoms: list[gemmi.NeighborSearch.Mark] = []
    artefact_atoms: list[gemmi.NeighborSearch.Mark] = []
    for pos, mark in residue_neighbours.items():
        # Image 0 is the identity i.e. part of the normal model
        cra = mark.to_cra(structure[0])
        logger.debug(f"### CRA: {cra}")
        logger.debug(f"Mark pos: {mark.pos()}")
        logger.debug(f"Canonical atom pos: {cra.atom.pos}")
        logger.debug(f"Image idx: {mark.image_idx}")

        nearest_image = structure.cell.find_nearest_pbc_image(
            mark.pos(), cra.atom.pos, mark.image_idx
        )
        logger.debug(f"{nearest_image}")
        logger.debug(f"{nearest_image.sym_idx}")
        logger.debug(f"{nearest_image.pbc_shift}")

        if nearest_image.sym_idx != 0:
            artefact_atoms.append(mark)
        else:
            model_atoms.append(mark)

    return model_atoms, artefact_atoms


def get_ligand_neighbourhood(
    structure: Structure,
    ns: gemmi.NeighborSearch,
    fragment: gemmi.Residue,
    min_dist: float = 0.01,
    max_dist: float = 5.0,
) -> LigandNeighbourhood:
    # For each atom, get the neighbouring atoms, and filter them on their
    # real space position
    residue_neighbours: dict[
        tuple[float, float, float], gemmi.NeighborSearch.Mark
    ] = {}
    _artefact_atoms = []
    _model_atoms = []
    model_atoms: dict[AtomID, Atom] = {}
    artefact_atoms: dict[AtomID, Atom] = {}

    for atom in fragment:
        atom_neighbours: list[gemmi.NeighborSearch.Mark] = ns.find_neighbors(
            atom,
            min_dist=min_dist,
            max_dist=max_dist,
        )
        for neighbour in atom_neighbours:
            cra = neighbour.to_cra(structure[0])
            nearest_image = structure.cell.find_nearest_pbc_image(
                atom.pos, cra.atom.pos, neighbour.image_idx
            )
            logger.debug(f"{nearest_image}")
            logger.debug(f"{nearest_image.sym_idx}")
            logger.debug(f"{nearest_image.pbc_shift}")

            fpos = structure.cell.fractionalize(cra.atom.pos)
            logger.debug(f"--FPos: {fpos}")

            ftransform = ns.get_image_transformation(neighbour.image_idx)
            logger.debug(f"--Transform: {ftransform}")

            fpos_transformed = ftransform.apply(fpos)
            fpos.x = fpos_transformed + nearest_image.pbc_shift[0]
            fpos.y = fpos_transformed + nearest_image.pbc_shift[1]
            fpos.z = fpos_transformed + nearest_image.pbc_shift[2]

            logger.debug(f"--Transformed FPos: {fpos_transformed}")

            pos = structure.cell.orthogonalize(fpos_transformed)
            # logger.debug(
            #     f"--Transformed pos: {pos} vs Original pos: {atom.pos}"
            # )

            # nearest_image_dist = nearest_image.dist()
            # dist = atom.pos.dist(pos)
            # logger.debug(f"--Distance: {dist} vs nid {nearest_image_dist}")

            rounded_pos = (
                round(pos.x, 1),
                round(pos.y, 1),
                round(pos.z, 1),
            )

            residue_neighbours[rounded_pos] = cra

            # if nearest_image.sym_idx != 0:
            #     artefact_atom_id: AtomID = AtomID(
            #         chain=cra.chain.name,
            #         residue=cra.residue.seqid.num,
            #         atom=cra.atom.name,
            #     )
            #     artefact_atoms[artefact_atom_id] = Atom(
            #         element=cra.atom.element.name,
            #         atom_id=artefact_atom_id,
            #         x=pos.x,
            #         y=pos.y,
            #         z=pos.z,
            #     )
            #     # artefact_atoms[rounded_pos] = pos
            # else:
            #     # _model_atoms.append(neighbour)
            #     model_atom_id: AtomID = AtomID(
            #         chain=cra.chain.name,
            #         residue=cra.residue.seqid.num,
            #         atom=cra.atom.name,
            #     )
            #     model_atoms[model_atom_id] = Atom(
            #         element=atom.element.name,
            #         atom_id=model_atom_id,
            #         x=atom.x,
            #         y=atom.y,
            #         z=atom.z,
            #     )

    # exit()
    # logger.debug(f"Found {len(residue_neighbours)} atoms near residue")

    # # Seperate out model and artefact atoms
    # _model_atoms, _artefact_atoms = get_model_and_artefact_atoms(
    #     residue_neighbours, structure
    # )
    # logger.debug(f"Got {len(_model_atoms)} model atoms")
    # logger.debug(f"Got {len(_artefact_atoms)} artefact atoms")

    # Model atoms
    model_atoms: dict[AtomID, Atom] = {}
    for atom in _model_atoms:
        cra = atom.to_cra(structure[0])
        model_atom_id: AtomID = AtomID(
            chain=cra.chain.name,
            residue=cra.residue.seqid.num,
            atom=cra.atom.name,
        )
        model_atoms[model_atom_id] = Atom(
            element=atom.element.name,
            atom_id=model_atom_id,
            x=atom.x,
            y=atom.y,
            z=atom.z,
        )

    # Artefact atoms
    artefact_atoms: dict[AtomID, Atom] = {}
    for atom in _artefact_atoms:
        artefact_cra = atom.to_cra(structure[0])
        artefact_atom_id: AtomID = AtomID(
            chain=artefact_cra.chain.name,
            residue=artefact_cra.residue.seqid.num,
            atom=artefact_cra.atom.name,
        )
        artefact_atoms[artefact_atom_id] = Atom(
            element=atom.element.name,
            atom_id=artefact_atom_id,
            x=atom.x,
            y=atom.y,
            z=atom.z,
        )

    # Cosntruct the neighbourhood
    logger.debug(model_atoms)
    logger.debug(artefact_atoms)
    ligand_neighbourhood: LigandNeighbourhood = LigandNeighbourhood(
        atoms=model_atoms,
        artefact_atoms=artefact_atoms,
    )

    return ligand_neighbourhood


def get_dataset_neighbourhoods(
    dataset: Dataset, max_radius: float = 5.0
) -> dict[LigandID, LigandNeighbourhood]:
    # Load the structure
    structure: Structure = gemmi.read_structure(dataset.pdb)
    logger.debug(f"{structure.cell}")

    # Get the bound fragments
    fragments: dict[LigandID, gemmi.Residue] = get_structure_fragments(
        dataset, structure
    )
    logger.debug(f"Get {len(fragments)} fragment neighbourhoods")
    logger.debug(fragments)

    # Construct the neighbourhood search
    ns: gemmi.NeighborSearch = gemmi.NeighborSearch(
        structure[0], structure.cell, max_radius
    ).populate()

    # For each bound fragment, identify the neighbourhood atoms and
    # partition them into model and artefact
    fragment_neighbourhoods: dict[LigandID, LigandNeighbourhood] = {}
    for ligand_id, fragment in fragments.items():
        fragment_neighbourhoods[ligand_id] = get_ligand_neighbourhood(
            structure, ns, fragment
        )

    return fragment_neighbourhoods


def get_ligand_neighbourhoods(
    system_data: SystemData,
) -> dict[LigandID, LigandNeighbourhood]:
    # Iterate over data, loading in structures, getting ligands for each
    # structure and finding their neighbourhoods
    ligand_neighbourhoods: dict[LigandID, LigandNeighbourhood] = {}
    for dataset in system_data.dataset:
        dataset_ligand_neighbourhoods: dict[
            LigandID, LigandNeighbourhood
        ] = get_dataset_neighbourhoods(dataset)
        ligand_neighbourhoods.update(dataset_ligand_neighbourhoods)

    return ligand_neighbourhoods
