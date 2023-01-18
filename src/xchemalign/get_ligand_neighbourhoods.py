import gemmi

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
) -> tuple[list[gemmi.NeighborSearch.Mark], list[gemmi.NeighborSearch.Mark]]:
    # Check each mark for its image and partition them on this
    model_atoms: list[gemmi.NeighborSearch.Mark] = []
    artefact_atoms: list[gemmi.NeighborSearch.Mark] = []
    for pos, mark in residue_neighbours.items():
        # Image 0 is the identity i.e. part of the normal model
        if mark.image_idx == 0:
            model_atoms.append(mark)
        else:
            artefact_atoms.append(mark)

    return model_atoms, artefact_atoms


def get_ligand_neighbourhood(
    structure: Structure,
    ns: gemmi.NeighborSearch,
    fragment: gemmi.Residue,
    min_dist: float = 0.01,
    max_dist: float = 3.0,
) -> LigandNeighbourhood:
    # For each atom, get the neighbouring atoms, and filter them on their
    # real space position
    residue_neighbours: dict[
        tuple[float, float, float], gemmi.NeighborSearch.Mark
    ] = {}
    for atom in fragment:
        atom_neighbours: list[gemmi.NeighborSearch.Mark] = ns.find_neighbors(
            atom,
            min_dist=min_dist,
            max_dist=max_dist,
        )
        for neighbour in atom_neighbours:
            residue_neighbours[
                round(neighbour.x, 2),
                round(neighbour.y, 2),
                round(neighbour.z, 2),
            ] = neighbour

    # Seperate out model and artefact atoms
    _model_atoms, _artefact_atoms = get_model_and_artefact_atoms(
        residue_neighbours
    )

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
    ligand_neighbourhood: LigandNeighbourhood = LigandNeighbourhood(
        model_atoms=model_atoms,
        artefact_atoms=artefact_atoms,
    )

    return ligand_neighbourhood


def get_dataset_neighbourhoods(
    dataset: Dataset, max_radius: float = 5.0
) -> dict[LigandID, LigandNeighbourhood]:
    # Load the structure
    structure: Structure = gemmi.read_structure(dataset.pdb)

    # Get the bound fragments
    fragments: dict[LigandID, gemmi.Residue] = get_structure_fragments(
        dataset, structure
    )

    # Construct the neighbourhood search
    ns: gemmi.NeighborSearch = gemmi.NeighborSearch(
        structure[0], structure.cell, max_radius
    )

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
