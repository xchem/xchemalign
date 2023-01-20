import os

import gemmi
import networkx as nx
from loguru import logger

from xchemalign.data import AlignableSite
from xchemalign.matching import match_atom


def get_transforms(reference_neighbourhood, neighbourhood):

    alignable_cas = []
    for (
        ligand_1_atom_id,
        ligand_1_atom,
    ) in reference_neighbourhood.atoms.items():
        for (
            ligand_2_atom_id,
            ligand_2_atom,
        ) in neighbourhood.atoms.items():
            if ligand_1_atom_id.atom == "CA":
                if match_atom(ligand_1_atom, ligand_2_atom):
                    alignable_cas.append(
                        (
                            gemmi.Position(
                                ligand_1_atom.x,
                                ligand_1_atom.y,
                                ligand_1_atom.z,
                            ),
                            gemmi.Position(
                                ligand_2_atom.x,
                                ligand_2_atom.y,
                                ligand_2_atom.z,
                            ),
                        )
                    )

    if len(alignable_cas) < 3:
        return gemmi.Transform()

    sup = gemmi.superpose_positions(
        [alignable_ca[0] for alignable_ca in alignable_cas],
        [alignable_ca[1] for alignable_ca in alignable_cas],
    )
    # logger.debug(f"Superposition: rmsd {sup.rmsd} n {len(alignable_cas)}")

    return sup.transform


def superpose_structure(transform, structure):
    new_structure = structure.clone()

    for model in new_structure:
        for chain in model:
            span = chain.whole()
            span.transform_pos_and_adp(transform)

    return new_structure


def generate_aligned_structures(
    output_dir, ligand_neighbourhoods, system_data, sites: list[AlignableSite]
):
    # Get structures
    structures = {}
    for dataset in system_data.dataset:
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        structures[dataset.dtag] = structure

    #
    for site in sites:
        # logger.debug(f"Site id is: {site.id}")
        site_dir = output_dir / f"{site.id}"
        if not site_dir.exists():
            os.mkdir(site_dir)
        ligand_ids = site.ligand_ids
        neighbourhoods = [ligand_neighbourhoods[n] for n in ligand_ids]
        neigbhourhood_structures = [
            structures[ligand_id.dtag] for ligand_id in ligand_ids
        ]
        transforms = [
            get_transforms(neighbourhoods[0], neighbourhood)
            for neighbourhood in neighbourhoods
        ]
        superposed_structures = [
            superpose_structure(transform, structure)
            for transform, structure in zip(
                transforms, neigbhourhood_structures
            )
        ]
        for ligand_id, structure in zip(ligand_ids, superposed_structures):
            structure.write_pdb(
                str(site_dir / f"{ligand_id.dtag}_{ligand_id.id}.pdb")
            )


def generate_aligned_structures_connected_components(
    output_dir,
    ligand_neighbourhoods,
    system_data,
    sites: list[AlignableSite],
    g,
):
    # Get structures
    structures = {}
    for dataset in system_data.dataset:
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        structures[dataset.dtag] = structure

    # Iterate sites
    for site in sites:
        logger.debug(f"Site id is: {site.id}")
        site_dir = output_dir / f"{site.id}"
        if not site_dir.exists():
            os.mkdir(site_dir)
        ligand_ids = site.ligand_ids

        # Select the alignment reference ligand_id
        reference_ligand_id = ligand_ids[0]

        # For each other ligand
        for moving_ligand_id in ligand_ids[:-1]:
            # Get the shortest alignment path to the reference
            shortest_path = nx.shortest_path(
                g, moving_ligand_id, reference_ligand_id
            )
            logger.debug(f"Shortest path: {shortest_path}")

            # Initial structure
            structure = structures[moving_ligand_id.dtag].clone()

            # Walk the path, iteratively applying transforms
            previous_ligand_id = moving_ligand_id
            for next_ligand_id in shortest_path:
                # Get the transform from previous frame to new one
                transform = get_transforms(
                    ligand_neighbourhoods[next_ligand_id],
                    ligand_neighbourhoods[previous_ligand_id],
                )

                # Apply the translation to the new frame
                structure = superpose_structure(transform, structure)

            # Write the fully aligned structure
            out_path = (
                site_dir / f"{moving_ligand_id.dtag}_{moving_ligand_id.id}.pdb"
            )
            structure.write_pdb(str(out_path))

        # neighbourhoods = [ligand_neighbourhoods[n] for n in ligand_ids]
        # neigbhourhood_structures = [
        #     structures[ligand_id.dtag] for ligand_id in ligand_ids
        # ]
        # transforms = [
        #     get_transforms(neighbourhoods[0], neighbourhood)
        #     for neighbourhood in neighbourhoods
        # ]
        # superposed_structures = [
        #     superpose_structure(transform, structure)
        #     for transform, structure in zip(
        #         transforms, neigbhourhood_structures
        #     )
        # ]
        # for ligand_id, structure in zip(ligand_ids, superposed_structures):
        #     structure.write_pdb(
        #         str(site_dir / f"{ligand_id.dtag}_{ligand_id.id}.pdb")
        #     )
