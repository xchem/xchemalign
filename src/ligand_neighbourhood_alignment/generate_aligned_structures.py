# import os
from pathlib import Path

import gemmi
import networkx as nx
from loguru import logger

from ligand_neighbourhood_alignment.data import (  # Transform,; AlignableSite,; XtalForms,
    AssignedXtalForms,
    CanonicalSites,
    ConformerSites,
    LigandNeighbourhoods,
    Output,
    SiteTransforms,
    Transforms,
    XtalForms,
    gemmi_to_transform,
    transform_to_gemmi,
)


def superpose_structure(transform, structure):
    new_structure = structure.clone()

    for model in new_structure:
        for chain in model:
            span = chain.whole()
            span.transform_pos_and_adp(transform)

    return new_structure


# def generate_aligned_structures(
#     output_dir, ligand_neighbourhoods, system_data,
# sites: list[AlignableSite]
# ):
#     # Get structures
#     structures = {}
#     for dataset in system_data.dataset:
#         structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
#         structures[dataset.dtag] = structure

#     #
#     for site in sites:
#         # logger.debug(f"Site id is: {site.id}")
#         site_dir = output_dir / f"{site.id}"
#         if not site_dir.exists():
#             os.mkdir(site_dir)
#         ligand_ids = site.ligand_ids
#         neighbourhoods = [ligand_neighbourhoods[n] for n in ligand_ids]
#         neigbhourhood_structures = [
#             structures[ligand_id.dtag] for ligand_id in ligand_ids
#         ]
#         transforms = [
#             get_transforms(neighbourhoods[0], neighbourhood)
#             for neighbourhood in neighbourhoods
#         ]
#         superposed_structures = [
#             superpose_structure(transform, structure)
#             for transform, structure in zip(
#                 transforms, neigbhourhood_structures
#             )
#         ]
#         for ligand_id, structure in zip(ligand_ids, superposed_structures):
#             structure.write_pdb(
#                 str(site_dir / f"{ligand_id.dtag}_{ligand_id.id}.pdb")
#             )


# def generate_aligned_structures_connected_components(
#     output_dir,
#     ligand_neighbourhoods,
#     system_data,
#     sites: Sites,
#     g,
# ):
#     # Get structures
#     structures = get_structures(system_data)

#     # Iterate sites
#     for site in sites:
#         logger.debug(f"Site id is: {site.id}")
#         site_dir = output_dir / f"{site.id}"
#         if not site_dir.exists():
#             os.mkdir(site_dir)
#         ligand_ids = site.ligand_ids

#         # Select the alignment reference ligand_id
#         reference_ligand_id = site.reference

#         # For each other ligand
#         for moving_ligand_id in ligand_ids[:-1]:
#             # Get the shortest alignment path to the reference
#             shortest_path = nx.shortest_path(
#                 g, moving_ligand_id, reference_ligand_id
#             )
#             logger.debug(f"Shortest path: {shortest_path}")

#             # Initial structure
#             structure = structures[moving_ligand_id.dtag].clone()

#             # Walk the path, iteratively applying transforms
#             previous_ligand_id = moving_ligand_id
#             for next_ligand_id in shortest_path:
#                 # Get the transform from previous frame to new one
#                 transform, alignment_ids = get_transforms(
#                     ligand_neighbourhoods[next_ligand_id],
#                     ligand_neighbourhoods[previous_ligand_id],
#                 )
#                 logger.debug(
#                     [f"{lid.residue}/{lid.atom}" for lid in alignment_ids]
#                 )

#                 # Apply the translation to the new frame
#                 structure = superpose_structure(transform, structure)

#             # Write the fully aligned structure
#             out_path = (
#                 site_dir / f"{moving_ligand_id.dtag}_{
# moving_ligand_id.id}.pdb"
#             )
#             structure.write_pdb(str(out_path))

#         # neighbourhoods = [ligand_neighbourhoods[n] for n in ligand_ids]
#         # neigbhourhood_structures = [
#         #     structures[ligand_id.dtag] for ligand_id in ligand_ids
#         # ]
#         # transforms = [
#         #     get_transforms(neighbourhoods[0], neighbourhood)
#         #     for neighbourhood in neighbourhoods
#         # ]
#         # superposed_structures = [
#         #     superpose_structure(transform, structure)
#         #     for transform, structure in zip(
#         #         transforms, neigbhourhood_structures
#         #     )
#         # ]
#         # for ligand_id, structure in zip(ligand_ids, superposed_structures):
#         #     structure.write_pdb(
#         #         str(site_dir / f"{ligand_id.dtag}_{ligand_id.id}.pdb")
#         #     )


def expand_structure(_structure, xtalforms: AssignedXtalForms, moving_ligand_id):
    # TODO: Make this work
    return _structure
    ...
    # for xtalform in xtalforms.xtalforms:
    #     if moving_ligand_id in xtalform.members:
    #         # Add the relvant images
    #         if xtalform.images:
    #             ...


# def _align_structures(
#     structures,
#     sites: Sites,
#     transforms: Transforms,
#     neighbourhoods: LigandNeighbourhoods,
#     xtalforms: XtalForms,
#     g,
#     _output_dir: Path,
# ):
#     # Iterate sites
#     for site in sites.sites:
#         logger.debug(f"Site id is: {site.id}")
#         site_dir = _output_dir / f"{site.id}"
#         if not site_dir.exists():
#             os.mkdir(site_dir)
#         ligand_ids = site.members

#         # Select the alignment reference ligand_id
#         reference_ligand_id = ligand_ids[0]

#         # For each other ligand
#         for moving_ligand_id in ligand_ids[:-1]:
#             # Get the shortest alignment path to the reference
#             shortest_path = nx.shortest_path(
#                 g, moving_ligand_id, reference_ligand_id
#             )
#             logger.debug(f"Shortest path: {shortest_path}")

#             # Initial structure
#             _structure = structures[moving_ligand_id.dtag].clone()

#             # Expand structure
#             structure = expand_structure(
#                 _structure, xtalforms, moving_ligand_id
#             )

#             # Walk the path, iteratively applying transforms
#             previous_ligand_id = moving_ligand_id
#             for next_ligand_id in shortest_path:
#                 # Get the transform from previous frame to new one
#                 transform = transforms.get_transform(
#                     (previous_ligand_id, next_ligand_id)
#                 )

#                 # logger.debug(
#                 #     [f"{lid.residue}/{lid.atom}" for lid in alignment_ids]
#                 # )

#                 # Apply the translation to the new frame
#                 structure = superpose_structure(transform, structure)

#             # Write the fully aligned structure
#             out_path = (
#
# site_dir / f"{moving_ligand_id.dtag}_{moving_ligand_id.id}.pdb"
#             )
#             structure.write_pdb(str(out_path))


def align_structure(
        _structure,
        moving_ligand_id,
        reference_ligand_id,
        g,
        transforms,
        site_transforms: SiteTransforms,
        canonical_site_id,
        conformer_site_id,
        out_path,
):
    shortest_path = nx.shortest_path(g, moving_ligand_id, reference_ligand_id)
    logger.debug(f"Shortest path: {shortest_path}")

    previous_ligand_id = moving_ligand_id
    running_transform = gemmi.Transform()
    for next_ligand_id in shortest_path:
        # Get the transform from previous frame to new one
        # Transform is 2 onto 1
        if next_ligand_id != previous_ligand_id:
            transform = transforms.get_transform(
                (
                    next_ligand_id,
                    previous_ligand_id,
                ),
            )
            running_transform = transform.combine(running_transform)

        # Apply the translation to the new frame
        previous_ligand_id = next_ligand_id

    # Subsite alignment transform
    subsite_transform = transform_to_gemmi(
        site_transforms.get_conformer_site_transform(canonical_site_id, conformer_site_id)
    )
    subsite_transform.combine(running_transform)

    # Site alignment transform
    site_transform = transform_to_gemmi(site_transforms.get_canonical_site_transform(canonical_site_id))
    site_transform.combine(running_transform)

    logger.debug(f"Transform from native frame to reference frame is: {gemmi_to_transform(running_transform)}")

    _structure = superpose_structure(running_transform, _structure)

    # Write the fully aligned structure
    _structure.write_pdb(str(out_path))


from ligand_neighbourhood_alignment import dt


def _align_structure(
        _structure,
        moving_ligand_id: tuple[str, str, str],
        reference_ligand_id: tuple[str, str, str],
        g,
        neighbourhood_transforms: dict[tuple[tuple[str, str, str], tuple[str, str, str]], dt.Transform],
        conformer_site_transforms: dict[tuple[str, str], dt.Transform],
        # canonical_site_transforms: dict[str, dt.Transform],
        canonical_site_id: str,
        conformer_site_id: str,
        out_path: Path,
):
    shortest_path: list[tuple[str, str, str]] = nx.shortest_path(g, moving_ligand_id, reference_ligand_id)
    logger.debug(f"Shortest path: {shortest_path}")

    previous_ligand_id = moving_ligand_id
    running_transform = gemmi.Transform()
    for next_ligand_id in shortest_path:
        # Get the transform from previous frame to new one
        # Transform is 2 onto 1
        if next_ligand_id != previous_ligand_id:
            transform = transform_to_gemmi(neighbourhood_transforms[(next_ligand_id, previous_ligand_id,)])
            running_transform = transform.combine(running_transform)

        # Apply the translation to the new frame
        previous_ligand_id = next_ligand_id

    # Subsite alignment transform
    confomer_site_transform = transform_to_gemmi(
        conformer_site_transforms[(canonical_site_id, conformer_site_id)]
    )
    running_transform = confomer_site_transform.combine(running_transform)

    # Site alignment transform
    # canonical_site_transform = transform_to_gemmi(canonical_site_transforms[canonical_site_id])
    # running_transform = canonical_site_transform.combine(running_transform)

    logger.debug(f"Transform from native frame to reference frame is: {gemmi_to_transform(running_transform)}")

    _structure = superpose_structure(running_transform, _structure)

    # Write the fully aligned structure
    _structure.write_pdb(str(out_path))


def _align_reference_structure(
        _structure,
        dtag: str,
        # moving_ligand_id: tuple[str,str,str],
        # reference_ligand_id: tuple[str,str,str],
        # g,
        # neighbourhood_transforms: dict[tuple[tuple[str,str,str], tuple[str,str,str]], dt.Transform],
        # conformer_site_transforms: dict[tuple[str, str], dt.Transform],
        reference_structure_transforms: dict[tuple[str, str], dt.Transform],
        # canonical_site_transforms: dict[str, dt.Transform],
        canonical_site_id: str,
        # conformer_site_id: str,
        out_path: Path,
):
    running_transform = transform_to_gemmi(reference_structure_transforms[(dtag, canonical_site_id)])

    # Site alignment transform
    # canonical_site_transform = transform_to_gemmi(canonical_site_transforms[canonical_site_id])
    # canonical_site_transform.combine(running_transform)

    logger.debug(f"Transform from native frame to reference frame is: {gemmi_to_transform(running_transform)}")

    _structure = superpose_structure(running_transform, _structure)

    # Write the fully aligned structure
    _structure.write_pdb(str(out_path))


# def align_artefacts():
#     # Transform the artefacts
#     xtalform = xtalforms.get_xtalform(
#         assigned_xtalforms.get_xtalform_id(moving_ligand_id.dtag)
#     )
#     assembly = generate_assembly(xtalform, _structure)
#     neighbourhood = neighbourhoods.get_neighbourhood(moving_ligand_id)
#     neighbourhood = neighbourhood
#     remove_non_contact_chains(assembly, neighbourhood)

#     _structure = superpose_structure(running_transform, assembly)

#     # Write the fully aligned structure
#     out_name = "{dtag}_{residue}_assembly.pdb"
#     out_path = subsite_dir / out_name.format(
#         dtag=moving_ligand_id.dtag,
#         residue=moving_ligand_id.residue,
#     )
#     _structure.write_pdb(str(out_path))


def _align_structures_from_sites(
        structures,
        canonical_sites: CanonicalSites,
        conformer_sites: ConformerSites,
        transforms: Transforms,
        neighbourhoods: LigandNeighbourhoods,
        xtalforms: XtalForms,
        assigned_xtalforms: AssignedXtalForms,
        g,
        site_transforms: SiteTransforms,
        # _output_dir: Path,
        output: Output,
):
    # asd = _output_dir / "aligned"
    # if not asd.exists():
    #     os.mkdir(asd)
    # Iterate sites
    for canonical_site_id, canonical_site in canonical_sites.iter():
        logger.debug(f"Canonical Site id is: {canonical_site_id}")
        # site_dir = asd / f"{site.id}"
        # if not site_dir.exists():
        #     os.mkdir(site_dir)

        for conformer_site_id, conformer_site in conformer_sites.iter():
            if conformer_site_id not in canonical_site.subsite_ids:
                continue
            # subsite_id: int = subsite.id
            # subsite_dir = site_dir / f"{canonical_site_id}"
            logger.debug(f"Conformer Site id is: {canonical_site_id}")

            # if not subsite_dir.exists():
            #     os.mkdir(subsite_dir)

            ligand_ids = conformer_site.members
            logger.debug(f"Conformer Site members: {len(conformer_site.members)}")

            # Select the alignment reference ligand_id
            # TODO: Track reference properly
            reference_ligand_id = conformer_site.reference_ligand_id

            # For each other ligand
            for moving_ligand_id in ligand_ids:
                # if moving_ligand_id.dtag != "Mpro-J0055":
                #     continue

                logger.info(f"Alligning ligand: {moving_ligand_id}")
                # Get the shortest alignment path to the reference
                # Initial structure
                _structure = structures[moving_ligand_id.dtag].clone()

                # Get output path
                aod = Path(output.source_dir)
                output_path = (
                        aod
                        / output.dataset_output[moving_ligand_id.dtag][moving_ligand_id.chain][
                            moving_ligand_id.residue
                        ].aligned_structures[canonical_site_id]
                )
                # Align the ligand
                align_structure(
                    _structure,
                    moving_ligand_id,
                    reference_ligand_id,
                    g,
                    transforms,
                    site_transforms,
                    canonical_site_id,
                    conformer_site_id,
                    output_path,
                )

                # Align the artefacts
                # align_artefacts(
                #     _structure, ligand_output.aligned_structures[cs_id]
                # )

                # Expand structure
                # structure = expand_structure(
                #     _structure, assigned_xtalforms, moving_ligand_id
                # )

                # Walk the path, iteratively applying transforms
