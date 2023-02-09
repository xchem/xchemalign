import os
from pathlib import Path

import gemmi
import networkx as nx
from loguru import logger

from xchemalign.data import (  # Transform,; AlignableSite,; XtalForms,
    AssignedXtalForms,
    LigandNeighbourhoods,
    Sites,
    SiteTransforms,
    Transforms,
    XtalForms,
    transform_to_gemmi,
)
from xchemalign.structures import (
    generate_assembly,
)  # remove_non_contact_chains,


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


def expand_structure(
    _structure, xtalforms: AssignedXtalForms, moving_ligand_id
):
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


def _align_structures_from_sites(
    structures,
    sites: Sites,
    transforms: Transforms,
    neighbourhoods: LigandNeighbourhoods,
    xtalforms: XtalForms,
    assigned_xtalforms: AssignedXtalForms,
    g,
    site_transforms: SiteTransforms,
    _output_dir: Path,
):
    asd = _output_dir / "aligned_structures"
    if not asd.exists():
        os.mkdir(asd)
    # Iterate sites
    for site_id, site in zip(sites.site_ids, sites.sites):
        logger.debug(f"Site id is: {site.id}")
        site_dir = asd / f"{site.id}"
        if not site_dir.exists():
            os.mkdir(site_dir)

        for subsite in site.subsites:
            subsite_id: int = subsite.id
            subsite_dir = site_dir / f"{subsite_id}"
            logger.debug(f"SubSite id is: {subsite_id}")

            if not subsite_dir.exists():
                os.mkdir(subsite_dir)

            ligand_ids = subsite.members
            logger.debug(f"Subsite members: {len(subsite.members)}")

            # Select the alignment reference ligand_id
            # TODO: Track reference properly
            reference_ligand_id = subsite.reference_ligand_id

            # For each other ligand
            for moving_ligand_id in ligand_ids:
                logger.info(f"Alligning ligand: {moving_ligand_id}")
                # Get the shortest alignment path to the reference
                # Initial structure
                _structure = structures[moving_ligand_id.dtag].clone()

                # Expand structure
                # structure = expand_structure(
                #     _structure, assigned_xtalforms, moving_ligand_id
                # )

                # Walk the path, iteratively applying transforms
                shortest_path = nx.shortest_path(
                    g, moving_ligand_id, reference_ligand_id
                )
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
                        running_transform = transform.combine(
                            running_transform
                        )

                    # Apply the translation to the new frame
                    previous_ligand_id = next_ligand_id

                # Subsite alignment transform
                subsite_transform = transform_to_gemmi(
                    site_transforms.get_subsite_transform(site_id, subsite_id)
                )
                subsite_transform.combine(running_transform)

                # Site alignment transform
                site_transform = transform_to_gemmi(
                    site_transforms.get_site_transform(site_id)
                )
                site_transform.combine(running_transform)

                _structure = superpose_structure(running_transform, _structure)

                # Write the fully aligned structure
                out_path = (
                    subsite_dir
                    / f"{moving_ligand_id.dtag}_{moving_ligand_id.residue}.pdb"
                )
                _structure.write_pdb(str(out_path))

                # Transform the artefacts
                xtalform = xtalforms.get_xtalform(
                    assigned_xtalforms.get_xtalform_id(moving_ligand_id.dtag)
                )
                assembly = generate_assembly(xtalform, _structure)
                neighbourhood = neighbourhoods.get_neighbourhood(
                    moving_ligand_id
                )
                neighbourhood = neighbourhood
                # remove_non_contact_chains(assembly, neighbourhood)

                _structure = superpose_structure(running_transform, assembly)

                # Write the fully aligned structure
                out_name = "{dtag}_{residue}_assembly.pdb"
                out_path = subsite_dir / out_name.format(
                    dtag=moving_ligand_id.dtag,
                    residue=moving_ligand_id.residue,
                )
                _structure.write_pdb(str(out_path))
