import json
import os

# import os
import subprocess

# import sys
from pathlib import Path

import fire
import numpy as np
import pandas as pd
import yaml
from loguru import logger
from rich import print

from ligand_neighbourhood_alignment import constants
from ligand_neighbourhood_alignment.align_xmaps import _align_xmaps

# from ligand_neighbourhood_alignment.get_system_sites import get_system_sites
from ligand_neighbourhood_alignment.build_alignment_graph import build_alignment_graph
from ligand_neighbourhood_alignment.data import (  # save_xtalforms,
    Assemblies,
    AssignedXtalForms,
    CanonicalSites,
    ChainOutput,
    ConformerSites,
    Dataset,
    DatasetID,
    DatasetOutput,
    Datasource,
    LigandID,
    LigandNeighbourhoods,
    LigandOutput,
    Options,
    Output,
    PanDDA,
    SiteTransforms,
    SystemData,
    Transforms,
    XtalForms,
    read_assigned_xtalforms,
    read_canonical_sites,
    read_graph,
    read_neighbourhoods,
    read_output,
    read_site_transforms,
    read_structures,
    read_system_data,
    read_transforms,
    read_xtalforms,
    save_assigned_xtalforms,
    save_canonical_sites,
    save_data,
    save_output,
)
from ligand_neighbourhood_alignment.generate_aligned_structures import _align_structures_from_sites
from ligand_neighbourhood_alignment.generate_sites_from_components import (  # get_xtalform_sites_from_canonical_sites,
    _generate_sites_from_components,
    get_components,
    get_conformer_sites_from_components,
    get_site_transforms,
    get_sites_from_conformer_sites,
    get_structures,
    get_subsite_transforms,
)
from ligand_neighbourhood_alignment.get_alignability import get_alignability
from ligand_neighbourhood_alignment.get_graph import get_graph
from ligand_neighbourhood_alignment.get_ligand_neighbourhoods import get_ligand_neighbourhoods
from ligand_neighbourhood_alignment.make_data_json import (
    get_ligand_binding_events_from_panddas,
    get_ligand_binding_events_from_structure,
    make_data_json_from_pandda_dir,
)


def cas_ligands():
    return "\tgraphics_to_ca_plus_ligands_sec_struct_representation(p) \n"


def _change_sites_reference(_source_dir: Path, site_id: int):
    sites: CanonicalSites = read_canonical_sites(_source_dir)
    sites.reference_site_id = site_id
    save_canonical_sites(sites, _source_dir)

    logger.info('Run "update" to generate xmaps and structures with new reference')


def _change_site_reference(_source_dir: Path, site_id: int, subsite_id: int):
    sites: CanonicalSites = read_canonical_sites(_source_dir)
    site = sites.get_site(site_id)
    site.reference_subsite_id = subsite_id
    new_reference_subsite = site.get_subsite(subsite_id)
    site.reference_subsite = new_reference_subsite
    site.reference_ligand_id = new_reference_subsite.reference_ligand_id
    save_canonical_sites(sites, _source_dir)

    logger.info('Run "update" to generate xmaps and structures with new reference')


def _change_subsite_reference(
        _source_dir: Path,
        site_id: int,
        subsite_id: int,
        dtag: int,
        chain: str,
        residue: int,
):
    sites: CanonicalSites = read_canonical_sites(_source_dir)
    site = sites.get_site(site_id)
    site.reference_subsite_id = subsite_id
    subsite = site.get_subsite(subsite_id)
    new_lid = LigandID(dtag=dtag, chain=chain, residue=residue)
    if new_lid not in subsite.members:
        raise Exception(f"LigandID {new_lid} not in {subsite.members}")
    subsite.reference_ligand_id = new_lid

    save_canonical_sites(sites, _source_dir)

    logger.info('Run "update" to generate xmaps and structures with new reference')


def _add_model_building_dir_to_system_data(system_data: SystemData, _data_source_dir: Path):
    datasource = Datasource(path=str(_data_source_dir), datasource_type="model_building")

    if not system_data.datasources:
        logger.info("No Datasources: Creating new list!")
        system_data.datasources = [
            datasource,
        ]
    else:
        new_datasources = [
                              _datasource for _datasource in system_data.datasources if
                              _datasource.path != str(_data_source_dir)
                          ] + [
                              datasource,
                          ]
        system_data.datasources = new_datasources

    return system_data


def _add_model_building_dir(_source_dir: Path, _data_source_dir: Path):
    if not _source_dir.exists():
        raise Exception(f"No such dir: {_source_dir}")
    system_data = read_system_data(_source_dir)

    datasource_paths = [_datasource.path for _datasource in system_data.datasources]
    logger.info(f"Datasources are: {datasource_paths}")

    system_data = _add_model_building_dir_to_system_data(system_data, _data_source_dir)

    save_data(system_data, _source_dir)
    logger.info(f"Added dir {_data_source_dir} to datasources")
    datasource_paths = [_datasource.path for _datasource in system_data.datasources]
    logger.info(f"Datasources are now: {datasource_paths}")


def _add_manual_dir_to_system_data(system_data: SystemData, _data_source_dir: Path):
    datasource = Datasource(path=str(_data_source_dir), datasource_type="manual")

    if not system_data.datasources:
        system_data.datasources = [
            datasource,
        ]
    else:
        new_datasources = [
                              _datasource for _datasource in system_data.datasources if
                              _datasource.path != str(_data_source_dir)
                          ] + [
                              datasource,
                          ]
        system_data.datasources = new_datasources

    return system_data


def _add_manual_dir(_source_dir: Path, _data_source_dir: Path):
    system_data = read_system_data(_source_dir)

    if not _source_dir.exists():
        raise Exception(f"No such dir: {_source_dir}")

    system_data = _add_manual_dir_to_system_data(system_data, _data_source_dir)

    save_data(system_data, _source_dir)
    logger.info(f"Added dir {_data_source_dir} to datasources")
    datasource_paths = [_datasource.path for _datasource in system_data.datasources]
    logger.info(f"Datasources are: {datasource_paths}")


def _add_pandda_to_system_data(system_data: SystemData, _pandda_dir: Path):
    analyses_dir: Path = _pandda_dir / constants.PANDDA_ANALYSES_DIR
    event_table_path: Path = analyses_dir / constants.PANDDA_EVENTS_INSPECT_TABLE_PATH

    if event_table_path.exists():

        pandda = PanDDA(path=str(_pandda_dir), event_table_path=str(event_table_path))

        if not system_data.panddas:
            system_data.panddas = [
                pandda,
            ]
        else:
            new_panddas = [_pandda for _pandda in system_data.panddas if _pandda.path != str(_pandda_dir)] + [
                pandda,
            ]
            system_data.panddas = new_panddas

    else:
        raise Exception(f"No event table at: {event_table_path}")

    return system_data


def _add_pandda(_source_dir: Path, _pandda_dir: Path):
    system_data = read_system_data(_source_dir)

    system_data = _add_pandda_to_system_data(system_data, _pandda_dir)

    save_data(system_data, _source_dir)

    logger.info(f"Added PanDDA {_pandda_dir} to panddas")
    pandda_paths = [_pandda.path for _pandda in system_data.panddas]
    logger.debug(f"PanDDAs are: {pandda_paths}")


def _add_data_to_system_data(system_data):
    # Get the PanDDA event tables
    pandda_event_tables = {pandda.path: pd.read_csv(pandda.event_table_path) for pandda in system_data.panddas}
    logger.info(f"Read {len(pandda_event_tables)} PanDDA event tables")

    # Get the
    datasets = {}
    for datasource in system_data.datasources:
        logger.info(f"Parsing datasource: {datasource.path}")
        if datasource.datasource_type == "model_building":
            for model_dir in Path(datasource.path).glob("*"):
                dtag = model_dir.name
                dataset_id = DatasetID(dtag=dtag)
                if dataset_id in datasets:
                    st = f"Dataset ID {dataset_id} already found! Using new!"
                    logger.warning(st)
                    continue

                pdb = model_dir / constants.MODEL_DIR_PDB
                xmap = model_dir / constants.MODEL_DIR_XMAP
                mtz = model_dir / constants.MODEL_DIR_MTZ
                if not pdb.exists():
                    continue

                ligand_binding_events = get_ligand_binding_events_from_panddas(pandda_event_tables, pdb, dtag)
                if len(ligand_binding_events.ligand_ids) == 0:
                    logger.warning(f"Dataset {dtag} has no ligand binding events!")
                    continue
                dataset = Dataset(
                    dtag=dtag,
                    pdb=str(pdb),
                    xmap=str(xmap),
                    mtz=str(mtz),
                    ligand_binding_events=ligand_binding_events,
                )
                # dataset_ids.append(dataset_id)
                datasets[dataset_id] = dataset
                logger.debug(f"Added dataset: {dataset_id}")

        elif datasource.datasource_type == "manual":
            for model_dir in Path(datasource.path).glob("*"):
                dtag = model_dir.name
                dataset_id = DatasetID(dtag=dtag)

                if dataset_id in datasets:
                    st = f"Dataset ID {dataset_id} already found! Using new!"
                    logger.warning(st)
                try:
                    pdb = next(model_dir.glob("*.pdb"))
                except Exception:
                    raise Exception(f"Could not find pdb in dir: {model_dir}")
                try:
                    xmap = next(model_dir.glob("*.ccp4"))
                except Exception as e:
                    print(e)
                    xmap = None
                    logger.warning("No xmap!")
                try:
                    mtz = next(model_dir.glob("*.mtz"))
                except Exception as e:
                    print(e)
                    mtz = None
                    logger.warning("No mtz!")

                ligand_binding_events = get_ligand_binding_events_from_structure(pdb, xmap, dtag)
                if len(ligand_binding_events.ligand_ids) == 0:
                    logger.warning(f"Dataset {dtag} has no ligand binding events!")
                    continue
                dataset = Dataset(
                    dtag=dtag,
                    pdb=str(pdb),
                    xmap=str(xmap),
                    mtz=str(mtz),
                    ligand_binding_events=ligand_binding_events,
                )
                datasets[dataset_id] = dataset
                logger.debug(f"Added dataset: {dataset_id}")
        else:
            raise Exception(f"Source type {datasource.datasource_type} unknown!")

    system_data.dataset_ids = list(datasets.keys())
    system_data.datasets = list(datasets.values())

    return system_data


def _parse_data_sources(_source_dir: Path):
    system_data = read_system_data(_source_dir)

    system_data = _add_data_to_system_data(system_data)

    save_data(system_data, _source_dir)

    logger.info(f"Found {len(system_data.dataset_ids)} datasets!")


def save_schema(model, path):
    with open(path / model.__name__, "w") as f:
        f.write(model.schema_json(indent=2))


def get_closest_xtalform(xtalforms: XtalForms, structures, dataset_id):
    structure = structures[dataset_id]
    structure_spacegroup = structure.spacegroup_hm
    structure_cell = structure.cell

    xtalform_deltas = {}

    for xtalform_id, xtalform in xtalforms.iter():
        ref_structure = structures[xtalform.reference]
        ref_spacegroup = ref_structure.spacegroup_hm
        ref_structure_cell = ref_structure.cell

        if ref_spacegroup != structure_spacegroup:
            continue

        deltas = np.array(
            [
                structure_cell.a / ref_structure_cell.a,
                structure_cell.b / ref_structure_cell.b,
                structure_cell.c / ref_structure_cell.c,
                structure_cell.alpha / ref_structure_cell.alpha,
                structure_cell.beta / ref_structure_cell.beta,
                structure_cell.gamma / ref_structure_cell.gamma,
            ]
        )
        xtalform_deltas[xtalform_id] = deltas

    if len(xtalform_deltas) == 0:
        return None, None

    closest_xtalform = min(
        xtalform_deltas,
        key=lambda _xtalform_id: np.sum(np.abs(xtalform_deltas[_xtalform_id] - 1)),
    )

    return closest_xtalform, xtalform_deltas[closest_xtalform]


from ligand_neighbourhood_alignment import dt
import gemmi


def _get_assigned_xtalforms(system_data, xtalforms):
    structures = read_structures(system_data)

    dataset_ids = []
    xtalform_ids = []
    for dataset_id, dataset in system_data.iter():

        closest_xtalform_id, deltas = get_closest_xtalform(xtalforms, structures, dataset_id)

        if (closest_xtalform_id is None) & (deltas is None):
            logger.info(f"No reference in same spacegroup for: {dataset_id}")
            logger.info(f"Structure path is: {dataset.pdb}")
            raise Exception()

        if np.any(deltas > 1.1) | np.any(deltas < 0.9):
            logger.info(f"No reference for dataset: {dataset_id}")
            logger.info(f"Deltas to closest unit cell are: {deltas}")
            logger.info(f"Structure path is: {dataset.pdb}")

            raise Exception()

        dataset_ids.append(dataset_id)
        xtalform_ids.append(closest_xtalform_id)

    assigned_xtalforms = AssignedXtalForms(dataset_ids=dataset_ids, xtalform_ids=xtalform_ids)

    return assigned_xtalforms


def _assign_xtalforms(
        _source_dir: Path,
        assemblies: Assemblies,
        xtalforms: XtalForms,
        system_data: SystemData,
):
    assigned_xtalforms = _get_assigned_xtalforms(system_data, xtalforms)

    save_assigned_xtalforms(_source_dir, assigned_xtalforms)

    return assigned_xtalforms


def _get_structures(datasets):
    structures = {}
    for dtag, dataset in datasets.items():
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        structures[dataset.dtag] = structure

    return structures


def _get_closest_xtalform(xtalforms: dict[str, dt.XtalForm], structure, structures):
    structure_spacegroup = structure.spacegroup_hm
    structure_cell = structure.cell

    xtalform_deltas = {}

    for xtalform_id, xtalform in xtalforms.items():
        ref_structure = structures[xtalform.reference]
        ref_spacegroup = ref_structure.spacegroup_hm
        ref_structure_cell = ref_structure.cell

        if ref_spacegroup != structure_spacegroup:
            continue

        deltas = np.array(
            [
                structure_cell.a / ref_structure_cell.a,
                structure_cell.b / ref_structure_cell.b,
                structure_cell.c / ref_structure_cell.c,
                structure_cell.alpha / ref_structure_cell.alpha,
                structure_cell.beta / ref_structure_cell.beta,
                structure_cell.gamma / ref_structure_cell.gamma,
            ]
        )
        xtalform_deltas[xtalform_id] = deltas

    if len(xtalform_deltas) == 0:
        return None, None

    closest_xtalform = min(
        xtalform_deltas,
        key=lambda _xtalform_id: np.sum(np.abs(xtalform_deltas[_xtalform_id] - 1)),
    )

    return closest_xtalform, xtalform_deltas[closest_xtalform]


def _assign_dataset(dataset, assemblies, xtalforms, structure, structures):
    closest_xtalform_id, deltas = _get_closest_xtalform(
        xtalforms,
        structure,
        structures,
    )

    if (closest_xtalform_id is None) & (deltas is None):
        logger.info(f"No reference in same spacegroup for")
        logger.info(f"Structure path is: {dataset.pdb}")
        raise Exception()

    if np.any(deltas > 1.1) | np.any(deltas < 0.9):
        logger.info(f"No reference for dataset")
        logger.info(f"Deltas to closest unit cell are: {deltas}")
        logger.info(f"Structure path is: {dataset.pdb}")

        raise Exception()

    return closest_xtalform_id


def _save_assignments(fs_model: dt.FSModel, dataset_assignments: dict[str, str]):
    with open(fs_model.dataset_assignments, 'w') as f:
        yaml.safe_dump(dataset_assignments, f)


def _generate_assembly(xtalform: dt.XtalForm, structure, assemblies: dict[str, dt.Assembly]):
    full_st = structure.clone()
    chains_to_delete = []
    for model in full_st:
        for chain in model:
            chains_to_delete.append((model.name, chain.name))

    for model_name, chain_name in chains_to_delete:
        del full_st[model_name][chain_name]

    for xtalform_assembly_id, xtalform_assembly in xtalform.assemblies.items():
        assembly = assemblies[xtalform_assembly.assembly]
        chains = xtalform_assembly.chains
        reference = assembly.reference

        for generator in assembly.generators:
            op = gemmi.Op(generator.triplet)
            chain_clone = structure[0][generator.chain].clone()
            for residue in chain_clone:
                for atom in residue:
                    atom_frac = structure.cell.fractionalize(atom.pos)
                    new_pos_frac = op.apply_to_xyz([atom_frac.x, atom_frac.y, atom_frac.z])
                    new_pos_orth = structure.cell.orthogonalize(gemmi.Fractional(*new_pos_frac))

                    atom.pos = gemmi.Position(*new_pos_orth)
            chain_clone.name = generator.reference_chain
            full_st[0].add_chain(chain_clone)

    num_chains = 0
    for model in full_st:
        for chain in model:
            num_chains += 1
    logger.debug(f"Generated {num_chains} assembly chains")

    return full_st


def _get_structure_fragments(dataset: dt.Dataset, structure):
    fragments: dict[tuple[str, str, str], gemmi.Residue] = {}
    # lig_number: int = 0
    for model in structure:
        for chain in model:
            for residue in chain.get_ligands():
                for lbe in dataset.ligand_binding_events.ligand_binding_events:
                    # if (
                    #     (residue.name == "LIG")
                    #     & (lbe.chain == chain.name)
                    #     & (lbe.residue == residue.seqid.num)
                    # ):
                    if (lbe.chain == chain.name) & (lbe.residue == residue.seqid.num):
                        ligand_id = (dataset.dtag, str(chain.name), str(lbe.residue),)
                        fragments[ligand_id] = residue
                    # lig_number = lig_number + 1

    return fragments


from ligand_neighbourhood_alignment.get_ligand_neighbourhoods import _get_ligand_neighbourhood


def _get_dataset_neighbourhoods(
        dataset: dt.Dataset, xtalform: dt.XtalForm, assemblies: dict[str, dt.Assembly], max_radius: float = 7.0
) -> dict[tuple[str, str, str], dt.Neighbourhood]:
    # Load the structure
    logger.debug(dataset.pdb)
    structure = gemmi.read_structure(dataset.pdb)
    logger.debug(f"{structure.cell}")

    # Get the rest of the assembly
    assembly = _generate_assembly(xtalform, structure, assemblies)

    # Get the bound fragments
    fragments: dict[tuple[str, str, str], gemmi.Residue] = _get_structure_fragments(dataset, assembly)
    logger.debug(f"Get {len(fragments)} fragment neighbourhoods")
    logger.debug(fragments)

    # Construct the neighbourhood search
    ns: gemmi.NeighborSearch = gemmi.NeighborSearch(assembly[0], assembly.cell, max_radius).populate()

    # For each bound fragment, identify the neighbourhood atoms and
    # partition them into model and artefact
    fragment_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood] = {}
    for ligand_id, fragment in fragments.items():
        fragment_neighbourhoods[ligand_id] = _get_ligand_neighbourhood(assembly, ns, fragment, max_dist=max_radius)

    return fragment_neighbourhoods


def _get_neighbourhoods(dataset: dt.Dataset, xtalform: dt.XtalForm, assemblies: dict[str, dt.Assembly]):
    dataset_ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood] = _get_dataset_neighbourhoods(
        dataset, xtalform, assemblies
    )
    return dataset_ligand_neighbourhoods


def _save_neighbourhoods(
        fs_model: dt.FSModel,
        ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood],
):
    with open(fs_model.ligand_neighbourhoods, 'w') as f:
        dic = {}
        for ligand_id, neighbourhood in ligand_neighbourhoods.items():
            dic["/".join(ligand_id)] = neighbourhood.to_dict()
        yaml.safe_dump(dic, f)


def _update(
        fs_model: dt.FSModel,
        datasets: dict[str, dt.Dataset],
        reference_datasets: dict[str, dt.Dataset],
        new_datasets: dict[str, dt.Dataset],
        assemblies: dict[str, dt.Assembly],
        xtalforms: dict[str, dt.XtalForm],
        dataset_assignments: dict[str, str],
        ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood],
        alignability_graph: dt.AlignabilityGraph,
        ligand_neighbourhood_transforms: dict[tuple[tuple[str, str, str], tuple[str, str, str]], dt.Transform],
        conformer_sites: dict[str, dt.ConformerSite],
        conformer_site_transforms: dict[tuple[str, str], dt.Transform],
        canonical_sites: dict[str, dt.CanonicalSite],
        canonical_site_transforms: dict[str, dt.Transform],
        xtalform_sites: dict[str, dt.XtalFormSite],
):
    # Get the structures
    structures: dict = _get_structures(datasets)

    # Assign datasets
    for dtag, dataset in new_datasets.items():
        dataset_assignments[dtag] = _assign_dataset(dataset, assemblies, xtalforms, structures[dtag], structures)
    _save_assignments(fs_model, dataset_assignments)

    # Get neighbourhoods
    for dtag, dataset in new_datasets.items():
        xtalform = xtalforms[dataset_assignments[dtag]]
        neighborhoods = _get_neighbourhoods(dataset, xtalform, assemblies)
        for lid, neighbourhood in neighborhoods.items():
            ligand_neighbourhoods[lid] = neighbourhood
    _save_neighbourhoods(fs_model, ligand_neighbourhoods)

    # Update graph
    for dtag, dataset in new_datasets.items():
        for lid, neighbourhood in ligand_neighbourhoods.items():
            if lid[0] in datasets:
                continue
            alignments, transforms = _get_alignments(neighbourhood, structures)
            for target_lid, transform in transforms.items():
                ligand_neighbourhood_transforms[(lid, target_lid)] = transform
            _update_graph(alignability_graph, alignments)
    _save_graph(fs_model, alignability_graph)

    # Update conformer sites
    connected_components = _get_connected_components(alignability_graph)
    for connected_component in connected_components:
        # Match new component to old ones by membership, and expand old ones if available otherwise create new one
        _update_conformer_sites(conformer_sites, connected_component)
    _save_conformer_sites(fs_model, conformer_sites)

    # Update canonical sites
    for conformer_site_id, conformer_site in conformer_sites.items():
        # If conformer site in a canonical site, replace with new data, otherwise
        # Check if residues match as usual, otherwise create a new canon site for it
        _update_canonical_sites(conformer_site)
    _save_canonical_sites(fs_model, canonical_sites)

    # Update crystalform sites
    for xtalform_site_id, xtalform_site in xtalform_sites.items():
        # If canonical site in a xtalform site, replace with new data, otherwise
        # Check if residues match as usual, otherwise create a new canon site for it
        _update_xtalform_sites(canonical_site, dataset_assignments)
    _save_xtalform_sites(fs_model, xtalform_sites)

    # Get conformer site transforms
    for canonical_site_id, canonical_site in canonical_sites.items():
        for conformer_site_id, conformer_site in canonical_site.conformer_sites.items():
            _update_conformer_site_transform(
                conformer_site_transforms,
                canonical_site,
                conformer_site,
            )
    _save_conformer_site_transforms(fs_model, conformer_site_transforms)

    # Get canonical site tranforms
    for canonical_site_id, canonical_site in canonical_sites.items():
        _update_canonical_site_transforms(
            canonical_site_transforms,
            canonical_site,
            canonical_sites,
        )
    _save_canonical_site_transforms(fs_model, canonical_site_transforms)

    # Update output: check if aligned data for each lid in canon site is already there and if not add it
    _update_fs_model(
        fs_model,
        canonical_sites,
    )
    _save_fs_model(fs_model)

    # Generate new aligned structures
    for canonical_site_id, canonical_site in canonical_sites.items():
        for conformer_site_id, conformer_site in canonical_site.conformer_sites.items():
            for lid in conformer_site.ligand_ids:
                _update_aligned_structures()

    # Generate alignments of references to each canonical site
    for canonical_site_id, canonical_site in canonical_sites.items():
        for dtag, reference_dataset in reference_datasets.items():
            _update_reference_alignments(

            )

    # Generate new aligned maps
    for canonical_site_id, canonical_site in canonical_sites.items():
        for conformer_site_id, conformer_site in canonical_site.conformer_sites.items():
            for lid in conformer_site.ligand_ids:
                _update_aligned_xmaps()


def _load_assemblies(assemblies_file, new_assemblies_yaml):
    assemblies = {}

    if assemblies_file.exists():

        with open(assemblies_file, 'r') as f:
            dic = yaml.safe_load(f)

        for assembly_id, assembly_info in dic.items():
            assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    # Load new info and update
    if new_assemblies_yaml.exists():
        with open(new_assemblies_yaml, 'r') as f:
            new_assemblies_dict = yaml.safe_load(f)
    else:
        new_assemblies_dict = {}

    for assembly_id, assembly_info in new_assemblies_dict.items():
        if assembly_id in assemblies:
            continue
        assemblies[assembly_id] = dt.Assembly.from_dict(assembly_info)

    return assemblies


def _load_xtalforms(xtalforms_file, new_xtalforms_yaml):
    xtalforms = {}

    if xtalforms_file.exists():

        with open(xtalforms_file, 'r') as f:
            dic = yaml.safe_load(f)

        for xtalform_id, xtalform_info in dic.items():
            xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    # Load new info and update
    if new_xtalforms_yaml.exists():
        with open(new_xtalforms_yaml, 'r') as f:
            new_xtalforms_dict = yaml.safe_load(f)
    else:
        new_xtalforms_dict = {}

    for xtalform_id, xtalform_info in new_xtalforms_dict.items():
        if xtalform_id in xtalforms:
            continue
        xtalforms[xtalform_id] = dt.XtalForm.from_dict(xtalform_info)

    return xtalforms


def _load_dataset_assignments(dataset_assignments_yaml):
    dataset_assignments = {}
    if dataset_assignments_yaml.exists():

        with open(dataset_assignments_yaml, 'r') as f:
            dic = yaml.safe_load(f)

        for dtag, assignment in dic.items():
            dataset_assignments[dtag] = assignment
    return dataset_assignments


def _load_ligand_neighbourhoods(ligand_neighbourhoods_yaml):
    ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood] = {}

    if ligand_neighbourhoods_yaml.exists():

        with open(ligand_neighbourhoods_yaml, 'r') as f:
            dic = yaml.safe_load(f)

        if dic:
            for ligand_id, neighbourhood_info in dic.items():
                dtag, chain, residue = ligand_id.split("/")
                neighbourhood = dt.Neighbourhood.from_dict(neighbourhood_info)
                ligand_neighbourhoods[(dtag, chain, residue)] = neighbourhood

    return ligand_neighbourhoods


def _load_alignability_graph(alignability_graph):
    if alignability_graph.exists():
        return nx.read_gml(
            str(alignability_graph),
            destringizer=lambda x: x.split("/"),
        )

    else:
        return None


def _load_ligand_neighbourhood_transforms(ligand_neighbourhood_transforms_yaml):
    ligand_neighbourhood_transforms = {}
    if ligand_neighbourhood_transforms_yaml.exists():

        with open(ligand_neighbourhood_transforms_yaml, 'r') as f:
            dic = yaml.safe_load(f)

        for ligand_transform_key, ligand_transform in dic.items():
            ligand_1_id, ligand_2_id = ligand_transform_key.split("~")
            dtag_1, chain_1, residue_1 = ligand_1_id.split("/")
            dtag_2, chain_2, residue_2 = ligand_1_id.split("/")
            ligand_neighbourhood_transforms[(
                (dtag_1, chain_1, residue_1),
                (dtag_2, chain_2, residue_2)
            )] = dt.Transform.from_dict(ligand_transform)

    return ligand_neighbourhood_transforms


def _load_conformer_sites(conformer_sites_yaml):
    conformer_sites = {}
    if conformer_sites_yaml.exists():
        with open(conformer_sites_yaml, 'r') as f:
            dic = yaml.safe_load(f)
        for conformer_site_id, conformer_site_info in dic.items():
            conformer_sites[conformer_site_id] = dt.ConformerSite.from_dict(conformer_site_info)

    return conformer_sites


def _load_conformer_site_transforms(conformer_site_transforms_yaml):
    conformer_site_transforms = {}
    if conformer_site_transforms_yaml.exists():
        with open(conformer_site_transforms_yaml, 'r') as f:
            dic = yaml.safe_load(f)

        for conformer_site_transform_id, conformer_site_transform_info in dic.items():
            conformer_site_1, conformer_site_2 = conformer_site_transform_id.split("/")

            conformer_site_transforms[(conformer_site_1, conformer_site_2)] = dt.Transform.from_dict(
                conformer_site_transform_info)

    return conformer_site_transforms


def _load_canonical_sites(canonical_sites_yaml):
    canonical_sites = {}
    if canonical_sites_yaml.exists():
        with open(canonical_sites_yaml, 'r') as f:
            dic = yaml.safe_load(f)

        for canonical_site_id, canonical_site_info in dic.items():
            canonical_sites[canonical_site_id] = dt.CanonicalSite.from_dict(canonical_site_info)

    return canonical_sites


def _load_canonical_site_transforms(canonical_site_transforms_yaml):
    canonical_site_transforms = {}
    if canonical_site_transforms_yaml.exists():
        with open(canonical_site_transforms_yaml, 'r') as f:
            dic = yaml.safe_load(f)

        for canonical_site_id, canonical_site_transform_info in dic.items():
            canonical_site_transforms[canonical_site_id] = dt.Transform.from_dict(canonical_site_transform_info)

    return canonical_site_transforms


def _load_xtalform_sites(xtalform_sites_yaml):
    xtalform_sites = {}
    if xtalform_sites_yaml.exists():
        with open(xtalform_sites_yaml, 'r') as f:
            dic = yaml.safe_load(f)

        for xtalform_site_id, xtalform_site_info in dic.items():
            xtalform_sites[xtalform_site_id] = dt.XtalFormSite.from_dict(xtalform_site_info)

    return xtalform_sites


class CLI:
    def schema(self, output_dir: str):
        _output_dir = Path(output_dir)

        if not _output_dir.exists():
            os.mkdir(_output_dir)

        save_schema(SystemData, _output_dir)
        save_schema(LigandNeighbourhoods, _output_dir)
        save_schema(CanonicalSites, _output_dir)
        save_schema(Transforms, _output_dir)
        save_schema(SiteTransforms, _output_dir)
        save_schema(SystemData, _output_dir)
        save_schema(SystemData, _output_dir)
        save_schema(SystemData, _output_dir)

    def update(self, options_json: str):

        options = Options.parse_file(options_json)
        # logger.info(f"Options json path is: {options}")
        logger.info(f"Output dir is: {options.source_dir}")

        fs_model = dt.FSModel.from_dir(options.source_dir)

        if not Path(options.source_dir).exists():
            os.mkdir(options.source_dir)

        source_data_model = dt.SourceDataModel.from_fs_model(
            fs_model,
            options.datasources,
            options.datasource_types,
            options.panddas
        )

        datasets, reference_datasets, new_datasets = source_data_model.get_datasets()

        # Get assemblies
        assemblies: dict[str, dt.Assembly] = _load_assemblies(fs_model.assemblies, Path(options.assemblies_json))
        for key, assembly in assemblies.items():
            print(assembly)
            for gen in assembly.generators:
                print([gen.chain, gen.reference_chain, gen.triplet])
        # Get xtalforms
        xtalforms: dict[str, dt.XtalForm] = _load_xtalforms(fs_model.xtalforms, Path(options.xtalforms_json))
        for key, xtalform in xtalforms.items():
            print(xtalform)
            for ass, xtalform_ass in xtalform.assemblies.items():
                for chn, trns in zip(xtalform_ass.chains, xtalform_ass.transforms):
                    print([chn, trns])

        # Get the dataset assignments
        dataset_assignments = _load_dataset_assignments(Path(fs_model.dataset_assignments))

        # Get Ligand neighbourhoods
        ligand_neighbourhoods: dict[tuple[str, str, str], dt.Neighbourhood] = _load_ligand_neighbourhoods(
            fs_model.ligand_neighbourhoods)

        # Get alignability graph
        alignability_graph = _load_alignability_graph(fs_model.alignability_graph)

        #
        ligand_neighbourhood_transforms: dict[
            tuple[tuple[str, str, str], tuple[str, str, str]], dt.Transform] = _load_ligand_neighbourhood_transforms(
            fs_model.ligand_neighbourhood_transforms)

        # Get conformer sites
        conformer_sites: dict[str, dt.ConformerSite] = _load_conformer_sites(fs_model.conformer_sites)

        #
        conformer_site_transforms: dict[tuple[str, str], dt.Transform] = _load_conformer_site_transforms(
            fs_model.conformer_site_transforms)

        # Get canonical sites
        canonical_sites: dict[str, dt.CanonicalSite] = _load_canonical_sites(fs_model.canonical_sites)

        #
        canonical_site_transforms: dict[str, dt.Transform] = _load_canonical_site_transforms(
            fs_model.conformer_site_transforms)

        # Get xtalform sites
        xtalform_sites: dict[str, dt.XtalFormSite] = _load_xtalform_sites(fs_model.xtalform_sites)

        # Run the update
        _update(
            fs_model,
            datasets,
            reference_datasets,
            new_datasets,
            assemblies,
            xtalforms,
            dataset_assignments,
            ligand_neighbourhoods,
            alignability_graph,
            ligand_neighbourhood_transforms,
            conformer_sites,
            conformer_site_transforms,
            canonical_sites,
            canonical_site_transforms,
            xtalform_sites,
        )

    def process_all(self, option_json: str):
        options = Options.parse_file(option_json)

        # Initialize the output directory and create empty
        # jsons in it
        system_data = self.init(options.source_dir)

        # Add the datasources in the options json and add them to
        # the datasource json
        for datasource_dir, datasource_type in zip(options.datasources, options.datasource_types):
            if datasource_type == "model_building":
                _add_model_building_dir_to_system_data(system_data, Path(datasource_dir))
            elif datasource_type == "manual":
                _add_manual_dir_to_system_data(system_data, Path(datasource_dir))

        # Add the PanDDAs in the options json and add them to the pandda json
        for pandda_dir in options.panddas:
            _add_pandda_to_system_data(system_data, Path(pandda_dir))

        # Copy the assembly json into the source directory (checking validity)
        # assemblies = Assemblies.read(Path(options.assemblies_json))
        # assemblies.save(Path(options.source_dir) / constants.ASSEMBLIES_FILE_NAME)

        # Copy the xtalform json into the source directory (checking validity)
        xtalforms = XtalForms.read(Path(options.xtalforms_json))
        # xtalforms.save(Path(options.source_dir) / constants.XTALFORMS_FILE_NAME)

        # Parse the data sources and PanDDAs, matching ligands up to events
        system_data = _add_data_to_system_data(system_data)

        # Assign each dataset to the clsoest xtalform and fail if this
        # is not possible
        assigned_xtalforms = _get_assigned_xtalforms(system_data, xtalforms)

        # Build the alignment graph
        ligand_neighbourhoods: LigandNeighbourhoods = get_ligand_neighbourhoods(
            system_data,
            xtalforms,
            assigned_xtalforms,
        )

        num_neighbourhoods = len(ligand_neighbourhoods.ligand_neighbourhoods)
        logger.info(f"Found {num_neighbourhoods} ligand neighbourhoods")

        # Get alignability
        logger.info("Getting alignbaility matrix...!")
        alignability_matrix, transforms = get_alignability(ligand_neighbourhoods, system_data)
        logger.info("Got alignability matrix!")

        # logger.debug(alignability_matrix)
        logger.debug("Alignability matrix shape: {alignability_matrix.shape}")

        # Generate the graph
        logger.info("Getting alignability graph...")
        g = get_graph(alignability_matrix, ligand_neighbourhoods)

        # Generate canonical, conformer and xtalform sites from the
        # alignment graph

        # Get the connected components
        logger.info("Getiting connected components...")
        connected_components = get_components(g)
        logger.info(f"Number of connected components: {len(connected_components)}")

        # Get the subsites from the connected components with overlap
        logger.info("Geting sites...")
        conformer_sites: ConformerSites = get_conformer_sites_from_components(
            connected_components, ligand_neighbourhoods
        )
        logger.info(f"Number of subsites: {len(conformer_sites.conformer_sites)}")

        # Merge the connected components with shared residues into sites
        logger.info("Getting sites...")
        _sites = get_sites_from_conformer_sites(conformer_sites, ligand_neighbourhoods)
        logger.info(f"Number of sites: {len(_sites)}")

        canonical_sites: CanonicalSites = CanonicalSites(
            site_ids=[s.id for s in _sites],
            sites=_sites,
            reference_site=_sites[0],
            reference_site_id=_sites[0].id,
        )

        # Get the xtalform sites
        # xtalform_sites = get_xtalform_sites_from_canonical_sites(
        #     canonical_sites,
        #     assigned_xtalforms,
        #     xtalforms,
        #     # assemblies,
        # )

        # Get the subsite transforms
        logger.info("Getting transfroms between subsites...")
        structures = get_structures(system_data)
        subsite_transforms = get_subsite_transforms(canonical_sites, structures)

        # Get the site transforms
        logger.info("Getting transforms between sites...")
        site_transforms = get_site_transforms(canonical_sites, structures)
        site_transforms = SiteTransforms(
            canonical_site_transform_ids=[key for key in site_transforms.keys()],
            canonical_site_transforms=[tr for tr in site_transforms.values()],
            conformer_site_transform_ids=[key for key in subsite_transforms.keys()],
            conformer_site_transforms=[tr for tr in subsite_transforms.values()],
        )
        # Fully specify the output now that the sites are known
        output = read_output(Path(options.source_dir))
        dataset_output_dict = {}
        for ligand_id in ligand_neighbourhoods.ligand_ids:
            dtag, chain, residue = (
                ligand_id.dtag,
                ligand_id.chain,
                ligand_id.residue,
            )

            if dtag not in dataset_output_dict:
                dataset_output = DatasetOutput(aligned_chain_output={})
                dataset_output_dict[dtag] = dataset_output
            else:
                dataset_output = dataset_output_dict[dtag]

            if chain not in dataset_output.aligned_chain_output:
                chain_output = ChainOutput(
                    aligned_ligands={},
                )
                dataset_output_dict[dtag].aligned_chain_output[chain] = chain_output
            else:
                chain_output = dataset_output_dict[dtag].aligned_chain_output[chain]

            chain_output.aligned_ligands[residue] = LigandOutput(
                aligned_structures={}, aligned_artefacts={}, aligned_xmaps={}, aligned_event_maps={}
            )

            # Add output for each canonical site that the ligand is aligned to
            for site_id, site in canonical_sites.iter():
                if ligand_id not in site.members:
                    continue

                chain_output.aligned_ligands[residue].aligned_structures[site_id] = (
                        output.aligned_dir
                        + "/"
                        + constants.ALIGNED_STRUCTURE_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

                chain_output.aligned_ligands[residue].aligned_artefacts[site_id] = (
                        output.aligned_dir
                        + "/"
                        + constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

                chain_output.aligned_ligands[residue].aligned_xmaps[site_id] = (
                        output.aligned_dir
                        + "/"
                        + constants.ALIGNED_XMAP_TEMPLATE.format(dtag=dtag, chain=chain, residue=residue, site=site_id)
                )

                chain_output.aligned_ligands[residue].aligned_event_maps[site_id] = (
                        output.aligned_dir
                        + "/"
                        + constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

        # Save the output file
        output.dataset_output = dataset_output_dict

        # Align structures to each canonical site
        _align_structures_from_sites(
            structures,
            canonical_sites,
            conformer_sites,
            transforms,
            ligand_neighbourhoods,
            xtalforms,
            assigned_xtalforms,
            g,
            site_transforms,
            output,
        )

        # Align xmaps to each canonical site
        _align_xmaps(
            system_data,
            structures,
            canonical_sites,
            conformer_sites,
            ligand_neighbourhoods,
            g,
            transforms,
            site_transforms,
            output,
        )

    def process(self, option_json: str):
        options = Options.parse_file(option_json)

        # Initialize the output directory and create empty
        # jsons in it
        self.init(options.source_dir)

        # Add the datasources in the options json and add them to
        # the datasource json
        for datasource_dir, datasource_type in zip(options.datasources, options.datasource_types):
            if datasource_type == "model_building":
                self.add_data_source(
                    options.source_dir,
                    datasource_dir,
                    source_type=datasource_type,
                )
            elif datasource_type == "manual":
                self.add_data_source(
                    options.source_dir,
                    datasource_dir,
                    source_type=datasource_type,
                )

        # Add the PanDDAs in the options json and add them to the pandda json
        for pandda_dir in options.panddas:
            self.add_pandda(options.source_dir, pandda_dir)

        # Copy the assembly json into the source directory (checking validity)
        assemblies = Assemblies.read(Path(options.assemblies_json))
        assemblies.save(Path(options.source_dir) / constants.ASSEMBLIES_FILE_NAME)

        # Copy the xtalform json into the source directory (checking validity)
        xtalforms = XtalForms.read(Path(options.xtalforms_json))
        xtalforms.save(Path(options.source_dir) / constants.XTALFORMS_FILE_NAME)

        # Parse the data sources and PanDDAs, matching ligands up to events
        self.parse_data_sources(options.source_dir)

        # Assign each dataset to the clsoest xtalform and fail if this
        # is not possible
        self.assign_xtalforms(options.source_dir)

        # Build the alignment graph
        self.build_graph(options.source_dir)

        # Generate canonical, conformer and xtalform sites from the
        # alignment graph
        self.generate_sites_from_components(options.source_dir)

        # Fully specify the output now that the sites are known
        neighbourhoods = read_neighbourhoods(Path(options.source_dir))
        canonical_sites = CanonicalSites.read(Path(options.source_dir) / constants.CANONICAL_SITE_FILE)
        output = read_output(Path(options.source_dir))
        dataset_output_dict = {}
        for ligand_id in neighbourhoods.ligand_ids:
            dtag, chain, residue = (
                ligand_id.dtag,
                ligand_id.chain,
                ligand_id.residue,
            )

            if dtag not in dataset_output_dict:
                dataset_output = DatasetOutput(aligned_chain_output={})
                dataset_output_dict[dtag] = dataset_output
            else:
                dataset_output = dataset_output_dict[dtag]

            if chain not in dataset_output.aligned_chain_output:
                chain_output = ChainOutput(
                    aligned_ligands={},
                )
                dataset_output_dict[dtag].aligned_chain_output[chain] = chain_output
            else:
                chain_output = dataset_output_dict[dtag].aligned_chain_output[chain]

            chain_output.aligned_ligands[residue] = LigandOutput(
                aligned_structures={}, aligned_artefacts={}, aligned_xmaps={}, aligned_event_maps={}
            )

            # Add output for each canonical site that the ligand is aligned to
            for site_id, site in canonical_sites.iter():
                if ligand_id not in site.members:
                    continue

                chain_output.aligned_ligands[residue].aligned_structures[site_id] = (
                        output.aligned_dir
                        + "/"
                        + constants.ALIGNED_STRUCTURE_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

                chain_output.aligned_ligands[residue].aligned_artefacts[site_id] = (
                        output.aligned_dir
                        + "/"
                        + constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

                chain_output.aligned_ligands[residue].aligned_xmaps[site_id] = (
                        output.aligned_dir
                        + "/"
                        + constants.ALIGNED_XMAP_TEMPLATE.format(dtag=dtag, chain=chain, residue=residue, site=site_id)
                )

                chain_output.aligned_ligands[residue].aligned_event_maps[site_id] = (
                        output.aligned_dir
                        + "/"
                        + constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

        # Save the output file
        output.dataset_output = dataset_output_dict
        save_output(output, Path(options.source_dir))

        # Align structures to each canonical site
        self.align_structures(options.source_dir)

        # Align xmaps to each canonical site
        self.align_xmaps(options.source_dir)

    def assign_xtalforms(self, source_dir: str):
        _source_dir = Path(source_dir)

        assemblies = Assemblies.read(_source_dir / constants.ASSEMBLIES_FILE_NAME)
        xtalforms = read_xtalforms(_source_dir)
        system_data = read_system_data(_source_dir)

        _assign_xtalforms(
            _source_dir,
            assemblies,
            xtalforms,
            system_data,  # dataset_xtalforms
        )

    def write_options_json(self, source_dir, options_json):
        _source_dir = Path(source_dir)
        _options_json = Path(options_json)

        system_data = read_system_data(_source_dir)
        options = Options(
            source_dir=source_dir,
            datasources=[ds.path for ds in system_data.datasources],
            panddas=[pandda.path for pandda in system_data.panddas],
        )
        with open(_options_json, "w") as f:
            f.write(options.json())

    def init(self, source_dir: str):
        _source_dir = Path(source_dir)

        if not _source_dir.exists():
            os.mkdir(_source_dir)

        system_data = SystemData(datasources=[], panddas=[], dataset_ids=[], datasets=[])

        save_data(system_data, _source_dir)

        output = Output(
            source_dir=str(_source_dir),
            system_data=str(constants.DATA_JSON_PATH),
            xtalforms=str(constants.XTALFORMS_FILE_NAME),
            assigned_xtalforms=str(constants.ASSIGNED_XTALFORMS_FILE_NAME),
            neighbourhoods=str(constants.NEIGHBOURHOODS_FILE_NAME),
            graph=str(constants.ALIGNABILITY_GRAPH_FILE_NAME),
            transforms=str(constants.TRANSFORMS_FILE_NAME),
            sites=str(constants.SITES_FILE_NAME),
            site_transforms=str(constants.SITES_TRANSFORMS_FILE_NAME),
            aligned_dir=str(constants.ALIGNED_STRUCTURES_DIR),
            dataset_output={},
        )
        if not (_source_dir / constants.ALIGNED_STRUCTURES_DIR).exists():
            os.mkdir(_source_dir / constants.ALIGNED_STRUCTURES_DIR)
        save_output(output, _source_dir)

        return system_data

    def add_data_source(
            self,
            source_dir: str,
            data_source_dir: str,
            source_type: str = "model_building",
    ):
        _source_dir = Path(source_dir)
        _data_source_dir = Path(data_source_dir)

        if source_type == "model_building":
            _add_model_building_dir(_source_dir, _data_source_dir)

        elif source_type == "manual":
            _add_manual_dir(_source_dir, _data_source_dir)

        else:
            raise Exception()

    def add_pandda(self, source_dir: str, pandda_dir: str):
        _source_dir = Path(source_dir)
        _pandda_dir = Path(pandda_dir)

        _add_pandda(_source_dir, _pandda_dir)

    def parse_data_sources(self, source_dir: str):
        _source_dir = Path(source_dir)

        _parse_data_sources(_source_dir)

    def open_site(
            self,
            option_json: str,
            site_id: int,
    ):
        options = Options.parse_file(option_json)
        output = Output.read(Path(options.source_dir) / constants.OUTPUT_JSON_PATH)

        _source_dir = Path(options.source_dir)
        script_path = _source_dir / "coot_script.py"
        script = ""
        script += 'if __name__ == "__main__": \n'
        script += '\tset_nomenclature_errors_on_read("ignore")\n'
        script += "\tset_recentre_on_read_pdb(0) \n"

        # str_dir = _source_dir / constants.ALIGNED_STRUCTURES_DIR

        # for site_dir in str_dir.glob("*"):
        #     if site_dir.name != f"{site_id}":
        #         continue

        #     for subsite_dir in site_dir.glob("*"):
        #         for pdb in subsite_dir.glob("*"):
        #             script += f'\tp = read_pdb("{pdb}")\n '
        #             script += cas_ligands()

        for dtag, dataset_output in output.dataset_output.items():
            for chain, chain_output in dataset_output.aligned_chain_output.items():
                for residue, residue_output in chain_output.aligned_ligands.items():
                    for _site_id, pdb in residue_output.aligned_structures.items():
                        logger.debug(_site_id)
                        if _site_id == site_id:
                            script += f'\tp = read_pdb("{options.source_dir}/{pdb}")\n '
                            script += cas_ligands()

        # for dataset in output.dataset_output

        logger.debug(script)

        with open(script_path, "w") as f:
            f.write(script)

        p = subprocess.Popen(f"coot --script {script_path}", shell=True)
        p.communicate()
        # os.remove(script_path)

    def merge_clusters(self, cluster_1: int, cluster_2: int, sites_path: str = "."):
        # sites = read_sites(sites_path)
        ...

    def suggest_merges(self):
        ...

    def pretty_print_dataset(self, source_dir: str):
        _source_dir = Path(source_dir)
        system_data = read_system_data(_source_dir)
        print(system_data)

    def align(
            self,
            system_data_dir: str,
            source_dir: str,
    ):
        self.build_system_data(system_data_dir, source_dir)
        self.build_graph(source_dir)
        self.generate_sites_from_components(source_dir)
        self.align_structures(source_dir)
        self.align_xmaps(source_dir)

    # def update(self, system_data_dir: str, source_dir: str):
    #
    #     # _source_dir: Path = Path(source_dir)
    #
    #     self.build_system_data(system_data_dir, source_dir)
    #     self.build_graph(source_dir)
    #     # self.update_sites(source_dir)
    #     self.align_structures(source_dir)
    #     self.align_xmaps(source_dir)

    def build_graph(
            self,
            source_dir: str,
    ):
        _source_dir: Path = Path(source_dir)

        build_alignment_graph(_source_dir)

    def build_system_data(self, system_data_dir: str, output_dir: str):
        print(output_dir)
        _system_data_dir: Path = Path(system_data_dir).resolve()
        _output_dir: Path = Path(output_dir).resolve()

        make_data_json_from_pandda_dir(_system_data_dir, _output_dir)

    def change_sites_reference(self, source_dir: str, site_id: int):
        _source_dir: Path = Path(source_dir)

        _change_sites_reference(_source_dir, site_id)

    def change_site_reference(self, source_dir: str, site_id: int, subsite_id: int):
        _source_dir: Path = Path(source_dir)

        _change_site_reference(_source_dir, site_id, subsite_id)

    def change_subsite_reference(
            self,
            source_dir: str,
            site_id: int,
            subsite_id: int,
            dtag: int,
            chain: str,
            residue: int,
    ):
        _source_dir: Path = Path(source_dir)

        _change_subsite_reference(_source_dir, site_id, subsite_id, dtag, chain, residue)

    def align_structures(self, source_dir: str):
        _source_dir: Path = Path(source_dir)
        # _output_dir: Path = Path(output_dir)

        g = read_graph(_source_dir)
        transforms: Transforms = read_transforms(_source_dir)
        neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
        # xtalforms: XtalForms = read_xtalforms(_source_dir)
        assigned_xtalforms = read_assigned_xtalforms(_source_dir)
        canonical_sites: CanonicalSites = CanonicalSites.read(_source_dir / constants.CANONICAL_SITE_FILE)
        conformer_sites: ConformerSites = ConformerSites.read(_source_dir / constants.CONFORMER_SITE_FILE)
        system_data: SystemData = read_system_data(_source_dir)
        site_transforms = read_site_transforms(_source_dir)
        xtalforms = read_xtalforms(_source_dir)
        output = Output.read(_source_dir / constants.OUTPUT_JSON_PATH)

        # get Structures
        structures = read_structures(system_data)

        # Align structures
        _align_structures_from_sites(
            structures,
            canonical_sites,
            conformer_sites,
            transforms,
            neighbourhoods,
            xtalforms,
            assigned_xtalforms,
            g,
            site_transforms,
            # _source_dir,
            output,
        )

    def align_xmaps(self, source_dir: str):
        _source_dir: Path = Path(source_dir)
        # _output_dir: Path = Path(output_dir)

        g = read_graph(_source_dir)
        transforms: Transforms = read_transforms(_source_dir)
        neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
        # xtalforms: XtalForms = read_xtalforms(_source_dir)
        # xtalforms = XtalForms(xtalforms=[])
        sites: CanonicalSites = CanonicalSites.read(_source_dir / constants.CANONICAL_SITE_FILE)
        conformer_sites = ConformerSites.read(_source_dir / constants.CONFORMER_SITE_FILE)
        system_data: SystemData = read_system_data(_source_dir)
        site_transforms = read_site_transforms(_source_dir)
        output = Output.read(_source_dir / constants.OUTPUT_JSON_PATH)

        # get Structures
        structures = read_structures(system_data)

        _align_xmaps(
            system_data,
            structures,
            sites,
            conformer_sites,
            neighbourhoods,
            g,
            transforms,
            site_transforms,
            output,
        )

    def generate_sites_from_components(self, source_dir: str):
        _source_dir: Path = Path(source_dir)

        _generate_sites_from_components(_source_dir)


if __name__ == "__main__":
    fire.Fire(CLI)
