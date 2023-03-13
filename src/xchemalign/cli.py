import os

# import os
import subprocess

# import sys
from pathlib import Path

import fire
import numpy as np
import pandas as pd
from loguru import logger
from rich import print

from xchemalign import constants
from xchemalign.align_xmaps import _align_xmaps

# from xchemalign.get_system_sites import get_system_sites
from xchemalign.build_alignment_graph import build_alignment_graph
from xchemalign.data import (  # save_xtalforms,
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
from xchemalign.generate_aligned_structures import _align_structures_from_sites
from xchemalign.generate_sites_from_components import _generate_sites_from_components
from xchemalign.make_data_json import (
    get_ligand_binding_events_from_panddas,
    get_ligand_binding_events_from_structure,
    make_data_json_from_pandda_dir,
)

# logger.add(sys.stdout, colorize=True, format="{time} {level} \n{message}")


# def _update_sites(source_dir: Path):

#     ...


def _suggest_merges(sites: CanonicalSites):
    ...


def align_structures():
    ...


# def _align_xmaps():
#     ...


def report_site_update(sites, updated_sites):
    ...


def report_merges(suggested_merges):
    ...


def report_alignments():
    ...


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


# def _ligand_binding_events_from_pdb(pdb: Path, xmap: Path):
#     st = gemmi.read_structure(str(pdb))
#     for model in st:
#         for chain in model:
#             for residue in chain:
#                 if residue.name == "LIG":
#                     lbe = LigandBindingEvent(id=)


def _add_model_building_dir(_source_dir: Path, _data_source_dir: Path):
    system_data = read_system_data(_source_dir)

    if not _source_dir.exists():
        raise Exception(f"No such dir: {_source_dir}")

    datasource_paths = [_datasource.path for _datasource in system_data.datasources]
    logger.info(f"Datasources are: {datasource_paths}")

    # for model_dir in _source_dir.glob("*"):
    #     dtag = model_dir.name
    #     mtz = model_dir / constants.MODEL_DIR_MTZ
    #     pdb = model_dir / constants.MODEL_DIR_PDB
    #     xmap = model_dir / constants.MODEL_DIR_XMAP
    #     # ligand_binding_events = _ligand_binding_events_from_pdb(pdb, xmap)
    #     # dataset = Dataset(dtag=dtag, pdb=pdb,
    # ligand_binding_events=ligand_binding_events)
    #     initial_dataset = InitialDataset(dtag=dtag, mtz, pdb, xmap)

    datasource = Datasource(path=str(_data_source_dir), datasource_type="model_building")

    if not system_data.datasources:
        logger.info("No Datasources: Creating new list!")
        system_data.datasources = [
            datasource,
        ]
    else:
        new_datasources = [
            _datasource for _datasource in system_data.datasources if _datasource.path != str(_data_source_dir)
        ] + [
            datasource,
        ]
        system_data.datasources = new_datasources

    save_data(system_data, _source_dir)
    logger.info(f"Added dir {_data_source_dir} to datasources")
    datasource_paths = [_datasource.path for _datasource in system_data.datasources]
    logger.info(f"Datasources are now: {datasource_paths}")


def _add_manual_dir(_source_dir: Path, _data_source_dir: Path):
    system_data = read_system_data(_source_dir)

    if not _source_dir.exists():
        raise Exception(f"No such dir: {_source_dir}")

    datasource = Datasource(path=str(_data_source_dir), datasource_type="manual")

    if not system_data.datasources:
        system_data.datasources = [
            datasource,
        ]
    else:
        new_datasources = [
            _datasource for _datasource in system_data.datasources if _datasource.path != str(_data_source_dir)
        ] + [
            datasource,
        ]
        system_data.datasources = new_datasources

    save_data(system_data, _source_dir)
    logger.info(f"Added dir {_data_source_dir} to datasources")
    datasource_paths = [_datasource.path for _datasource in system_data.datasources]
    logger.info(f"Datasources are: {datasource_paths}")


def _add_pandda(_source_dir: Path, _pandda_dir: Path):
    system_data = read_system_data(_source_dir)

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

        save_data(system_data, _source_dir)
    else:
        raise Exception(f"No event table at: {event_table_path}")

    logger.info(f"Added PanDDA {_pandda_dir} to panddas")
    pandda_paths = [_pandda.path for _pandda in system_data.panddas]
    logger.debug(f"PanDDAs are: {pandda_paths}")


def _parse_data_sources(_source_dir: Path):
    system_data = read_system_data(_source_dir)

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

                # mtz = model_dir / constants.MODEL_DIR_MTZ
                pdb = model_dir / constants.MODEL_DIR_PDB
                xmap = model_dir / constants.MODEL_DIR_XMAP
                mtz = model_dir / constants.MODEL_DIR_MTZ
                if not pdb.exists():
                    continue
                # if not xmap.exists():
                #     continue
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
                # dataset_ids.append(dataset_id)
                datasets[dataset_id] = dataset
                logger.debug(f"Added dataset: {dataset_id}")
        else:
            raise Exception(f"Source type {datasource.datasource_type} unknown!")

    system_data.dataset_ids = list(datasets.keys())
    system_data.datasets = list(datasets.values())

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


def _assign_xtalforms(
    _source_dir: Path,
    assemblies: Assemblies,
    xtalforms: XtalForms,
    system_data: SystemData,
    # dataset_xtalforms: DatasetXtalforms,
):
    # system_data: SystemData = read_system_data(_source_dir)
    # ligand_neighbourhoods = read_neighbourhoods(_source_dir)
    # assemblies = read_assemblies(_source_dir)
    # xtalforms = read_xtalforms(_source_dir)
    structures = read_structures(system_data)

    dataset_ids = []
    xtalform_ids = []
    for dataset_id, dataset in system_data.iter():
        # if dataset_id.dtag in dataset_xtalforms.dataset_xtalforms:
        #     closest_xtalform_id = dataset_xtalforms.dataset_xtalforms[
        #         dataset_id.dtag
        #     ]

        # else:
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

    save_assigned_xtalforms(_source_dir, assigned_xtalforms)


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

        # xtalforms = XtalForms(
        #     xtalform_ids=[_xtalform.id for _xtalform in options.xtalforms],
        #     xtalforms=[_xtalform for _xtalform in options.xtalforms],
        # )
        # Copy the assembly json into the source directory (checking validity)
        assemblies = Assemblies.read(Path(options.assemblies_json))
        assemblies.save(Path(options.source_dir) / constants.ASSEMBLIES_FILE_NAME)

        # Copy the xtalform json into the source directory (checking validity)
        # xtalforms = read_xtalforms(Path(options.xtalforms_json))
        # save_xtalforms(Path(options.source_dir), xtalforms)
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
        # sites = read_sites(Path(options.source_dir))
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
                aligned_structures={},
                aligned_artefacts={},
                aligned_xmaps={},
            )

            # Add output for each canonical site that the ligand is aligned to
            for site_id, site in canonical_sites.iter():
                if ligand_id not in site.members:
                    continue

            chain_output.aligned_ligands[residue].aligned_structures[
                site_id
            ] = constants.ALIGNED_STRUCTURE_TEMPLATE.format(dtag=dtag, chain=chain, residue=residue, site=site_id)

            chain_output.aligned_ligands[residue].aligned_artefacts[
                site_id
            ] = constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                dtag=dtag, chain=chain, residue=residue, site=site_id
            )

            chain_output.aligned_ligands[residue].aligned_xmaps[site_id] = constants.ALIGNED_XMAP_TEMPLATE.format(
                dtag=dtag, chain=chain, residue=residue, site=site_id
            )

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
        # dataset_xtalforms = read_dataset_xtalforms(_source_dir)

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
            system_data=str(_source_dir / constants.DATA_JSON_PATH),
            xtalforms=str(_source_dir / constants.XTALFORMS_FILE_NAME),
            assigned_xtalforms=str(_source_dir / constants.ASSIGNED_XTALFORMS_FILE_NAME),
            neighbourhoods=str(_source_dir / constants.NEIGHBOURHOODS_FILE_NAME),
            graph=str(_source_dir / constants.ALIGNABILITY_GRAPH_FILE_NAME),
            transforms=str(_source_dir / constants.TRANSFORMS_FILE_NAME),
            sites=str(_source_dir / constants.SITES_FILE_NAME),
            site_transforms=str(_source_dir / constants.SITES_TRANSFORMS_FILE_NAME),
            dataset_output={},
        )
        save_output(output, _source_dir)

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

    def open_site(self, source_dir: str, site_id: int):
        _source_dir = Path(source_dir)
        script_path = _source_dir / "coot_script.py"
        script = ""
        script += 'if __name__ == "__main__": \n'
        script += '\tset_nomenclature_errors_on_read("ignore")\n'
        script += "\tset_recentre_on_read_pdb(0) \n"

        str_dir = _source_dir / constants.ALIGNED_STRUCTURES_DIR

        for site_dir in str_dir.glob("*"):
            if site_dir.name != f"{site_id}":
                continue

            for subsite_dir in site_dir.glob("*"):
                for pdb in subsite_dir.glob("*"):
                    script += f'\tp = read_pdb("{pdb}")\n '
                    script += cas_ligands()

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

    def update(self, system_data_dir: str, source_dir: str):

        # _source_dir: Path = Path(source_dir)

        self.build_system_data(system_data_dir, source_dir)
        self.build_graph(source_dir)
        # self.update_sites(source_dir)
        self.align_structures(source_dir)
        self.align_xmaps(source_dir)

    # def update_sites(self, path: str = "."):
    #     _path: Path = Path(path)

    #     # Read input data
    #     g = read_graph(_path)
    #     neighbourhoods = read_neighbourhoods(_path)
    #     sites = read_sites(_path)

    #     # Update sites
    #     updated_sites = _update_sites(g, neighbourhoods, sites)

    #     _update_sites(_path)
    #     # # Check for possible merges and report
    #     # suggested_merges = _suggest_merges(sites)

    #     # # Report
    #     # report_site_update(sites, updated_sites)
    #     # report_merges(suggested_merges)

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
            _source_dir,
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
        system_data: SystemData = read_system_data(_source_dir)
        site_transforms = read_site_transforms(_source_dir)

        # get Structures
        structures = read_structures(system_data)

        _align_xmaps(
            system_data,
            structures,
            sites,
            neighbourhoods,
            g,
            transforms,
            site_transforms,
            _source_dir,
        )

    def generate_sites_from_components(self, source_dir: str):
        _source_dir: Path = Path(source_dir)

        _generate_sites_from_components(_source_dir)


if __name__ == "__main__":
    fire.Fire(CLI)
