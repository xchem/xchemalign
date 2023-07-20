import json
import re
from pathlib import Path

import pandas as pd
from loguru import logger

from ligand_neighbourhood_alignment import constants
from ligand_neighbourhood_alignment.make_data_json import (
    get_ligand_binding_events_from_panddas,
    get_ligand_binding_events_from_structure,
)


class LigandNeighbourhoodOutput:
    def __init__(self,
                 aligned_structures: dict[int, str],
                 aligned_artefacts: dict[int, str],

                 aligned_xmaps: dict[int, str],
                 aligned_event_maps: dict[int, str],
                 ):
        self.aligned_structures = aligned_structures
        self.aligned_artefacts: dict[int, str] = aligned_artefacts
        self.aligned_xmaps: dict[int, str] = aligned_xmaps
        self.aligned_event_maps: dict[int, str] = aligned_event_maps

    @staticmethod
    def from_dict(dic):
        return LigandNeighbourhoodOutput(
            aligned_structures=dic["aligned_structures"],
            aligned_artefacts=dic["aligned_artefacts"],
            aligned_xmaps=dic["aligned_xmaps"],
            aligned_event_maps=dic["aligned_event_maps"],
        )


class FSModel:
    def __init__(
            self,
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
            alignments,
            reference_alignments
    ):
        self.assemblies = assemblies
        self.xtalforms = xtalforms
        self.dataset_assignments = dataset_assignments
        self.ligand_neighbourhoods = ligand_neighbourhoods
        self.alignability_graph = alignability_graph
        self.ligand_neighbourhood_transforms = ligand_neighbourhood_transforms
        self.conformer_sites = conformer_sites
        self.conformer_site_transforms = conformer_site_transforms
        self.canonical_sites = canonical_sites
        self.canonical_site_transforms = canonical_site_transforms
        self.xtalform_sites = xtalform_sites
        self.alignments = alignments
        self.reference_alignments = reference_alignments

    @staticmethod
    def from_dir(source_dir: str):
        source_dir = Path(source_dir)

        fs_model = source_dir / constants.FS_MODEL_FILE_NAME
        if fs_model.exists():
            return FSModel.from_file(fs_model)

        else:
            assemblies = source_dir / constants.ASSEMBLIES_FILE_NAME
            xtalforms = source_dir / constants.ASSEMBLIES_FILE_NAME
            dataset_assignments = source_dir / constants.ASSIGNED_XTALFORMS_FILE_NAME
            ligand_neighbourhoods = source_dir / constants.NEIGHBOURHOODS_FILE_NAME
            alignability_graph = source_dir / constants.ALIGNABILITY_GRAPH_FILE_NAME
            ligand_neighbourhood_transforms = source_dir / constants.TRANSFORMS_FILE_NAME
            conformer_sites = source_dir / constants.CONFORMER_SITE_FILE
            conformer_site_transforms = source_dir / constants.CONFORMER_SITES_TRANSFORMS_FILE_NAME
            canonical_sites = source_dir / constants.CANONICAL_SITE_FILE
            canonical_site_trasnforms = source_dir / constants.CANONICAL_SITES_TRANSFORMS_FILE_NAME
            xtalform_sites = source_dir / constants.XTALFORM_SITE_FILE
            alignments = {}
            reference_alignments = {}

            return FSModel(
                assemblies,
                xtalforms,
                dataset_assignments,
                ligand_neighbourhoods,
                alignability_graph,
                ligand_neighbourhood_transforms,
                conformer_sites,
                conformer_site_transforms,
                canonical_sites,
                canonical_site_trasnforms,
                xtalform_sites,
                alignments,
                reference_alignments
            )

    @staticmethod
    def from_dict(dic):
        alignments = {}
        for dtag, dataset_alignments in alignments["alignments"].items():
            alignments[dtag] = {}
            for chain, chain_alignments in dataset_alignments.items():
                alignments[dtag][chain] = {}
                for ligand_neighbourhood, ligand_neighbourhood_alignments in chain_alignments.items():
                    alignments[dtag][chain][ligand_neighbourhood] = LigandNeighbourhoodOutput.from_dict(
                        ligand_neighbourhood_alignments)

        reference_alignments = {}
        for dtag, dataset_alignments in alignments["reference_alignments"].items():
            alignments[dtag] = {}
            for chain, chain_alignments in dataset_alignments.items():
                alignments[dtag][chain] = {}
                for ligand_neighbourhood, ligand_neighbourhood_alignments in chain_alignments.items():
                    alignments[dtag][chain][ligand_neighbourhood] = LigandNeighbourhoodOutput.from_dict(
                        ligand_neighbourhood_alignments)

        return FSModel(
            assemblies=Path(dic['assemblies']),
            xtalforms=Path(dic['xtalforms']),
            dataset_assignments=Path(dic['dataset_assignments']),
            ligand_neighbourhoods=Path(dic['ligand_neighbourhoods']),
            alignability_graph=Path(dic['alignability_graph']),
            ligand_neighbourhood_transforms=Path(dic['ligand_neighbourhood_transforms']),
            conformer_sites=Path(dic['conformer_sites']),
            conformer_site_transforms=Path(dic['conformer_site_transforms']),
            canonical_sites=Path(dic['canonical_sites']),
            canonical_site_transforms=Path(dic['canonical_site_transforms']),
            xtalform_sites=Path(dic['xtalform_sites']),
            alignments=alignments,
            reference_alignments=reference_alignments,
        )

    def to_dict(self, ):
        ...


class Datasource:
    path: str
    datasource_type: str


class PanDDA:
    path: str
    event_table_path: str


class Dataset:
    def __init__(self,
                 dtag,
                 pdb,
                 xmap,
                 mtz,
                 ligand_binding_events,
                 ):
        ...


class SourceDataModel:

    def __init__(self,
                 fs_model: FSModel,
                 datasources: list[Datasource],
                 panddas: list[PanDDA],
                 ):
        self.fs_model = fs_model
        self.datasources = datasources
        self.panddas = panddas

    @staticmethod
    def from_fs_model(
            fs_model: FSModel,
            datasources,
            datasource_types,
            panddas,
    ):
        ...

    def get_datasets(self):

        datasets = {}
        reference_datasets = {}
        new_datasets = {}

        # Get the pandda tables
        pandda_event_tables = {pandda.path: pd.read_csv(pandda.event_table_path) for pandda in self.panddas}

        # Get all the datasets attested in the data sources
        for datasource in self.datasources:
            logger.info(f"Parsing datasource: {datasource.path}")
            if datasource.datasource_type == "model_building":
                for model_dir in Path(datasource.path).glob("*"):
                    dtag = model_dir.name
                    if dtag in datasets:
                        st = f"Dataset ID {dtag} already found! Using new!"
                        logger.warning(st)
                        continue

                    pdb = model_dir / constants.MODEL_DIR_PDB
                    xmap = model_dir / constants.MODEL_DIR_XMAP
                    mtz = model_dir / constants.MODEL_DIR_MTZ
                    if not pdb.exists():
                        continue

                    ligand_binding_events = get_ligand_binding_events_from_panddas(
                        pandda_event_tables,
                        pdb,
                        dtag,
                    )
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
                    datasets[dtag] = dataset
                    logger.debug(f"Added dataset: {dtag}")

            elif datasource.datasource_type == "manual":
                for model_dir in Path(datasource.path).glob("*"):
                    dtag = model_dir.name

                    if dtag in datasets:
                        st = f"Dataset ID {dtag} already found! Using new!"
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
                    datasets[dtag] = dataset
                    reference_datasets[dtag] = dataset
                    logger.debug(f"Added dataset: {dtag}")
            else:
                raise Exception(f"Source type {datasource.datasource_type} unknown!")

        # Determine which of these are new using the fs_model output
        for dtag, dataset in datasets.items():
            if dtag not in self.fs_model.alignments:
                new_datasets[dtag] = dataset

        return datasets, reference_datasets, new_datasets

    @staticmethod
    def from_dict():
        ...

    def to_dict(self, path: Path):
        ...


class Generator:
    def __init__(
            self,
            reference_chain: str,
            chain: str,
            triplet: str
    ):
        self.reference_chain: str = reference_chain
        self.chain: str = chain
        self.triplet: str = triplet


class Assembly:
    def __init__(self,
                 reference: str,
                 generators: list[Generator]
                 ):
        self.reference = reference
        self.generators = generators

    @staticmethod
    def from_dict(dic):
        reference = dic['reference']
        biomol = dic['biomol']
        chains = dic['chains']

        # Split biomol on commas and strip whitespace
        biomol_matches = re.findall(
            '([A-Z]+)',
            biomol
        )

        # Split chains on commas that do not follow a number, x,y or z and strip whitespace
        chain_matches = re.findall(
            '([A-Z]+([(]+[^()]+[)]+)*)',
            chains
        )

        # Make generators
        generators = []
        for biomol_match, chain_match in zip(biomol_matches, chain_matches):
            if len(chain_match[1]) == 0:
                xyz = 'x,y,z'
            else:
                xyz = chain_match[1]
            generators.append(
                Generator(
                    biomol_match,
                    chain_match[0],
                    xyz
                )
            )

        return Assembly(
            reference,
            generators
        )


class XtalFormAssembly:
    def __init__(
            self,
            assembly: str,
            chains: list[str],
    ):
        self.assembly = assembly
        self.chains = chains


class XtalForm:
    def __init__(self,
                 reference: str,
                 assemblies: dict[str, XtalFormAssembly]
                 ):
        self.reference = reference
        self.assemblies = assemblies

    @staticmethod
    def from_dict(dic):
        reference = dic['reference']
        assemblies = dic['assemblies']

        _assemblies = {}
        for xtalform_assembly_id, xtalform_assembly_info in assemblies.items():
            assemblies[xtalform_assembly_id] = {}
            assembly = xtalform_assembly_info['assembly']
            chains = xtalform_assembly_info['chains']
            chains_matches = re.findall(
                '([A-Z]+)',
                chains
            )
            _assemblies[xtalform_assembly_id] = XtalFormAssembly(
                assembly,
                [x for x in chains_matches]
            )
        return XtalForm(reference, _assemblies)


class Transform:
    def __init__(self,
                 vec,
                 mat):
        self.vec: list[float] = vec
        self.mat: list[list[float]] = mat

    @staticmethod
    def from_dict(dic):
        return Transform(
            dic['vec'],
            dic['mat']
        )


class Atom:
    def __init__(
            self,
            element: str,
            x: float,
            y: float,
            z: float,
            image: Transform,
    ):
        self.element: str = element
        self.x: float = x
        self.y: float = y
        self.z: float = z
        self.image: Transform = image

    @staticmethod
    def from_dict(dic):
        return Atom(
            dic["element"],
            dic["x"],
            dic["y"],
            dic["z"],
            Transform.from_dict(dic["image"])
        )


class Neighbourhood:
    def __init__(self,
                 atoms,
                 artefact_atoms
                 ):
        self.atoms = atoms
        self.artefact_atoms = artefact_atoms

    @staticmethod
    def from_dict(dic):
        atoms = {}
        artefact_atoms = {}

        _atoms = dic['atoms']
        _artefact_atoms = dic['artefact_atoms']
        for atom_id, atom_info in _atoms.items():
            chain, residue, atom = atom_id.split("/")
            atoms[(chain, residue, atom)] = Atom.from_dict(atom_info)
            ...

        for atom_id, atom_info in _artefact_atoms.items():
            chain, residue, atom = atom_id.split("/")
            artefact_atoms[(chain, residue, atom)] = Atom.from_dict(atom_info)

        return Neighbourhood(
            atoms,
            artefact_atoms
        )


class AlignabilityGraph:
    ...


class ConformerSite:
    def __init__(
            self,
            residues: list[tuple[str, str]],
            members: list[tuple[str, str, str]],
            reference_ligand_id: tuple[str, str, str]
    ):
        self.residues: list[tuple[str, str]] = residues
        self.members: list[tuple[str, str, str]] = members
        self.reference_ligand_id: tuple[str, str, str] = reference_ligand_id

    @staticmethod
    def from_dict(dic):
        residues = []
        for res in dic['residues']:
            chain, residue = res.split("/")
            residues.append((chain, residue))
        members = []
        for member in dic['members']:
            dtag, chain, residue = member.split("/")
            members.append((dtag, chain, residue))
        ref_dtag, ref_chain, ref_residue = dic["reference_ligand_id"].split("/")
        return ConformerSite(
            residues,
            members,
            (ref_dtag, ref_chain, ref_residue)
        )


class CanonicalSite:
    def __init__(
            self,
            conformer_site_ids: list[str],
            residues: list[tuple[str, str]],
            reference_conformer_site_id: str,
            global_reference_dtag: str
    ):
        self.conformer_site_ids: list[str] = conformer_site_ids
        self.residues: list[tuple[str, str]] = residues
        self.reference_conformer_site_id: str = reference_conformer_site_id
        self.global_reference_dtag: str = global_reference_dtag

    @staticmethod
    def from_dict(dic):
        residues = []
        for res in dic['residues']:
            chain, residue = res.split("/")
            residues.append((chain, residue))

        return CanonicalSite(
            dic['conformer_site_ids'],
            residues,
            dic['reference_conformer_site_id'],
            dic['global_reference_dtag']

        )


class XtalFormSite:
    def __init__(self,
                 xtalform_id: str,
                 crystallographic_chain: str,
                 canonical_site_id: str,
                 members: list[tuple[str, str, str]]
                 ):
        self.xtalform_id: str = xtalform_id
        self.crystallographic_chain: str = crystallographic_chain
        self.canonical_site_id: str = canonical_site_id
        self.members: list[tuple[str, str, str]] = members

    @staticmethod
    def from_dict(dic):
        members = []
        for member in dic['members']:
            dtag, chain, residue = member.split("/")
            members.append((dtag, chain, residue))
        return XtalFormSite(
            dic['xtalform_id'],
            dic['crystallographic_chain'],
            dic['canonical_site_id'],
            members
        )
