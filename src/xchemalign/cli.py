from pathlib import Path

import fire

from xchemalign.align_xmaps import _align_xmaps

# from xchemalign.get_system_sites import get_system_sites
from xchemalign.build_alignment_graph import build_alignment_graph

# from xchemalign import constants
from xchemalign.data import (
    LigandNeighbourhoods,
    Sites,
    SystemData,
    Transforms,
    XtalForms,
)
from xchemalign.generate_aligned_structures import _align_structures
from xchemalign.generate_sites_from_components import (
    _generate_sites_from_components,
)
from xchemalign.make_data_json import make_data_json_from_pandda_dir


def _update_sites(g, neighbourhoods, sites):
    ...


def _suggest_merges(sites: Sites):
    ...


def read_graph(pathpath: Path):
    ...


def read_sites(path: Path):
    ...


def read_neighbourhoods(path: Path):
    ...


def read_transforms(path: Path):
    ...


def read_structures(system_data):
    ...


def read_xtalforms(path: Path):
    ...


def align_structures():
    ...


# def _align_xmaps():
#     ...


def read_system_data(path: Path):
    ...


def report_site_update(sites, updated_sites):
    ...


def report_merges(suggested_merges):
    ...


def report_alignments():
    ...


class CLI:
    def merge_clusters(
        self, cluster_1: int, cluster_2: int, sites_path: str = "."
    ):
        # sites = read_sites(sites_path)
        ...

    def suggest_merges(self):
        ...

    def update_sites(self, path: str = "."):
        _path: Path = Path(path)

        # Read input data
        g = read_graph(_path)
        neighbourhoods = read_neighbourhoods(_path)
        sites = read_sites(_path)

        # Update sites
        updated_sites = _update_sites(g, neighbourhoods, sites)

        # Check for possible merges and report
        suggested_merges = _suggest_merges(sites)

        # Report
        report_site_update(sites, updated_sites)
        report_merges(suggested_merges)

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

    def change_site_reference(self):
        ...

    def align_structures(self):
        ...

    def align_xmaps(self):
        ...

    def generate_sites_from_components(self, source_dir: str):
        _source_dir: Path = Path(source_dir)

        _generate_sites_from_components(_source_dir)

    def align_all(self, source_dir: str, output_dir: str):
        _source_dir: Path = Path(source_dir)
        _output_dir: Path = Path(output_dir)

        g = read_graph(_source_dir)
        transforms: Transforms = read_transforms(_source_dir)
        neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
        xtalforms: XtalForms = read_xtalforms(_source_dir)
        sites: Sites = read_sites(_source_dir)
        system_data: SystemData = read_system_data(_source_dir)

        # get Structures
        structures = read_structures(system_data)

        # Align structures
        _align_structures(
            structures,
            sites,
            transforms,
            neighbourhoods,
            xtalforms,
            g,
            _output_dir,
        )

        # Align xmaps
        _align_xmaps(
            system_data, sites, neighbourhoods, transforms, _output_dir
        )

        # Report
        report_alignments()


if __name__ == "__main__":
    fire.Fire(CLI)
