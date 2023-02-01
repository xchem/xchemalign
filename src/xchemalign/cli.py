# import os
import subprocess
from pathlib import Path

import fire

from xchemalign import constants
from xchemalign.align_xmaps import _align_xmaps

# from xchemalign.get_system_sites import get_system_sites
from xchemalign.build_alignment_graph import build_alignment_graph
from xchemalign.data import (
    LigandNeighbourhoods,
    Sites,
    SystemData,
    Transforms,
    XtalForms,
    read_graph,
    read_neighbourhoods,
    read_site_transforms,
    read_sites,
    read_structures,
    read_system_data,
    read_transforms,
)
from xchemalign.generate_aligned_structures import _align_structures_from_sites
from xchemalign.generate_sites_from_components import (
    _generate_sites_from_components,
)
from xchemalign.make_data_json import make_data_json_from_pandda_dir


def _update_sites(g, neighbourhoods, sites):
    ...


def _suggest_merges(sites: Sites):
    ...


def read_xtalforms(path: Path):
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


class CLI:
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

    def merge_clusters(
        self, cluster_1: int, cluster_2: int, sites_path: str = "."
    ):
        # sites = read_sites(sites_path)
        ...

    def suggest_merges(self):
        ...

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

    def align_structures(self, source_dir: str):
        _source_dir: Path = Path(source_dir)
        # _output_dir: Path = Path(output_dir)

        g = read_graph(_source_dir)
        transforms: Transforms = read_transforms(_source_dir)
        neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
        # xtalforms: XtalForms = read_xtalforms(_source_dir)
        xtalforms = XtalForms(xtalforms=[])
        sites: Sites = read_sites(_source_dir)
        system_data: SystemData = read_system_data(_source_dir)
        site_transforms = read_site_transforms(_source_dir)

        # get Structures
        structures = read_structures(system_data)

        # Align structures
        _align_structures_from_sites(
            structures,
            sites,
            transforms,
            neighbourhoods,
            xtalforms,
            g,
            site_transforms,
            _source_dir,
        )

    def align_xmaps(self, source_dir: str):
        _source_dir: Path = Path(source_dir)
        # _output_dir: Path = Path(output_dir)

        g = read_graph(_source_dir)
        transforms: Transforms = read_transforms(_source_dir)
        neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
        # xtalforms: XtalForms = read_xtalforms(_source_dir)
        # xtalforms = XtalForms(xtalforms=[])
        sites: Sites = read_sites(_source_dir)
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

        # def align_all(self, source_dir: str, output_dir: str):
        #     _source_dir: Path = Path(source_dir)
        #     _output_dir: Path = Path(output_dir)

        #     g = read_graph(_source_dir)
        #     transforms: Transforms = read_transforms(_source_dir)
        #     neighbourhoods = read_neighbourhoods(_source_dir)
        #     xtalforms: XtalForms = read_xtalforms(_source_dir)
        #     sites: Sites = read_sites(_source_dir)
        #     system_data: SystemData = read_system_data(_source_dir)

        #     # get Structures
        #     structures = read_structures(system_data)

        # Align structures
        # _align_structures_from_sites(
        #     structures,
        #     sites,
        #     transforms,
        #     neighbourhoods,
        #     xtalforms,
        #     g,
        #     _output_dir,
        # )

        # # Align xmaps
        # _align_xmaps(
        #     system_data, sites, neighbourhoods, transforms, _output_dir
        # )

        # Report
        report_alignments()


if __name__ == "__main__":
    fire.Fire(CLI)
