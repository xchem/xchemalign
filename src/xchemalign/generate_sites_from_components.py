from pathlib import Path

import networkx as nx
import numpy as np
from loguru import logger

from xchemalign.data import (
    AtomID,
    LigandNeighbourhood,
    LigandNeighbourhoods,
    ResidueID,
    Site,
    Sites,
    SubSite,
    read_graph,
    read_neighbourhoods,
)
from xchemalign.save_sites import save_sites


def get_components(g):
    cliques = list(nx.connected_components(g))
    logger.debug(f"Cliques are: {cliques}")
    return cliques


def get_residues_from_neighbourhood(n: LigandNeighbourhood):
    rids = []
    for atom in n.atoms:
        aid: AtomID = atom.atom_id
        rid = ResidueID(chain=aid.chain, residue=aid.residue)
        rids.append(rid)

    return list(set(rids))


def get_subsites_from_components(
    components, neighbourhoods: LigandNeighbourhoods
):
    ss = []
    j = 0
    for component in components:
        rs = []
        for lid in component:
            n: LigandNeighbourhood = neighbourhoods.get_neighbourhood(lid)
            lrs: list[ResidueID] = get_residues_from_neighbourhood(n)
            rs += lrs
        s = SubSite(id=j, name="", residues=list(set(rs)), members=component)
        j += 1
        ss.append(s)

    return ss


def get_sites_from_subsites(
    subsites: list[SubSite], neighbourhoods: LigandNeighbourhoods
):
    g = nx.Graph()

    # Add the nodes
    for ss in subsites:
        g.add_node(ss.id)

    # Form the site overlap matrix
    arr = np.zeros((len(subsites), len(subsites)))
    for ss1 in subsites:
        for ss2 in subsites:
            if ss1.id == ss2.id:
                continue
            v = set(ss1.residues).intersection(set(ss2.residues))
            logger.debug(f"{ss1.id} {ss2.id} {len(v)}")
            if len(v) > 5:
                arr[ss1.id, ss2.id] = 1

    # Complete the graph
    for idx, conn in np.ndenumerate(arr):
        x, y = idx

        if x == y:
            continue
        if conn:
            g.add_edge(x, y)

    logger.debug(f"Number of sites sharing residues: {len(g.edges)}")

    # Get the connected components
    cc = get_components(g)

    # Form the sites
    sites = []
    j = 0
    for component in cc:
        s = Site(
            id=j,
            subsites=[subsites[j] for j in component],
            members=list(
                set(sum([subsites[j].members for j in component], start=[]))
            ),
            residues=list(
                set(sum([subsites[j].residues for j in component], start=[]))
            ),
        )
        j += 1
        sites.append(s)

    return sites


def _generate_sites_from_components(_source_dir: Path):

    logger.info(f"Source dir: {_source_dir}")
    g = read_graph(_source_dir)
    neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
    logger.info(f"Number of neighbourhoods: {len(neighbourhoods.ligand_ids)}")

    # Get the connected components
    logger.info("Getiting connected components...")
    connected_components = get_components(g)
    logger.info(f"Number of connected components: {len(connected_components)}")

    # Get the subsites from the connected components with overlap
    logger.info("Geting sites...")
    subsites: list[SubSite] = get_subsites_from_components(
        connected_components, neighbourhoods
    )
    logger.info(f"Number of subsites: {len(subsites)}")

    # Merge the connected components with shared residues into sites
    logger.info("Getting sites...")
    _sites: list[Site] = get_sites_from_subsites(subsites, neighbourhoods)
    logger.info(f"Number of sites: {len(_sites)}")

    sites: Sites = Sites(site_ids=[s.id for s in _sites], sites=_sites)

    save_sites(sites, _source_dir)

    return sites
