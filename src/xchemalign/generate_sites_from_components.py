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
    SiteTransforms,
    SubSite,
    Transform,
    read_graph,
    read_neighbourhoods,
    read_system_data,
)
from xchemalign.save_sites import save_sites
from xchemalign.structures import get_structures, get_transform_from_residues


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
            subsite_ids=[subsites[j].id for j in component],
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


def get_subsite_transforms(sites: Sites, structures):

    transforms = {}
    for site_id, site in zip(sites.site_ids, sites.sites):
        rss = site.subsites[0].members[0].dtag
        rs = site.residues
        srs = structures[rss]

        for ssid, ss in zip(site.subsite_ids, site.subsites):
            ssr = ss.members[0].dtag
            ssrs = structures[ssr]
            transform = get_transform_from_residues(rs, srs, ssrs)
            transforms[(site_id, rss, ssid)] = Transform(
                vec=transform.vec.tolist(), mat=transform.mat.tolist()
            )

    return transforms


def get_site_transforms(sites: Sites, structures):
    transforms = {}
    rs = sites.sites[0]
    rsid = sites.site_ids[0]

    rss = structures[rs.subsites[0].members[0].dtag]
    ref_site_all_ress = [
        ResidueID(chain=chain.name, residue=res.seqid.num)
        for model in rss
        for chain in model
        for res in chain
    ]

    for site_id, site in zip(sites.site_ids, sites.sites):
        srs = site.members[0].dtag
        site_structure = structures[srs]

        transform = get_transform_from_residues(
            ref_site_all_ress, rss, site_structure
        )
        transforms[rsid, site_id] = Transform(
            vec=transform.vec.tolist(), mat=transform.mat.tolist()
        )

    return transforms


def _generate_sites_from_components(_source_dir: Path):

    logger.info(f"Source dir: {_source_dir}")
    g = read_graph(_source_dir)
    neighbourhoods: LigandNeighbourhoods = read_neighbourhoods(_source_dir)
    logger.info(f"Number of neighbourhoods: {len(neighbourhoods.ligand_ids)}")

    system_data = read_system_data(_source_dir)

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

    # Get the subsite transforms
    logger.info("Getting transfroms between subsites...")
    structures = get_structures(system_data)
    subsite_transforms = get_subsite_transforms(sites, structures)

    # Get the site transforms
    logger.info("Getting transforms between sites...")
    site_transforms = get_site_transforms(sites, structures)
    site_transforms = SiteTransforms(
        site_transform_ids=[key for key in site_transforms.keys()],
        site_transforms=[tr for tr in site_transforms.values()],
        subsite_transform_ids=[key for key in subsite_transforms.keys()],
        subsite_transforms=[tr for tr in subsite_transforms.values],
    )

    return sites
