from xchemalign.data import AlignableSite


def get_alignable_sites(connected_components):

    sites = []
    for site_ligand_ids in connected_components:
        sites.append(AlignableSite(id=0, name="", ligand_ids=site_ligand_ids))

    return sites
    ...
