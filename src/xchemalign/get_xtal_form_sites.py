# import loguru
from xchemalign.data import (
    CanonicalSite,
    LigandID,
    LigandNeighbourhood,
    SystemSites,
    XtalFormSite,
)
from xchemalign.matching import match_atoms


def get_xtal_form_sites(
    initial_system_sites: SystemSites | None,
    ligand_neighbourhoods: dict[LigandID, LigandNeighbourhood],
    canonical_sites: dict[int, CanonicalSite],
) -> dict[int, XtalFormSite]:
    # Get the current xtalform sites
    xtal_form_sites: dict[int, XtalFormSite] = {}
    if initial_system_sites:
        xtal_form_sites = initial_system_sites.xtal_form_site

    # Get the number of xtalform sites
    if len(xtal_form_sites) == 0:
        xtal_form_sites_num = 0
    else:
        xtal_form_sites_num = max(canonical_sites)

    # Iterate through canonical sites, partitioning their members according to
    # whether their artefact atoms can be aligned
    for canonical_site_id, canonical_site in canonical_sites.items():
        # Determine whether each member aligns to an existing xtalform site
        for xtal_form_id, xtal_form_site in xtal_form_sites.items():
            for ligand_id in canonical_site.members:
                ligand_neighbourhood: LigandNeighbourhood = (
                    ligand_neighbourhoods[ligand_id]
                )
                match: bool = match_atoms(
                    xtal_form_site.artefact_atoms,
                    ligand_neighbourhood.artefact_atoms,
                )
                if match:
                    if ligand_id not in xtal_form_site.members:
                        xtal_form_site.members.append(ligand_id)

                # Otherwise create a new xtalform site and add the ligand to it
                else:
                    xtal_form_sites_num = xtal_form_sites_num + 1
                    xtal_form_sites[xtal_form_sites_num] = XtalFormSite(
                        canon_site_id=canonical_site_id,
                        xtal_form_id=0,
                        code="",
                        refpdb="",
                        atoms=ligand_neighbourhood.atoms,
                        artefact_atoms=ligand_neighbourhood.artefact_atoms,
                        members=[
                            ligand_id,
                        ],
                    )

    return xtal_form_sites
