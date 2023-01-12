from xchemalign.data import CanonicalSite, LigandNeighbourhood


def match_atom(canonical_site_atom, ligand_neighbourhood_atom) -> bool:
    # TODO: Makre sure that these have sensible types
    return False


def match_atoms(
    site_1_atoms,
    site_2_atoms,
    min_alignable_atoms: int = 5,
) -> bool:
    # Check if there is an alignable number of atoms shared between the
    num_alignable_atoms: int = 0
    for canonical_site_atom in site_1_atoms:
        for ligand_neighbourhood_atom in site_2_atoms:
            if match_atom(canonical_site_atom, ligand_neighbourhood_atom):
                num_alignable_atoms += 1

    if num_alignable_atoms > min_alignable_atoms:
        return True
    else:
        return False


def match_neighbourhood_to_site(
    canonical_site: CanonicalSite,
    ligand_neighbourhood: LigandNeighbourhood,
    min_alignable_atoms: int = 5,
) -> bool:
    # Check if there is an alignable number of atoms shared between the
    return match_atoms(
        canonical_site.atoms, ligand_neighbourhood.atoms, min_alignable_atoms
    )


def match_neighbourhood_to_sites(
    canonical_sites: dict[int, CanonicalSite],
    ligand_neighbourhood: LigandNeighbourhood,
) -> int | None:
    for canonical_site_id, canonical_site in canonical_sites.items():
        match: bool = match_neighbourhood_to_site(
            canonical_site, ligand_neighbourhood
        )
        if match:
            return match

    return None
