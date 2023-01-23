import gemmi
import numpy as np
from loguru import logger

from xchemalign.data import (
    LigandNeighbourhood,
    LigandNeighbourhoods,
    SystemData,
)
from xchemalign.matching import match_atom


def match_cas(
    ligand_1_neighbourhood: LigandNeighbourhood,
    ligand_2_neighbourhood: LigandNeighbourhood,
    min_alignable_atoms: int = 5,
    max_alignable_rmsd: float = 2.0,
):
    alignable_cas = []
    for (
        ligand_1_atom_id,
        ligand_1_atom,
    ) in ligand_1_neighbourhood.atoms.items():
        for (
            ligand_2_atom_id,
            ligand_2_atom,
        ) in ligand_2_neighbourhood.atoms.items():
            if ligand_1_atom_id.atom == "CA":
                if match_atom(ligand_1_atom, ligand_2_atom, ignore_chain=True):
                    alignable_cas.append(
                        (
                            gemmi.Position(
                                ligand_1_atom.x,
                                ligand_1_atom.y,
                                ligand_1_atom.z,
                            ),
                            gemmi.Position(
                                ligand_2_atom.x,
                                ligand_2_atom.y,
                                ligand_2_atom.z,
                            ),
                        )
                    )

    if len(alignable_cas) > min_alignable_atoms:
        sup = gemmi.superpose_positions(
            [alignable_ca[0] for alignable_ca in alignable_cas],
            [alignable_ca[1] for alignable_ca in alignable_cas],
        )

        rmsd = sup.rmsd
        if rmsd < max_alignable_rmsd:
            return True
        else:
            return False
    else:
        return False


def get_alignability(
    ligand_neighbourhoods: LigandNeighbourhoods,
    system_data: SystemData,
):

    # Get structures
    structures = {}
    for dataset in system_data.dataset:
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        structures[dataset.dtag] = structure

    # Get connectivity matrix
    connectivity = []
    for (ligand_1_id, ligand_1_neighbourhood) in zip(
        ligand_neighbourhoods.ligand_ids,
        ligand_neighbourhoods.ligand_neighbourhoods.values(),
    ):
        connectivities = []
        for (ligand_2_id, ligand_2_neighbourhood,) in zip(
            ligand_neighbourhoods.ligand_ids,
            ligand_neighbourhoods.ligand_neighbourhoods.values(),
        ):
            # See if atoms match
            ca_match = match_cas(
                ligand_1_neighbourhood, ligand_2_neighbourhood
            )

            if ca_match:
                connectivities.append(1)
            else:
                connectivities.append(0)

        connectivity.append(connectivities)

    logger.debug(connectivity)

    return np.array(connectivity)
