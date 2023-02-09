import gemmi
from loguru import logger

from xchemalign.data import (
    LigandNeighbourhood,
    ResidueID,
    SystemData,
    XtalForm,
)
from xchemalign.matching import match_atom


def get_structures(system_data: SystemData):
    structures = {}
    for dataset in system_data.datasets:
        structure: gemmi.Structure = gemmi.read_structure(dataset.pdb)
        structures[dataset.dtag] = structure

    return structures


def _get_transform():
    ...


def get_transform_from_residues(rs: list[ResidueID], srs, ssrs):
    # Transform from ssrs to srs
    acs = []
    for resid in rs:
        chain, num = resid.chain, resid.residue
        try:
            srsr = srs[0][chain][f"{num}"][0]
            ssrsr = ssrs[0][chain][f"{num}"][0]

            srsca = srsr["CA"][0]
            ssrsca = ssrsr["CA"][0]
            acs.append((srsca, ssrsca))
        except Exception as e:
            # print(f"{chain} : {num}: : {srsr.name} {ssrsr.name} {e}")
            print(f"{chain} : {num} : {e}")

            continue

    logger.debug(f"{len(acs)}")
    if len(acs) < 3:
        print("####### SRS")

        for model in srs:
            for c in model:
                for r in c:
                    print(f"{c.name}: {r.seqid.num}")

        print("####### SSRS")
        for model in ssrs:
            for c in model:
                for r in c:
                    print(f"{c.name}: {r.seqid.num}")

        print("####### RESIDS")
        for rid in rs:
            print(f"{rid.chain} {rid.residue}")

        raise Exception()

    sup = gemmi.superpose_positions(
        [x[0].pos for x in acs], [x[1].pos for x in acs]
    )

    return sup.transform


def get_transforms(
    reference_neighbourhood: LigandNeighbourhood,
    neighbourhood: LigandNeighbourhood,
):

    alignable_cas = {}
    for (
        ligand_1_atom_id,
        ligand_1_atom,
    ) in zip(reference_neighbourhood.atom_ids, reference_neighbourhood.atoms):
        for (
            ligand_2_atom_id,
            ligand_2_atom,
        ) in zip(neighbourhood.atom_ids, neighbourhood.atoms):
            if ligand_1_atom_id.atom == "CA":
                if match_atom(ligand_1_atom, ligand_2_atom, ignore_chain=True):
                    alignable_cas[ligand_1_atom_id] = (
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

    if len(alignable_cas) < 3:
        return gemmi.Transform(), []

    sup = gemmi.superpose_positions(
        [alignable_ca[0] for alignable_ca in alignable_cas.values()],
        [alignable_ca[1] for alignable_ca in alignable_cas.values()],
    )
    # logger.debug(f"Superposition: rmsd {sup.rmsd} n {len(alignable_cas)}")

    return sup.transform, [ligand_id for ligand_id in alignable_cas.keys()]


def generate_assembly(xtalform: XtalForm, structure):
    assembly = structure.clone()
    for model in assembly:
        for chain in model:
            del assembly[model.name][chain.name]

    for j, generator in enumerate(xtalform.generators):
        op = gemmi.Op(generator.triplet)
        chain_clone = structure[0][generator.chain].clone()
        for residue in chain_clone:
            for atom in residue:
                atom_frac = structure.cell.fractionalize(atom.pos)
                new_pos_frac = op.apply_to_xyz(
                    [atom_frac.x, atom_frac.y, atom_frac.z]
                )
                new_pos_orth = structure.cell.orthogonalize(
                    gemmi.Position(*new_pos_frac)
                )

                atom.pos = gemmi.Position(*new_pos_orth)
        chain_clone.name = f"{generator.chain}{j}"
        assembly[0].add_chain(chain_clone)

    num_chains = 0
    for model in assembly:
        for chain in model:
            num_chains += 1
    logger.debug(f"Generated {num_chains} assembly chains")

    return assembly


def remove_non_contact_chains(assembly, neighbourhood: LigandNeighbourhood):

    contact_chains_list = []
    for model in assembly:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    for atom_id, art_atom in zip(
                        neighbourhood.artefact_atom_ids,
                        neighbourhood.artefact_atoms,
                    ):
                        if (
                            atom.pos.dist(
                                gemmi.Position(
                                    art_atom.x, art_atom.y, art_atom.z
                                )
                            )
                            < 0.1
                        ):
                            contact_chains_list.append(chain.name)

    contact_chains = list(set(contact_chains_list))
    logger.debug(f"Num contact chains: {len(contact_chains)}")
    for model in assembly:
        for chain in model:
            if chain.name not in contact_chains:
                del model[chain.name]
