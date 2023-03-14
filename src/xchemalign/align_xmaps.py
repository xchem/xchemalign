from math import ceil, floor
from pathlib import Path

import gemmi
import networkx as nx
import numpy as np
from loguru import logger

from xchemalign.data import (
    Block,
    CanonicalSites,
    ConformerSites,
    Dataset,
    DatasetID,
    LigandBindingEvent,
    LigandID,
    LigandNeighbourhood,
    LigandNeighbourhoods,
    Output,
    SiteTransforms,
    SystemData,
    Transform,
    Transforms,
    gemmi_to_transform,
    read_xmap,
    transform_to_gemmi,
    write_xmap,
)


def get_coord_array(neighbourhood: LigandNeighbourhood):
    coords = []
    for a in neighbourhood.atoms:
        coords.append([a.x, a.y, a.z])
    for a in neighbourhood.artefact_atoms:
        coords.append([a.x, a.y, a.z])

    return np.array(coords)


def get_bounds(coord_array):
    lb = np.min(coord_array, axis=0)
    up = np.max(coord_array, axis=0)

    return lb, up


def get_transformed_bounds(rlb, rub):
    arr = np.array([[rlb.x, rlb.y, rlb.z], [rub.x, rub.y, rub.z]])

    lb, ub = get_bounds(arr)
    return lb, ub


def get_grid_bounds(tlb, tub, xmap):

    cell = xmap.unit_cell

    tlbf = cell.fractionalize(gemmi.Position(*tlb))
    tubf = cell.fractionalize(gemmi.Position(*tub))

    tlbg = (
        floor(xmap.nu * tlbf.x),
        floor(xmap.nv * tlbf.y),
        floor(xmap.nw * tlbf.z),
    )
    tubg = (
        ceil(xmap.nu * tubf.x),
        ceil(xmap.nv * tubf.y),
        ceil(xmap.nw * tubf.z),
    )

    return tlbg, tubg


def get_subblocks(xrm, xr):
    xb = []
    xbt = []
    xit = xr[0]
    xbi = []
    for x, xm in zip(xr, xrm):
        if (xm == 0) & (len(xbt) > 0):
            xb.append(xbt)
            xbi.append(xit)
            xit = x
            xbt = []

        else:
            xbt.append(x)
    if len(xbt) > 0:
        xb.append(xbt)
        xbi.append(xit)
    return xb, xbi


def get_blocks(rglb, rgub, xmap):

    cell = xmap.unit_cell

    xr = np.arange(rglb[0], rgub[0] + 1)
    yr = np.arange(rglb[1], rgub[1] + 1)
    zr = np.arange(rglb[2], rgub[2] + 1)
    logger.debug(f"X range: {xr[0]} : {xr[-1]}")

    xrm = np.mod(xr, xmap.nu)
    yrm = np.mod(yr, xmap.nv)
    zrm = np.mod(zr, xmap.nw)
    logger.debug(f"X mod range: {xrm[0]} : {xrm[-1]}")

    xb, xbi = get_subblocks(xrm, xr)
    yb, ybi = get_subblocks(yrm, yr)
    zb, zbi = get_subblocks(zrm, zr)
    logger.debug(f"X subblock 0 range: {xb[0][0]} : {xb[0][-1]}")
    logger.debug(f"X subblock initial: {xbi[0]}")

    blocks = []
    for xsb, xsbi in zip(xb, xbi):
        for ysb, ysbi in zip(yb, ybi):
            for zsb, zsbi in zip(zb, zbi):
                # Transform jkl -> fractional -> orthogonal -> initial point
                # Fractional pos of block 0,0,0
                transform_vec = xmap.unit_cell.orthogonalize(
                    gemmi.Fractional(
                        xsbi / xmap.nu,
                        ysbi / xmap.nv,
                        zsbi / xmap.nw,
                    )
                )
                # Transform mat: fractional to orth
                orth_arr = np.array(cell.orth.mat.tolist())
                # Transform mat: grid to frac
                # frac_arr = np.diag([1 / cell.a, 1 / cell.b, 1 / cell.c])
                frac_arr = np.diag([1 / xmap.nu, 1 / xmap.nv, 1 / xmap.nw])

                block = Block(
                    xi=xsb[0],
                    yi=ysb[0],
                    zi=zsb[0],
                    xmi=np.mod(xsbi, xmap.nu),
                    ymi=np.mod(ysbi, xmap.nv),
                    zmi=np.mod(zsbi, xmap.nw),
                    dx=len(xsb),
                    dy=len(ysb),
                    dz=len(zsb),
                    transform=Transform(
                        vec=[
                            transform_vec.x,
                            transform_vec.y,
                            transform_vec.z,
                        ],
                        mat=(orth_arr @ frac_arr).tolist(),
                    ),
                )

                blocks.append(block)

    return blocks


def get_interpolation_range(neighbourhood: LigandNeighbourhood, transform, reference_xmap):

    # Get the gemmi transform
    transform_gemmi = transform

    # Coord array
    coord_arr = get_coord_array(neighbourhood)

    # Find the bounds
    lb, ub = get_bounds(coord_arr)
    logger.debug(f"Neighbourhood bounds are: {lb} : {ub}")

    # Transform the bounds
    rlb, rub = transform_gemmi.apply(gemmi.Position(*lb)), transform_gemmi.apply(gemmi.Position(*ub))
    logger.debug(f"Reference bounds are: {lb} : {ub}")

    # Get the new bounds
    tlb, tub = get_transformed_bounds(rlb, rub)
    logger.debug(f"Transformed bounds are: {tlb} : {tub}")

    # Get grid bounds
    rglb, rgub = get_grid_bounds(tlb, tub, reference_xmap)
    logger.debug(f"Grid bounds are: {rglb} : {rgub}")

    # Get the blocks
    blocks: list[Block] = get_blocks(rglb, rgub, reference_xmap)
    logger.debug(f"Num blocks: {blocks}")
    for b in blocks:
        s = f"Block: {b.xi} {b.yi} {b.zi} {b.dx} {b.dy} {b.dz}"
        logger.debug(s)
        logger.debug(b.transform)

    return blocks


def interpolate_range(
    reference_xmap,
    xmap,
    interpolation_ranges: list[Block],
    transform,
):
    # Make a xmap on reference template
    new_xmap = gemmi.FloatGrid(reference_xmap.nu, reference_xmap.nv, reference_xmap.nw)
    new_xmap.set_unit_cell(reference_xmap.unit_cell)
    new_xmap.spacegroup = gemmi.find_spacegroup_by_name("P1")
    grid_np = np.array(new_xmap, copy=False)
    logger.debug(f"Xmap shape is: {grid_np.shape}")

    # Interpolate values
    for interpolation_range in interpolation_ranges:
        arr = np.zeros(
            (
                interpolation_range.dx,
                interpolation_range.dy,
                interpolation_range.dz,
            ),
            dtype=np.float32,
        )
        logger.debug(f"Sample grid shape: {arr.shape}")

        range_transform = gemmi.Transform()
        range_transform.vec.fromlist(interpolation_range.transform.vec)
        range_transform.mat.fromlist(interpolation_range.transform.mat)

        xmap.interpolate_values(
            arr,
            transform.combine(transform_to_gemmi(interpolation_range.transform)),
        )

        # Assign to blocks
        rxi = interpolation_range.xmi
        rxf = interpolation_range.xmi + interpolation_range.dx
        ryi = interpolation_range.ymi
        ryf = interpolation_range.ymi + interpolation_range.dy
        rzi = interpolation_range.zmi
        rzf = interpolation_range.zmi + interpolation_range.dz
        logger.debug(f"Block X Range: {rxi} : {rxf}")
        logger.debug(f"Block Y Range: {ryi} : {ryf}")
        logger.debug(f"Block Z Range: {rzi} : {rzf}")

        grid_np[
            rxi:rxf,
            ryi:ryf,
            rzi:rzf,
        ] = arr

    return new_xmap


def align_xmap(
    neighbourhoods: LigandNeighbourhoods,
    g,
    transforms: Transforms,
    site_transforms: SiteTransforms,
    reference_xmap,
    subsite_reference_id: LigandID,
    site_id: int,
    subsite_id: int,
    lid: LigandID,
    xmap,
    output_path: Path,
):
    # Get the ligand neighbourhood
    neighbourhood: LigandNeighbourhood = neighbourhoods.get_neighbourhood(lid)

    # Get the xmap

    # Get the Transform to reference
    running_transform = gemmi.Transform()
    shortest_path = nx.shortest_path(g, lid, subsite_reference_id)
    logger.debug(f"Shortest path: {shortest_path}")

    previous_ligand_id = lid
    for next_ligand_id in shortest_path:
        # Get the transform from previous frame to new one
        # Transform is 2 onto 1
        if next_ligand_id != previous_ligand_id:

            transform = transforms.get_transform(
                (
                    next_ligand_id,
                    previous_ligand_id,
                ),
            )
            running_transform = transform.combine(running_transform)
        # Apply the translation to the new frame
        previous_ligand_id = next_ligand_id

    # Get the subsite transform
    subsite_transform = transform_to_gemmi(site_transforms.get_conformer_site_transform(site_id, subsite_id))

    # Get the site transform
    site_transform = transform_to_gemmi(site_transforms.get_canonical_site_transform(site_id))

    # Running transform
    running_transform = site_transform.combine(subsite_transform.combine(running_transform))

    logger.debug(f"Transform is: {gemmi_to_transform(running_transform)}")

    # Define the interpolation range
    interpolation_range = get_interpolation_range(neighbourhood, running_transform, reference_xmap)

    # Interpolate
    new_xmap = interpolate_range(
        reference_xmap,
        xmap,
        interpolation_range,
        running_transform.inverse(),
    )

    # Output the xmap
    write_xmap(
        new_xmap,
        output_path,
        neighbourhood,
        running_transform,
    )


def read_xmap_from_mtz(mtz_path: Path):
    mtz = gemmi.read_mtz_file(str(mtz_path))
    try:
        grid = mtz.transform_f_phi_to_map("2FOFCWT", "PH2FOFCWT", sample_rate=4)
        return grid
    except Exception:
        logger.warning("Trying FWT PHWT")
    try:
        grid = mtz.transform_f_phi_to_map("FWT", "PHWT", sample_rate=4)
        return grid
    except Exception:
        logger.warning("Failed to find structure factors!")

    raise Exception(f"Couldn't open {mtz_path}!")


def _align_xmaps(
    system_data: SystemData,
    structures,
    canonical_sites: CanonicalSites,
    conformer_sites: ConformerSites,
    neighbourhoods: LigandNeighbourhoods,
    g,
    transforms: Transforms,
    site_transforms: SiteTransforms,
    output: Output,
):

    # Get the global reference
    # reference_lid: LigandID = canonical_sites.reference_site.reference_ligand_id

    # Get that dataset
    # referance_ds: Dataset = system_data.get_dataset(DatasetID(dtag=reference_lid.dtag))
    # reference_binding_site = referance_ds.ligand_binding_events[
    # reference_lid]
    # logger.debug(f"PDB: {referance_ds.pdb}")

    # Reference_xmap_path
    # reference_xmap_path: Path = Path(reference_binding_site.xmap)
    # reference_mtz_path = Path(referance_ds.mtz)

    # Load the site reference xmap
    # reference_xmap = read_xmap(reference_xmap_path)
    # reference_xmap = read_xmap_from_mtz(reference_mtz_path)

    #
    # xmaps_dir = _output_dir / "aligned_xmaps"
    # if not xmaps_dir.exists():
    #     os.mkdir(xmaps_dir)

    for canonical_site_id, canonical_site in canonical_sites.iter():
        logger.debug(f"Aligning site: {canonical_site_id}")
        # site_reference_id = site.members[0]

        #
        # site_xmaps_dir = xmaps_dir / f"{site_id}"
        # if not site_xmaps_dir.exists():
        #     os.mkdir(site_xmaps_dir)
        reference_lid: LigandID = canonical_site.reference_ligand_id
        referance_ds: Dataset = system_data.get_dataset(DatasetID(dtag=reference_lid.dtag))
        logger.debug(f"PDB: {referance_ds.pdb}")
        reference_mtz_path = Path(referance_ds.mtz)
        reference_xmap = read_xmap_from_mtz(reference_mtz_path)

        for conformer_site_id, conformer_site in conformer_sites.iter():
            if conformer_site_id not in canonical_site.subsites:
                continue
            logger.debug(f"Aligning subsite: {conformer_site_id}")
            # Get the site reference
            # TODO: Make work
            conformer_site_reference_id = conformer_site.reference_ligand_id

            # subsite_xmaps_dir = site_xmaps_dir / f"{subsite_id}"
            # if not subsite_xmaps_dir.exists():
            #     os.mkdir(subsite_xmaps_dir)

            # For each ligand neighbourhood, find the xmap and transform
            for lid in conformer_site.members:
                # for lid, neighbourhood in zip(neighbourhoods.ligand_ids,
                # neighbourhoods.ligand_neighbourhoods):
                logger.debug(f"Aligning xmap: {lid}")

                dtag, chain, residue = (
                    lid.dtag,
                    lid.chain,
                    lid.residue,
                )

                # Get the ligand binding event
                dataset = system_data.get_dataset(DatasetID(dtag=lid.dtag))
                lbe: LigandBindingEvent = dataset.ligand_binding_events[lid]

                # Get the xmap path

                if lbe.xmap != "None":
                    xmap_path: Path = Path(lbe.xmap)
                    logger.debug(f"Xmap path: {xmap_path}")

                    xmap = read_xmap(xmap_path)
                else:
                    mtz_path: Path = Path(dataset.mtz)
                    logger.debug(f"Mtz path is: {mtz_path}")

                    xmap = read_xmap_from_mtz(mtz_path)

                aop = Path(output.source_dir) / output.aligned_dir

                # Align the event map
                output_path = aop / output.dataset_output[dtag][chain][residue].aligned_event_maps[canonical_site_id]
                align_xmap(
                    neighbourhoods,
                    g,
                    transforms,
                    site_transforms,
                    reference_xmap,
                    conformer_site_reference_id,
                    canonical_site_id,
                    conformer_site_id,
                    lid,
                    xmap,
                    output_path,
                )
                # Align the refined map
                output_path = aop / output.dataset_output[dtag][chain][residue].aligned_xmaps[canonical_site_id]

                if dataset.mtz:

                    xmap = read_xmap_from_mtz(Path(dataset.mtz))

                    align_xmap(
                        neighbourhoods,
                        g,
                        transforms,
                        site_transforms,
                        reference_xmap,
                        conformer_site_reference_id,
                        canonical_site_id,
                        conformer_site_id,
                        lid,
                        xmap,
                        output_path,
                    )
