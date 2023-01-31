from math import ceil, floor
from pathlib import Path

import gemmi
import networkx as nx
import numpy as np
from loguru import logger

from xchemalign.data import (
    Block,
    Dataset,
    DatasetID,
    LigandBindingEvent,
    LigandID,
    LigandNeighbourhood,
    LigandNeighbourhoods,
    Sites,
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

    tlbf = cell.fractionalize(tlb)
    tubf = cell.fractionalize(tub)

    tlbg = (
        floor(xmap.dx * tlbf.x),
        floor(xmap.dy * tlbf.y),
        floor(xmap.dz * tlbf.z),
    )
    tubg = (
        ceil(xmap.dx * tubf.x),
        ceil(xmap.dy * tubf.y),
        ceil(xmap.dz * tubf.z),
    )

    return tlbg, tubg


def get_subblocks(xrm, xr):
    xb = []
    xbt = []
    xit = xr[0]
    xbi = []
    for x, xm in zip(xr, xrm):
        if xm == 0:
            xb.append(xm)
            xbi.append(xit)
            xb = x

        else:
            xbt.append(xm)
    if len(xbt) > 0:
        xb.append(xbt)
        xbi.append(x)
    return xb, xbi


def get_blocks(rglb, rgub, xmap):

    cell = xmap.cell

    xr = np.arrange(rglb[0], rgub[0])
    yr = np.arrange(rglb[1], rgub[1])
    zr = np.arrange(rglb[2], rgub[2])

    xrm = np.mod(xr, xmap.dx)
    yrm = np.mod(yr, xmap.dy)
    zrm = np.mod(zr, xmap.dz)

    xb, xbi = get_subblocks(xrm)
    yb, ybi = get_subblocks(yrm)
    zb, zbi = get_subblocks(zrm)

    blocks = []
    for xsb, xsbi in zip(xb, xbi):
        for ysb, ysbi in zip(yb, ybi):
            for zsb, zsbi in zip(zb, zbi):
                # Transform jkl -> fractional -> orthogonal -> initial point
                transform_vec = (
                    xsbi * xmap.dz,
                    ysbi * xmap.dy,
                    zsbi * xmap.dx,
                )
                orth_arr = np.array(cell.orth.mat.tolist())
                frac_arr = np.diag([1 / cell.a, 1 / cell.b, 1 / cell.c])
                block = Block(
                    xi=xsb[0],
                    yi=ysb[0],
                    zi=zsb[0],
                    dx=len(xsb),
                    dy=len(ysb),
                    dz=len(zsb),
                    transform=Transform(
                        vec=transform_vec, mat=(orth_arr @ frac_arr).tolist()
                    ),
                )

                blocks.append(block)

    return blocks


def get_interpolation_range(
    neighbourhood: LigandNeighbourhood, transform, reference_xmap
):

    # Get the gemmi transform
    transform_gemmi = transform

    # Coord array
    coord_arr = get_coord_array(neighbourhood)

    # Find the bounds
    lb, ub = get_bounds(coord_arr)
    logger.debug(f"Neighbourhood bounds are: {lb} : {ub}")

    # Transform the bounds
    rlb, rub = transform_gemmi.apply(
        gemmi.Position(*lb)
    ), transform_gemmi.apply(gemmi.Position(*ub))
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
    new_xmap = gemmi.FloatGrid(
        reference_xmap.dx, reference_xmap.dy, reference_xmap.dz
    )
    new_xmap.set_unit_cell(reference_xmap.unit_cell)
    new_xmap.spacegroup = reference_xmap.spacegroup
    grid_np = np.array(new_xmap, copy=False)

    # Interpolate values
    for interpolation_range in interpolation_ranges:
        arr = np.zeros(
            interpolation_range.dx,
            interpolation_range.dy,
            interpolation_range.dz,
        )

        range_transform = gemmi.Transform()
        range_transform.vec.fromlist(interpolation_range.transform.vec)
        range_transform.mat.fromlist(interpolation_range.transform.mat)

        xmap.interpolate_values(
            arr, transform.combine(interpolation_range.transform)
        )

        # Assign to blocks
        rxi = interpolation_range.xi
        rxf = interpolation_range.xi + interpolation_range.dx
        ryi = interpolation_range.yi
        ryf = interpolation_range.yi + interpolation_range.dy
        rzi = interpolation_range.zi
        rzf = interpolation_range.zi + interpolation_range.dz
        grid_np[
            rxi:rxf,
            ryi:ryf,
            rzi:rzf,
        ] = arr

    return new_xmap


def _align_xmaps(
    system_data: SystemData,
    structures,
    sites: Sites,
    neighbourhoods: LigandNeighbourhoods,
    g,
    transforms: Transforms,
    site_transforms: SiteTransforms,
    _output_dir: Path,
):

    # Get the global reference
    reference_lid: LigandID = sites.sites[0].members[0]

    # Get that dataset
    referance_ds: Dataset = system_data.get_dataset(
        DatasetID(dtag=reference_lid.dtag)
    )
    reference_binding_site = referance_ds.ligand_binding_events[reference_lid]

    # Reference_xmap_path
    reference_xmap_path: Path = Path(reference_binding_site.xmap)

    # Load the site reference xmap
    reference_xmap = read_xmap(reference_xmap_path)

    for site_id, site in sites.iter():
        logger.debug(f"Aligning site: {site_id}")
        # site_reference_id = site.members[0]
        for subsite_id, subsite in site.iter():
            logger.debug(f"Aligning subsite: {subsite_id}")
            # Get the site reference
            # TODO: Make work
            subsite_reference_id = subsite.members[0]

            # For each ligand neighbourhood, find the xmap and transform
            for lid in subsite.members:
                # for lid, neighbourhood in zip(neighbourhoods.ligand_ids,
                # neighbourhoods.ligand_neighbourhoods):
                logger.debug(f"Aligning xmap: {lid}")

                # Get the ligand neighbourhood
                neighbourhood: LigandNeighbourhood = (
                    neighbourhoods.get_neighbourhood(lid)
                )

                # Get the ligand binding event
                lbe: LigandBindingEvent = system_data.get_dataset(
                    DatasetID(dtag=lid.dtag)
                ).ligand_binding_events[lid]

                # Get the xmap path
                xmap_path: Path = Path(lbe.xmap)
                logger.debug(f"Xmap path: {xmap_path}")

                # Get the xmap
                xmap = read_xmap(xmap_path)

                # Get the Transform to reference
                running_transform = gemmi.Transform()
                shortest_path = nx.shortest_path(g, lid, subsite_reference_id)

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
                        running_transform = transform.combine(
                            transform_to_gemmi(running_transform)
                        )

                # Get the subsite transform
                subsite_transform = transform_to_gemmi(
                    site_transforms.get_subsite_transform(site_id, subsite_id)
                )

                # Get the site transform
                site_transform = transform_to_gemmi(
                    site_transforms.get_site_transform(site_id)
                )

                # Running transform
                running_transform = site_transform.combine(
                    subsite_transform.combine(running_transform)
                )

                logger.debug(
                    f"Transform is: {gemmi_to_transform(running_transform)}"
                )

                # Define the interpolation range
                interpolation_range = get_interpolation_range(
                    neighbourhood, running_transform, reference_xmap
                )

                # Interpolate
                new_xmap = interpolate_range(
                    reference_xmap, xmap, interpolation_range, transform
                )

                # Output the xmap
                write_xmap(
                    new_xmap, _output_dir / "", neighbourhood, transform
                )
