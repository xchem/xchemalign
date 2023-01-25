from math import ceil, floor
from pathlib import Path

import gemmi
import numpy as np

from xchemalign.data import (
    Block,
    LigandNeighbourhood,
    LigandNeighbourhoods,
    Sites,
    SystemData,
    Transform,
    Transforms,
)


def transform_to_gemmi(transform: Transform):
    transform_gemmi = gemmi.Transform()
    transform_gemmi.vec.fromlist(transform.vec)
    transform_gemmi.mat.fromlist(transform.mat)

    return transform_gemmi


def read_xmap(path: Path):
    m = gemmi.read_ccp4_map(str(path), setup=True)
    return m.grid


def get_box(neighbourhood: LigandNeighbourhood, xmap, transform: Transform):

    transform_gemmi = transform_to_gemmi(transform)

    box = gemmi.FractionalBox()
    for atom in neighbourhood.atoms:
        box.extend(
            xmap.cell.fractionalize(
                transform_gemmi.apply(gemmi.Position(atom.x, atom.y, atom.z))
            )
        )

    for atom in neighbourhood.atoms:
        box.extend(
            xmap.cell.fractionalize(
                transform_gemmi.apply(gemmi.Position(atom.x, atom.y, atom.z))
            )
        )

    return box


def write_xmap(
    xmap, path: Path, neighbourhood: LigandNeighbourhood, transform: Transform
):

    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = xmap
    box = get_box(neighbourhood, xmap, transform)
    ccp4.set_extent(box)
    ccp4.update_ccp4_header()

    ccp4.write_ccp4_map(str(path))


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
    arr = np.array([rlb.x, rlb.y, rlb.z], [rub.x, rub.y, rub.z])

    lb, ub = get_bounds(arr)
    return lb, ub


def get_grid_bounds(tlb, tub, xmap):

    cell = xmap.cell

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
    neighbourhood: LigandNeighbourhood, transform: Transform, reference_xmap
):

    # Get the gemmi transform
    transform_gemmi = transform_to_gemmi(transform)

    # Coord array
    coord_arr = get_coord_array(neighbourhood)

    # Find the bounds
    lb, ub = get_bounds(coord_arr)

    # Transform the bounds
    rlb, rub = transform_gemmi.apply(
        gemmi.Position(*lb)
    ), transform_gemmi.apply(gemmi.Position(*ub))

    # Get the new bounds
    tlb, tub = get_transformed_bounds(rlb, rub)

    # Get grid bounds
    rglb, rgub = get_grid_bounds(tlb, tub, reference_xmap)

    # Get the blocks
    blocks: list[Block] = get_blocks(rglb, rgub, reference_xmap)

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
    sites: Sites,
    neighbourhoods: LigandNeighbourhoods,
    transforms: Transforms,
    _output_dir: Path,
):

    for site in sites.sites:
        # Get the site reference
        reference_ligand_id = site.reference

        # Get the refeence binding site
        rlbes = system_data.get_dataset(reference_ligand_id.dtag)
        reference_binding_site = rlbes.ligand_binding_events[
            reference_ligand_id
        ]

        # Reference_xmap_path
        reference_xmap_path: Path = Path(reference_binding_site.xmap)

        # Load the site reference xmap
        reference_xmap = read_xmap(reference_xmap_path)

        # For each ligand neighbourhood, find the xmap and transform
        for lid in site.ligand_ids:
            # for lid, neighbourhood in zip(neighbourhoods.ligand_ids,
            # neighbourhoods.ligand_neighbourhoods):

            # Get the ligand neighbourhood
            neighbourhood: LigandNeighbourhood = (
                neighbourhoods.get_neighbourhood(lid)
            )

            # Get the ligand binding event
            lbe = system_data.get_dataset(lid.dtag).ligand_binding_events[lid]

            # Get the xmap path
            xmap_path: Path = Path(lbe.xmap)

            # Get the xmap
            xmap = read_xmap(xmap_path)

            # Get the Transform to reference
            transform = transforms.get_transform((reference_ligand_id, lid))

            # Define the interpolation range
            interpolation_range = get_interpolation_range(
                neighbourhood, transform, reference_xmap
            )

            # Interpolate
            new_xmap = interpolate_range(
                reference_xmap, xmap, interpolation_range, transform
            )

            # Output the xmap
            write_xmap(new_xmap, _output_dir / "", neighbourhood, transform)

    # for dtag, dataset in zip(system_data.dataset_ids, system_data.datasets):
