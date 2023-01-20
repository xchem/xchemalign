import argparse
from pathlib import Path

# import gemmi
from loguru import logger

from xchemalign.data import (
    CanonicalSite,
    LigandID,
    LigandNeighbourhood,
    SiteObservation,
    SystemData,
    SystemSites,
    XtalFormSite,
)
from xchemalign.generate_aligned_structures import generate_aligned_structures
from xchemalign.get_alignability import get_alignability
from xchemalign.get_alignable_sites import get_alignable_sites
from xchemalign.get_canonical_sites import get_canonical_sites
from xchemalign.get_connected_components import get_connected_components
from xchemalign.get_ligand_neighbourhoods import get_ligand_neighbourhoods
from xchemalign.get_site_observations import get_site_observations
from xchemalign.get_xtal_form_sites import get_xtal_form_sites


def get_system_sites(
    system_sites_json_path: Path, data_json_path: Path, output_dir: Path
):
    logger.info(f"System sites json path is: {system_sites_json_path}")
    logger.info(f"Data json path is: {data_json_path}")

    # Load the input data
    initial_system_sites: SystemSites | None = None
    if system_sites_json_path.exists():
        initial_system_sites = SystemSites.parse_file(
            str(system_sites_json_path)
        )
    if initial_system_sites:
        logger.info(initial_system_sites)
    system_data: SystemData = SystemData.parse_file(str(data_json_path))

    # TODO: REMOVE
    # system_data.dataset = system_data.dataset[:20]

    if initial_system_sites:
        num_canonical_sites = len(initial_system_sites.canonical_site)
        num_xtal_form_sites = len(initial_system_sites.xtal_form_site)
        num_site_observations = len(initial_system_sites.site_observation)

        logger.debug(f"Got {num_canonical_sites} existinging canonical sites")
        logger.debug(f"Got {num_xtal_form_sites} existinging xtal form sites")
        logger.debug(
            f"Got {num_site_observations} existinging site observations"
        )

    num_input_datasets = len(system_data.dataset)

    logger.debug(f"Got {num_input_datasets} datasets")

    # Identify sites in the input data
    ligand_neighbourhoods: dict[
        LigandID, LigandNeighbourhood
    ] = get_ligand_neighbourhoods(system_data)

    num_neighbourhoods = len(ligand_neighbourhoods)
    logger.info(f"Found {num_neighbourhoods} ligand neighbourhoods")

    # Get alignability
    alignability_matrix = get_alignability(ligand_neighbourhoods, system_data)

    # logger.debug(alignability_matrix)
    logger.debug(alignability_matrix.shape)

    # Get connected components
    connected_components: list[list[LigandID]] = get_connected_components(
        alignability_matrix, ligand_neighbourhoods
    )
    logger.info(f"Found {len(connected_components)} connected components!")
    # logger.debug(connected_components)

    # Merge heavily connected components
    # connected_components = merge_components(connected_components)

    # Form sites
    sites = get_alignable_sites(connected_components, [])
    logger.info(f"Found {len(sites)} alignable sites!")
    # logger.debug(sites)

    # Generate aligned sites
    generate_aligned_structures(
        output_dir, ligand_neighbourhoods, system_data, sites
    )

    exit()

    # Get the canonical sites
    canonical_sites: dict[int, CanonicalSite] = get_canonical_sites(
        initial_system_sites, ligand_neighbourhoods
    )

    num_new_canon_sites = len(canonical_sites)
    logger.info(f"New number of canonical sites is {num_new_canon_sites}")
    logger.debug(
        [len(canon_site.atoms) for canon_site in canonical_sites.values()]
    )

    # Identify the xtalform sites
    xtal_form_sites: dict[int, XtalFormSite] = get_xtal_form_sites(
        initial_system_sites,
        ligand_neighbourhoods,
        canonical_sites,
    )

    new_num_xtal_form_sites = len(xtal_form_sites)
    logger.info(f"New number of xtal form sites is {new_num_xtal_form_sites}")
    # logger.debug([xtal_form_sites)

    # Construct the observed sites
    site_observations: dict[LigandID, SiteObservation] = get_site_observations(
        initial_system_sites,
        xtal_form_sites,
        ligand_neighbourhoods,
    )

    new_num_site_observations = len(site_observations)
    logger.info(
        f"New number of site observations is {new_num_site_observations}"
    )
    # logger.debug(site_observations)

    # Construct the new system sites
    system_sites: SystemSites = SystemSites(
        canonical_sites=canonical_sites,
        xtal_form_sites=xtal_form_sites,
        site_observations=site_observations,
    )

    # Output the new system_sites
    logger.info(f"Updating sites json at {system_sites_json_path}")
    with open(system_sites_json_path, "w") as f:
        f.write(system_sites.json())


# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    # Parse args
    logger.info("Parsing args...")

    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--system_sites_json_path",
        type=Path,
        help="",
    )
    parser.add_argument(
        "--data_json_path",
        type=Path,
        help="",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        help="",
    )

    args = parser.parse_args()

    # Run the program
    get_system_sites(
        args.system_sites_json_path, args.data_json_path, args.output_dir
    )
