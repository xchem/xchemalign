import argparse
from pathlib import Path

# import gemmi
# import loguru
from xchemalign.data import (
    CanonicalSite,
    LigandID,
    LigandNeighbourhood,
    SiteObservation,
    SystemData,
    SystemSites,
    XtalFormSite,
)
from xchemalign.get_canonical_sites import get_canonical_sites
from xchemalign.get_ligand_neighbourhoods import get_ligand_neighbourhoods
from xchemalign.get_site_observations import get_site_observations
from xchemalign.get_xtal_form_sites import get_xtal_form_sites


def get_system_sites(system_sites_json_path: Path, data_json_path: Path):
    # Load the input data
    initial_system_sites: SystemSites = SystemSites.parse_file(
        str(system_sites_json_path)
    )
    system_data: SystemData = SystemData.parse_file(str(data_json_path))

    # Identify sites in the input data
    ligand_neighbourhoods: dict[
        LigandID, LigandNeighbourhood
    ] = get_ligand_neighbourhoods(system_data)

    # Get the canonical sites
    canonical_sites: dict[int, CanonicalSite] = get_canonical_sites(
        initial_system_sites, ligand_neighbourhoods
    )

    # Identify the xtalform sites
    xtal_form_sites: dict[int, XtalFormSite] = get_xtal_form_sites(
        initial_system_sites,
        ligand_neighbourhoods,
        canonical_sites,
    )

    # Construct the observed sites
    site_observations: dict[LigandID, SiteObservation] = get_site_observations(
        initial_system_sites,
        xtal_form_sites,
        ligand_neighbourhoods,
    )

    # Construct the new system sites
    system_sites: SystemSites = SystemSites(
        canonical_sites=canonical_sites,
        xtal_form_sites=xtal_form_sites,
        site_observations=site_observations,
    )

    # Output the new system_sites
    with open(system_sites_json_path, "w") as f:
        f.write(system_sites.json())


# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    # Parse args
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

    args = parser.parse_args()

    # Run the program
    get_system_sites(
        args.system_sites_json_path,
        args.data_json_path,
    )
