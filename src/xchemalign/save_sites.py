from pathlib import Path

from xchemalign import constants
from xchemalign.data import CanonicalSites


def save_sites(sites: CanonicalSites, path: Path):
    with open(path / constants.SITES_FILE_NAME, "w") as f:
        f.write(sites.json())
